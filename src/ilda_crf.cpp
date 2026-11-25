#include <volt/ilda_crf.h>

#include <volt/core/frame_adapter.h>
#include <volt/core/progress.h>
#include <volt/analysis/cutoff_neighbor_finder.h>
#include <volt/ilda_log.h>
#include <volt/math/affine_decomposition.h>
#include <volt/math/affine_transformation.h>

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

namespace Volt {

using IldaCrf::outer;
using IldaCrf::rotMatToQuat;

void polarDecomp(const Matrix3& F, Matrix3& V, Matrix3& R) {
    const AffineTransformation tm(
        F(0, 0), F(0, 1), F(0, 2),
        F(1, 0), F(1, 1), F(1, 2),
        F(2, 0), F(2, 1), F(2, 2)
    );
    const AffineDecomposition ad{tm};

    const Matrix3 Rrot = Matrix3::rotation(ad.rotation);
    const Matrix3 Qs = Matrix3::rotation(ad.scaling.Q);
    const Matrix3 M = Rrot * Qs;
    const Vector3& S = ad.scaling.S;
    const Matrix3 D(
        S.x(), 0.0, 0.0,
        0.0, S.y(), 0.0,
        0.0, 0.0, S.z()
    );

    V = M * D * M.transposed();
    R = V.inverse() * F;
}

Matrix3 totalDeformationGradPTM(
    std::size_t particleIndex,
    StructureType type,
    const GrainRef& grain,
    IldaPtmNeighborQuery& query
) {
    const double a = query.scalePTMVec(type, grain.a);

    query.findNeighbors(particleIndex, &grain.orientation);

    Matrix3 W{Matrix3::Zero()};
    Matrix3 V{Matrix3::Zero()};
    for(int j = 0; j < query.count(); ++j) {
        const IldaPtmNeighborQuery::Neighbor neigh = query[j];
        const Vector3 idealA = neigh.idealVector * a;
        W = W + outer(neigh.delta, idealA);
        V = V + outer(idealA, idealA);
    }

    return W * V.inverse();
}

void estimateF(
    const LammpsParser::Frame& frame,
    const CisResult& cis,
    const std::vector<int>& grain,
    const std::vector<int>& structureType,
    const std::vector<std::uint8_t>& selection,
    double distF,
    IldaPtmNeighborQuery& query,
    GrainRef& gA,
    GrainRef& gB,
    IldaLog& log
) {
    const std::size_t n = static_cast<std::size_t>(frame.natoms);
    const std::size_t cisAtoms = static_cast<std::size_t>(
        std::count(cis.cis.begin(), cis.cis.end(), static_cast<std::uint8_t>(1)));

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions) {
        log.dump(IldaIssue::InterfaceSkinPositions);
        gA.F = Matrix3::Identity();
        gB.F = Matrix3::Identity();
        return;
    }
    ProgressSteps steps({
        {"Binning the atoms", 1.0},
        {"Collecting atoms near CIS sites", 1.0},
        {"Fitting deformation gradients", 2.0},
    });
    CutoffNeighborFinder finder;
    if(!finder.prepare(distF, positions.get(), frame.simulationCell)) {
        log.dump(IldaIssue::InterfaceSkinNeighborLists);
        gA.F = Matrix3::Identity();
        gB.F = Matrix3::Identity();
        return;
    }
    steps.next();

    std::set<std::size_t> neighborSet;
    {
        ProgressReporter gatherProgress("  Collecting atoms near CIS sites", cisAtoms);
        std::size_t visitedCis = 0;
        for(std::size_t i = 0; i < n; ++i) {
            if(!cis.cis[i]) {
                continue;
            }
            for(CutoffNeighborFinder::Query q(finder, i); !q.atEnd(); q.next()) {
                neighborSet.insert(q.current());
            }
            gatherProgress.update(++visitedCis);
        }
    }
    log.dump(InterfaceSkin{neighborSet.size()});

    steps.next();

    Matrix3 FsumA{Matrix3::Zero()};
    Matrix3 FsumB{Matrix3::Zero()};
    long NA = 0;
    long NB = 0;

    ProgressReporter fitProgress("  Fitting deformation gradients", neighborSet.size());
    std::size_t fitted = 0;
    for(const std::size_t nbri : neighborSet) {
        fitProgress.update(fitted++);
        if(nbri >= n) {
            continue;
        }
        if(cis.cis[nbri]) {
            continue;
        }
        const bool inA = (structureType[nbri] == static_cast<int>(gA.structureType)) && (grain[nbri] == gA.grainId);
        const bool inB = (structureType[nbri] == static_cast<int>(gB.structureType)) && (grain[nbri] == gB.grainId);
        if(!(inA || inB)) {
            continue;
        }
        if(!selection[nbri]) {
            continue;
        }

        const StructureType atomType = static_cast<StructureType>(structureType[nbri]);
        if(grain[nbri] == gA.grainId) {
            const Matrix3 Fi = totalDeformationGradPTM(nbri, atomType, gA, query);
            FsumA = FsumA + Fi;
            ++NA;
        } else if(grain[nbri] == gB.grainId) {
            const Matrix3 Fi = totalDeformationGradPTM(nbri, atomType, gB, query);
            FsumB = FsumB + Fi;
            ++NB;
        }
    }

    if(NA > 0) {
        gA.F = FsumA * (1.0 / static_cast<double>(NA));
    } else {
        gA.F = Matrix3::Identity();
    }
    if(NB > 0) {
        gB.F = FsumB * (1.0 / static_cast<double>(NB));
    } else {
        gB.F = Matrix3::Identity();
    }
    log.dump(DeformationGradientFit{NA, NB, gA, gB});
}

void computeF(
    const Vector3& xA, const Vector3& yA,
    const Vector3& xB, const Vector3& yB,
    const Matrix3& EcohA, const Matrix3& EcohB,
    double coveraA, double coveraB,
    GrainRef& gA, GrainRef& gB
) {
    Vector3 zA = xA.cross(yA);
    Vector3 zB = xB.cross(yB);

    const Vector3 xAn = xA.normalized();
    const Vector3 yAn = yA.normalized();
    const Vector3 zAn = zA.normalized();

    const Vector3 xBn = xB.normalized();
    const Vector3 yBn = yB.normalized();
    const Vector3 zBn = zB.normalized();

    const Matrix3 RA(
        xAn.x(), xAn.y(), xAn.z(),
        yAn.x(), yAn.y(), yAn.z(),
        zAn.x(), zAn.y(), zAn.z()
    );
    const Matrix3 RB(
        xBn.x(), xBn.y(), xBn.z(),
        yBn.x(), yBn.y(), yBn.z(),
        zBn.x(), zBn.y(), zBn.z()
    );

    const double idealCoverA = std::sqrt(8.0 / 3.0);
    Matrix3 ZA{Matrix3::Identity()};
    if(coveraA > 0.0) {
        ZA(2, 2) = coveraA / idealCoverA;
    }
    Matrix3 ZB{Matrix3::Identity()};
    if(coveraB > 0.0) {
        ZB(2, 2) = coveraB / idealCoverA;
    }

    const Matrix3 I{Matrix3::Identity()};
    gA.F = (I + EcohA) * (RA * ZA);
    gB.F = (I + EcohB) * (RB * ZB);

    gA.orientation = rotMatToQuat(RA);
    gB.orientation = rotMatToQuat(RB);
}

}
