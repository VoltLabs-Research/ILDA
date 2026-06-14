#include <volt/ilda_crf.h>

#include <volt/core/frame_adapter.h>
#include <volt/analysis/cutoff_neighbor_finder.h>
#include <volt/math/affine_decomposition.h>
#include <volt/math/affine_transformation.h>

#include <cmath>
#include <set>
#include <vector>

#include <spdlog/spdlog.h>

namespace Volt {

namespace {

inline Matrix3 outer(const Vector3& u, const Vector3& w) {
    return Matrix3(
        u.x() * w.x(), u.x() * w.y(), u.x() * w.z(),
        u.y() * w.x(), u.y() * w.y(), u.y() * w.z(),
        u.z() * w.x(), u.z() * w.y(), u.z() * w.z()
    );
}

// PARITY: returns quaternion in Volt {x,y,z,w} == ILDA (qx,qy,qz,qw) order
inline Quaternion rotMatToQuat(const Matrix3& R) {
    const double qw = std::sqrt(1.0 + R(0, 0) + R(1, 1) + R(2, 2)) / 2.0;
    const double qx = (R(2, 1) - R(1, 2)) / (4.0 * qw);
    const double qy = (R(0, 2) - R(2, 0)) / (4.0 * qw);
    const double qz = (R(1, 0) - R(0, 1)) / (4.0 * qw);
    return Quaternion(qx, qy, qz, qw);
}

}

// PARITY: rebuild V from stretches so trace(V)==sum(S) matches Python eigendecomp
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
    const double a = scalePTMVec(type, grain.a);

    query.findNeighbors(particleIndex, &grain.orientation);

    Matrix3 W{Matrix3::Zero()};
    Matrix3 V{Matrix3::Zero()};
    for(int j = 0; j < query.count(); ++j) {
        const IldaPtmNeighborQuery::Neighbor neigh = query[j];
        const Vector3 idealA = neigh.idealVector * a;
        // PARITY: neigh.delta is already PBC-wrapped; pbc_vec is a no-op, not re-applied
        W = W + outer(neigh.delta, idealA);
        V = V + outer(idealA, idealA);
    }

    return W * V.inverse();
}

// PARITY: sums accumulated in ascending atom-index order (float-summation parity).
void estimateF(
    const LammpsParser::Frame& frame,
    const CisResult& cis,
    const std::vector<int>& grain,
    const std::vector<int>& structureType,
    const std::vector<std::uint8_t>& selection,
    double distF,
    IldaPtmNeighborQuery& query,
    GrainRef& gA,
    GrainRef& gB
) {
    const std::size_t n = static_cast<std::size_t>(frame.natoms);

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions) {
        spdlog::error("estimateF: failed to create position property");
        gA.F = Matrix3::Identity();
        gB.F = Matrix3::Identity();
        return;
    }
    CutoffNeighborFinder finder;
    if(!finder.prepare(distF, positions.get(), frame.simulationCell)) {
        spdlog::error("estimateF: CutoffNeighborFinder::prepare failed");
        gA.F = Matrix3::Identity();
        gB.F = Matrix3::Identity();
        return;
    }

    // PARITY: std::set keeps neighbors sorted ascending like np.unique
    std::set<std::size_t> neighborSet;
    for(std::size_t i = 0; i < n; ++i) {
        if(!cis.cis[i]) {
            continue;
        }
        for(CutoffNeighborFinder::Query q(finder, i); !q.atEnd(); q.next()) {
            neighborSet.insert(q.current());
        }
    }

    Matrix3 FsumA{Matrix3::Zero()};
    Matrix3 FsumB{Matrix3::Zero()};
    long NA = 0;
    long NB = 0;

    spdlog::info("Estimating FA and FB...");
    for(const std::size_t nbri : neighborSet) {
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

        // PARITY: redundant grain check kept (phase filter above already guarantees it)
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
        spdlog::warn("estimateF: no atoms found to estimate FA");
        gA.F = Matrix3::Identity();
    }
    if(NB > 0) {
        gB.F = FsumB * (1.0 / static_cast<double>(NB));
    } else {
        spdlog::warn("estimateF: no atoms found to estimate FB");
        gB.F = Matrix3::Identity();
    }

    spdlog::info("Used {} and {} atoms to estimate FA and FB, respectively.", NA, NB);

    Matrix3 VA, RA, VB, RB;
    polarDecomp(gA.F, VA, RA);
    polarDecomp(gB.F, VB, RB);
    spdlog::info("Estimated coherent reference frame (with polar decomposition)");
    spdlog::info("FA trace(VA)={:.6f}", VA(0, 0) + VA(1, 1) + VA(2, 2));
    spdlog::info("FB trace(VB)={:.6f}", VB(0, 0) + VB(1, 1) + VB(2, 2));
}

void computeF(
    const Vector3& xA, const Vector3& yA,
    const Vector3& xB, const Vector3& yB,
    const Matrix3& EcohA, const Matrix3& EcohB,
    double coveraA, double coveraB,
    GrainRef& gA, GrainRef& gB
) {
    if(std::abs(xA.dot(yA)) > 1e-6) {
        spdlog::error("computeF: xA and yA are not orthogonal!");
    }
    if(std::abs(xB.dot(yB)) > 1e-6) {
        spdlog::error("computeF: xB and yB are not orthogonal!");
    }

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
