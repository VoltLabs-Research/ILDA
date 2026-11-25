#pragma once

#include <volt/ilda_types.h>
#include <volt/ilda_ptm_neighbor_query.h>
#include <volt/core/lammps_parser.h>
#include <volt/math/matrix3.h>
#include <volt/math/quaternion.h>
#include <volt/math/vector3.h>

#include <cmath>
#include <cstddef>
#include <vector>

namespace Volt {

class IldaLog;

namespace IldaCrf {

inline Matrix3 outer(const Vector3& u, const Vector3& w) {
    return Matrix3(
        u.x() * w.x(), u.x() * w.y(), u.x() * w.z(),
        u.y() * w.x(), u.y() * w.y(), u.y() * w.z(),
        u.z() * w.x(), u.z() * w.y(), u.z() * w.z()
    );
}

inline Quaternion rotMatToQuat(const Matrix3& R) {
    const double qw = std::sqrt(1.0 + R(0, 0) + R(1, 1) + R(2, 2)) / 2.0;
    const double qx = (R(2, 1) - R(1, 2)) / (4.0 * qw);
    const double qy = (R(0, 2) - R(2, 0)) / (4.0 * qw);
    const double qz = (R(1, 0) - R(0, 1)) / (4.0 * qw);
    return Quaternion(qx, qy, qz, qw);
}

}

void polarDecomp(const Matrix3& F, Matrix3& V, Matrix3& R);

Matrix3 totalDeformationGradPTM(
    std::size_t particleIndex,
    StructureType type,
    const GrainRef& grain,
    IldaPtmNeighborQuery& query
);

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
);

void computeF(
    const Vector3& xA, const Vector3& yA,
    const Vector3& xB, const Vector3& yB,
    const Matrix3& EcohA, const Matrix3& EcohB,
    double coveraA, double coveraB,
    GrainRef& gA, GrainRef& gB
);

}
