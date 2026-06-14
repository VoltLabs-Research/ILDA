#pragma once

#include <volt/ilda_types.h>
#include <volt/ilda_ptm_neighbor_query.h>
#include <volt/core/lammps_parser.h>
#include <volt/math/matrix3.h>
#include <volt/math/quaternion.h>
#include <volt/math/vector3.h>

#include <cstddef>
#include <vector>

namespace Volt {

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
    GrainRef& gB
);

void computeF(
    const Vector3& xA, const Vector3& yA,
    const Vector3& xB, const Vector3& yB,
    const Matrix3& EcohA, const Matrix3& EcohB,
    double coveraA, double coveraB,
    GrainRef& gA, GrainRef& gB
);

}
