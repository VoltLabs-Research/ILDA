#pragma once

#include <volt/ilda_types.h>
#include <volt/ilda_ptm_neighbor_query.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/simulation_cell.h>
#include <volt/math/matrix3.h>
#include <volt/math/quaternion.h>
#include <volt/math/vector3.h>

#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace Volt {

inline constexpr int ILDA_NTRY_MAX = 200;
inline constexpr double ILDA_CIRCUIT_TOL = 1e-6;

struct CircuitContext {
    const LammpsParser::Frame& frame;
    const SimulationCell& cell;
    const std::vector<Point3>& positions;
    const std::vector<int>& grain;
    const std::vector<int>& structureType;
    const std::vector<std::uint8_t>& selection;
    const CisResult& cis;
    const GrainRef& gA;
    const GrainRef& gB;
    Vector3 n{0.0, 0.0, 0.0};
    IldaPtmNeighborQuery& query;
};

std::optional<CircuitResult> burgersCircuit(
    int startIndex,
    int endIndex,
    const CircuitContext& ctx
);

std::optional<std::pair<Matrix3, Quaternion>> getHcpFccTrans(
    const IldaPtmNeighborQuery::Neighbor& fccNeigh,
    int fccIndex,
    IldaPtmNeighborQuery& fccQuery,
    int hcpIndex,
    const CircuitContext& ctx
);

}
