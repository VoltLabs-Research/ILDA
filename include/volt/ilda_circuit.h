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

class IldaLog;

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
    const CircuitContext& ctx,
    IldaLog& log
);

// The depth-first search for the path a Burgers circuit takes through one grain, and its state.
namespace IldaCircuit {

struct Neigh {
    double dist = 0.0;
    int index = -1;
    Vector3 delta{0.0, 0.0, 0.0};
    Vector3 a1{0.0, 0.0, 0.0};
    Matrix3 Ft2new{Matrix3::Identity()};
    Quaternion refOri2new{0.0, 0.0, 0.0, 0.0};
};

struct Frame {
    int atomIndex = -1;
    Vector3 spatialSum{0.0, 0.0, 0.0};
    Vector3 latticeSum{0.0, 0.0, 0.0};
    Matrix3 Ft2{Matrix3::Identity()};
    Quaternion refOri2{0.0, 0.0, 0.0, 0.0};
    int ntry = 0;
    bool entered = false;

    Matrix3 Ft{Matrix3::Identity()};
    int grainId = 0;
    StructureType gs = OTHER;
    Quaternion refOri{0.0, 0.0, 0.0, 1.0};
    double ptmScale = 1.0;

    std::vector<Neigh> neighs;
    std::size_t nextNeighbor = 0;
};

enum class Ret { DEADEND, ABORT, SUCCESS };

struct PathResult {
    Vector3 spatialSum{0.0, 0.0, 0.0};
    Vector3 latticeSum{0.0, 0.0, 0.0};
};

bool computeGrainProps(Frame& f, const CircuitContext& ctx);

void buildNeighs(Frame& f, int endIndex, const Vector3& target,
                 const std::vector<std::uint8_t>& visited,
                 const CircuitContext& ctx);

std::optional<PathResult> recursiveCircuitSearch(
    int startNbr, int startIndex, int endIndex, const Vector3& target,
    const CircuitContext& ctx
);

}

std::optional<std::pair<Matrix3, Quaternion>> getHcpFccTrans(
    const IldaPtmNeighborQuery::Neighbor& fccNeigh,
    int fccIndex,
    IldaPtmNeighborQuery& fccQuery,
    int hcpIndex,
    const CircuitContext& ctx
);

}
