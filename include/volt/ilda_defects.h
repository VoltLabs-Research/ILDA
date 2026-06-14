#pragma once

#include <volt/ilda_types.h>
#include <volt/ilda_circuit.h>
#include <volt/core/simulation_cell.h>
#include <volt/math/vector3.h>

#include <vector>

namespace Volt {

extern const Vector3 ILDA_COLORLIST[20];
inline constexpr int ILDA_NCOLORS = 20;

void findDefects(
    IldaMesh& mesh,
    const CircuitContext& ctx,
    double btol,
    double htol,
    double angtol,          // radians
    double bmagmax,
    std::vector<IldaSegment>& outSegments,
    std::vector<DisconnectionMode>& outModes,
    double& outTotalLength
);

}
