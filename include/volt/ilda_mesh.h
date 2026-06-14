#pragma once

#include <volt/ilda_types.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>
#include <volt/core/simulation_cell.h>

#include <cstdint>
#include <vector>

namespace Volt {

inline constexpr int ILDA_PRUNE_MAX = 10000;

bool buildCisSurfaceMesh(
    const SimulationCell& cell,
    const LammpsParser::Frame& frame,
    const std::vector<std::uint8_t>& cis,
    const std::vector<int>& grain,
    const std::vector<int>& structureType,
    int grainB,
    int g2s,
    double Rsphere,
    IldaMesh& outMesh
);

}
