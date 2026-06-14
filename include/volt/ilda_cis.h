#pragma once

#include <volt/ilda_types.h>
#include <volt/ilda_ptm_neighbor_query.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>
#include <volt/core/simulation_cell.h>

#include <cstdint>
#include <vector>

namespace Volt {

CisResult findCIS(
    const LammpsParser::Frame& frame,
    const std::vector<int>& grain,
    const std::vector<int>& structureType,
    const GrainRef& gA,
    const GrainRef& gB,
    const std::vector<std::uint8_t>& selection,
    IldaPtmNeighborQuery& query,
    double cisTol
);

}
