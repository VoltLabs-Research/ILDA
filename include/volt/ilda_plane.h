#pragma once

#include <volt/ilda_types.h>
#include <volt/core/simulation_cell.h>
#include <volt/math/vector3.h>

#include <utility>

namespace Volt {

std::pair<Vector3, double> findTerracePlane2(
    IldaMesh& mesh,
    IldaMesh::Vertex* startVertex,
    IldaMesh::Vertex* endVertex,
    const SimulationCell& cell,
    double btol
);

}
