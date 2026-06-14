#pragma once

#include <volt/helpers/half_edge_mesh.h>
#include <volt/math/vector3.h>
#include <volt/math/matrix3.h>
#include <volt/math/quaternion.h>
#include <volt/math/point3.h>
#include <volt/structures/crystal_structure_types.h>

#include <array>
#include <cstdint>
#include <vector>

namespace Volt {

struct IldaVertex {
    int atomIndex = -1;
    std::uint8_t isCIS = 0;
    Vector3 color{0.0, 0.0, 0.0};
};

struct IldaEdge {
    Vector3 burgers{0.0, 0.0, 0.0};
    double bnorm = 0.0;
    Vector3 plane{0.0, 0.0, 0.0};
    double step = 0.0;
    double stepIdealA = 0.0;
    double stepIdealB = 0.0;
    int bid = -1;
    int midedgeAtom = 0;
};

struct IldaFace {
    int region = 0;
};

using IldaMesh = Volt::HalfEdgeMesh<IldaEdge, IldaFace, IldaVertex>;

struct GrainRef {
    int grainId = 0;
    StructureType structureType = OTHER;
    Quaternion orientation{0.0, 0.0, 0.0, 1.0};
    Quaternion tableOrientation{0.0, 0.0, 0.0, 1.0};
    Matrix3 F{Matrix3::Identity()};
    double a = 1.0;
    double covera = 0.0;
    int grainIndex = 0;
};

struct IldaSegment {
    int atom1 = -1;
    int atom2 = -1;
    Point3 pos1{Point3::Origin()};
    Point3 pos2{Point3::Origin()};
    Vector3 burgers{0.0, 0.0, 0.0};
    double bnorm = 0.0;
    double step = 0.0;
    double stepIdealA = 0.0;
    double stepIdealB = 0.0;
    Vector3 plane{0.0, 0.0, 0.0};
    Vector3 color{0.0, 0.0, 0.0};
    int bid = -1;
};

struct IldaVariant {
    Vector3 bAvg{0.0, 0.0, 0.0};
    double bAvgNorm = 0.0;
    double bStdX = 0.0;
    double bStdY = 0.0;
    double bStdZ = 0.0;
    double length = 0.0;
    double hIdealA = 0.0;
    double hIdealB = 0.0;
    double hAvg = 0.0;
};

struct DisconnectionMode {
    int id = 0;
    std::vector<IldaVariant> variants;
};

struct CisResult {
    std::vector<std::uint8_t> cis;
    std::vector<std::array<int, 2>> cisNbrs;
    std::vector<std::array<Vector3, 2>> cisIdealVecs;
    std::vector<std::array<Vector3, 2>> cisSpatialVecs;
    std::vector<Vector3> color;
};

struct CircuitResult {
    Vector3 latticeSum{0.0, 0.0, 0.0};
    double hA = 0.0;
    double hB = 0.0;
};

}
