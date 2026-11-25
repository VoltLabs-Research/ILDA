#pragma once

#include <volt/helpers/half_edge_mesh.h>
#include <volt/math/vector3.h>
#include <volt/math/matrix3.h>
#include <volt/math/quaternion.h>
#include <volt/math/point3.h>
#include <volt/structures/crystal_structure_types.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
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

// The grain labels of the grain-segmentation atoms table, matched to the frame's atoms by particle identifier.
struct GrainLabels {
    std::size_t aligned = 0;
    std::size_t atoms = 0;
    std::string path;
};

struct GrainTable {
    std::size_t grains = 0;
    std::string path;
};

// The grains of the fiducial atoms A and B.
struct ResolvedGrains {
    const GrainRef& grainA;
    const GrainRef& grainB;
};

// The x and y axes of grains A and B the options give the coherent reference frame by; each pair must be orthogonal.
struct ReferenceFrameAxes {
    const Vector3& xA;
    const Vector3& yA;
    const Vector3& xB;
    const Vector3& yB;
};

// The coherent reference frame of grains A and B as the options give it, instead of estimated.
struct UserReferenceFrame {
    const GrainRef& grainA;
    const GrainRef& grainB;
};

// The grain A and grain B atoms next to a non-crystalline atom, among which the CIS atoms are searched.
struct InterfaceCandidates {
    std::size_t candidates = 0;
    std::size_t grainAAtoms = 0;
    std::size_t grainBAtoms = 0;
};

// Whether the first atom of a grain gave the reference lattice vectors the cis_tol extension matches against.
struct ReferenceLatticeVectors {
    int grainId = 0;
    bool found = false;
    double cisTol = 0.0;
};

// The atoms within distF of a CIS atom, whose deformation gradients give the coherent reference frame.
struct InterfaceSkin {
    std::size_t atoms = 0;
};

// The deformation gradients of grains A and B, averaged over this many atoms of each.
struct DeformationGradientFit {
    long atomsA = 0;
    long atomsB = 0;
    const GrainRef& grainA;
    const GrainRef& grainB;
};

// A Burgers circuit that could not be traced between two atoms, and why.
struct CircuitFailure {
    enum class Reason { InvalidIndices, NotCoincidentSites, SameAtom, NotClosed };

    Reason reason;
    int startIndex = -1;
    int endIndex = -1;
    // The spatial sum of a circuit that did not close.
    Vector3 spatialSum{0.0, 0.0, 0.0};
};

// The one Burgers circuit the options ask for, between two atoms given by particle identifier.
struct SingleCircuit {
    enum class Outcome { InvalidAtom, NotCoincidentSites, Traced };

    Outcome outcome;
    int atom1 = 0;
    int atom2 = 0;
    // A traced circuit that closed.
    std::optional<CircuitResult> circuit;
};

// The alpha-shape surface of grain B and the CIS atoms, before the non-CIS atoms are pruned from it.
struct CisSurface {
    const IldaMesh& mesh;
};

struct SurfacePruning {
    int prunedAtoms = 0;
};

// The CIS surface mesh once it is pruned, and whether its opposite half edges could all be connected.
struct CisSurfaceMesh {
    const IldaMesh& mesh;
    bool halfEdgesConnected = true;
};

// The mesh faces too large to be unwrapped across a periodic boundary, which the defect search skips.
struct PeriodicFaces {
    std::size_t skipped = 0;
};

// The Burgers circuits traced along the CIS mesh edges.
struct CircuitTally {
    std::size_t traced = 0;
    std::size_t aboveTolerance = 0;
    std::size_t failed = 0;
};

// The line defects found on the CIS mesh, grouped into disconnection modes. Without a terrace plane normal the ideal
// step heights are unknown.
struct DefectLines {
    const std::vector<IldaSegment>& segments;
    const std::vector<DisconnectionMode>& modes;
    double totalLength = 0.0;
    bool stepHeightsUnknown = false;
};

struct WrittenLineSegments {
    std::string path;
};

struct WrittenSummary {
    const std::vector<DisconnectionMode>& modes;
    double totalLength = 0.0;
    std::string path;
};

// An atom whose PTM template could not be remapped onto the reference orientation; the query then counts it as OTHER.
struct PtmTemplateRemapFailure {
    std::size_t atom = 0;
};

// A structure type whose PTM template vectors were not checked against a lattice constant; they keep its value as
// their scale.
struct UntestedPtmScale {
    StructureType type = OTHER;
};

// What made a pass stop early or warns about its input; the run goes on without what the pass would have built.
enum class IldaIssue {
    CisPositions,
    CisNeighborLists,
    InterfaceSkinPositions,
    InterfaceSkinNeighborLists,
    SurfaceInput,
    SurfacePositions,
    SurfaceConstruction,
    MixedCisMeshEdge,
};

}
