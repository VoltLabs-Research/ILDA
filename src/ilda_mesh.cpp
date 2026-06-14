#include <volt/ilda_mesh.h>

#include <volt/core/frame_adapter.h>
#include <volt/pipeline/delaunay_tessellation.h>
#include <volt/helpers/manifold_construction_helper.h>

#include <array>
#include <set>
#include <vector>

#include <spdlog/spdlog.h>

namespace Volt {

namespace {

void makeManifold(IldaMesh& mesh) {
    auto originalVertices = mesh.vertices();

    for(auto* vertex : originalVertices) {
        if(vertex->numEdges() < 3) continue;

        std::set<IldaMesh::Edge*> visitedEdges;
        IldaMesh::Edge* startEdge = vertex->edges();
        IldaMesh::Edge* currentEdge = startEdge;
        do {
            visitedEdges.insert(currentEdge);
            currentEdge = currentEdge->oppositeEdge()->nextFaceEdge();
        } while(currentEdge != startEdge);

        if(visitedEdges.size() == vertex->numEdges()) {
            continue;
        }

        while(visitedEdges.size() < vertex->numEdges()) {
            IldaMesh::Vertex* newVertex = mesh.createVertex(vertex->pos());
            newVertex->atomIndex = vertex->atomIndex;
            newVertex->isCIS = vertex->isCIS;
            newVertex->color = vertex->color;

            IldaMesh::Edge* fanStartEdge = nullptr;
            for(IldaMesh::Edge* e = vertex->edges(); e != nullptr; e = e->nextVertexEdge()) {
                if(visitedEdges.find(e) == visitedEdges.end()) {
                    fanStartEdge = e;
                    break;
                }
            }
            if(fanStartEdge == nullptr) break;

            std::vector<IldaMesh::Edge*> fanToTransfer;
            currentEdge = fanStartEdge;
            do {
                fanToTransfer.push_back(currentEdge);
                visitedEdges.insert(currentEdge);
                currentEdge = currentEdge->oppositeEdge()->nextFaceEdge();
            } while(currentEdge != fanStartEdge);

            for(IldaMesh::Edge* edgeToMove : fanToTransfer) {
                vertex->transferEdgeToVertex(edgeToMove, newVertex);
            }
        }
    }
}

}

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
) {
    outMesh.clear();

    const std::size_t n = static_cast<std::size_t>(frame.natoms);
    if(n == 0 || Rsphere <= 0.0) {
        spdlog::warn("buildCisSurfaceMesh: empty frame or non-positive Rsphere");
        return false;
    }

    std::vector<int> selected(n, 0);
    for(std::size_t i = 0; i < n; ++i) {
        const bool grainBmatch = (grain[i] == grainB) && (structureType[i] == g2s);
        if(grainBmatch || cis[i] == 1) {
            selected[i] = 1;
        }
    }

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions) {
        spdlog::error("buildCisSurfaceMesh: failed to create position property");
        return false;
    }

    // PARITY: CIS filter goes into `selected` (only_selected=True), NOT determineCellRegion (returns 1 for every alpha-solid tet).
    DelaunayTessellation tessellation;
    const double ghostLayerSize = 2.5 * Rsphere;
    tessellation.generateTessellation(
        cell,
        positions->constDataPoint3(),
        n,
        ghostLayerSize,
        /*coverDomainWithFiniteTets=*/false,
        selected.data()
    );

    // PARITY (alpha convention): alphaTest compares circumradius^2 < alpha, so alpha = Rsphere^2.
    const double alpha = Rsphere * Rsphere;
    ManifoldConstructionHelper<IldaMesh> helper{tessellation, outMesh, alpha, positions.get()};

    auto determineCellRegion = [](DelaunayTessellation::CellHandle) -> unsigned { return 1; };

    auto prepareMeshFace = [&](
        IldaMesh::Face* face,
        const std::array<int, 3>& vIdx,
        const std::array<DelaunayTessellation::VertexHandle, 3>& /*vH*/,
        DelaunayTessellation::CellHandle /*cellHandle*/
    ) {
        IldaMesh::Edge* e = face->edges();
        for(int i = 0; i < 3; ++i) {
            IldaMesh::Vertex* v = e->vertex1();
            const int atomIdx = vIdx[static_cast<std::size_t>(i)];
            v->atomIndex = atomIdx;
            if(atomIdx >= 0 && static_cast<std::size_t>(atomIdx) < n) {
                v->isCIS = cis[static_cast<std::size_t>(atomIdx)];
            }
            e = e->nextFaceEdge();
        }
    };

    if(!helper.construct(determineCellRegion, prepareMeshFace)) {
        spdlog::error("buildCisSurfaceMesh: ManifoldConstructionHelper::construct failed");
        return false;
    }

    makeManifold(outMesh);

    const std::size_t numVerts = outMesh.vertexCount();
    std::vector<Point3> vertexPos(numVerts);
    std::vector<int> vertexAtom(numVerts, -1);
    std::vector<std::uint8_t> vertexCis(numVerts, 0);
    for(std::size_t v = 0; v < numVerts; ++v) {
        IldaMesh::Vertex* vtx = outMesh.vertex(static_cast<int>(v));
        vertexPos[v] = vtx->pos();
        vertexAtom[v] = vtx->atomIndex;
        vertexCis[v] = vtx->isCIS;
    }

    std::vector<std::array<int, 3>> faces;
    faces.reserve(outMesh.faceCount());
    for(IldaMesh::Face* f : outMesh.faces()) {
        IldaMesh::Edge* e = f->edges();
        if(!e) continue;
        std::array<int, 3> tri{
            e->vertex1()->index(),
            e->vertex2()->index(),
            e->nextFaceEdge()->vertex2()->index()
        };
        faces.push_back(tri);
    }

    auto selOf = [&](int meshVertexIndex) -> int {
        const int atom = vertexAtom[static_cast<std::size_t>(meshVertexIndex)];
        return (atom >= 0) ? 1 : 0;
    };
    auto cisOf = [&](int meshVertexIndex) -> int {
        return vertexCis[static_cast<std::size_t>(meshVertexIndex)] ? 1 : 0;
    };

    int pruneCount = 0;
    while(true) {
        // PARITY: numpy logical_and 3-arg quirk drops vertex2 — only vertex0 AND vertex1 gate selection.
        std::vector<std::size_t> selectedFaces;
        selectedFaces.reserve(faces.size());
        for(std::size_t fi = 0; fi < faces.size(); ++fi) {
            if(selOf(faces[fi][0]) == 1 && selOf(faces[fi][1]) == 1) {
                selectedFaces.push_back(fi);
            }
        }

        // PARITY: scan ind=0,1,2, take the first hit (deterministic order matches Python).
        int ind = -1;
        std::size_t hitFace = 0;
        bool found = false;
        for(int candidateInd = 0; candidateInd < 3 && !found; ++candidateInd) {
            const int other1 = (candidateInd + 1) % 3;
            const int other2 = (candidateInd + 2) % 3;
            for(std::size_t k = 0; k < selectedFaces.size(); ++k) {
                const std::size_t fi = selectedFaces[k];
                if(cisOf(faces[fi][static_cast<std::size_t>(candidateInd)]) == 0 &&
                   cisOf(faces[fi][static_cast<std::size_t>(other1)]) == 1 &&
                   cisOf(faces[fi][static_cast<std::size_t>(other2)]) == 1) {
                    ind = candidateInd;
                    hitFace = fi;
                    found = true;
                    break;
                }
            }
        }

        if(!found) {
            break;
        }

        // PARITY (Python negative indexing): faces[j, ind-1] means ind=0 picks vertex2, ind=1 picks 0, ind=2 picks 1.
        const int rmVertex = faces[hitFace][static_cast<std::size_t>(ind)];
        const int pivotVertex = faces[hitFace][static_cast<std::size_t>((ind + 2) % 3)];

        std::vector<std::array<int, 3>> newFaces;
        newFaces.reserve(faces.size());
        for(const auto& tri : faces) {
            const bool hasRm = (tri[0] == rmVertex || tri[1] == rmVertex || tri[2] == rmVertex);
            const bool hasPivot = (tri[0] == pivotVertex || tri[1] == pivotVertex || tri[2] == pivotVertex);
            if(hasRm && hasPivot) {
                continue;
            }
            std::array<int, 3> t = tri;
            for(int c = 0; c < 3; ++c) {
                if(t[static_cast<std::size_t>(c)] == rmVertex) {
                    t[static_cast<std::size_t>(c)] = pivotVertex;
                }
            }
            newFaces.push_back(t);
        }
        faces.swap(newFaces);

        ++pruneCount;
        if(pruneCount > ILDA_PRUNE_MAX) {
            spdlog::warn("Number of surface mesh pruning events exceeded {} - pruning was aborted",
                         ILDA_PRUNE_MAX);
            break;
        }
    }

    spdlog::info("Pruned {} non-CIS atoms from the surface mesh.", pruneCount);

    outMesh.clear();
    outMesh.reserveVertices(numVerts);
    for(std::size_t v = 0; v < numVerts; ++v) {
        IldaMesh::Vertex* vtx = outMesh.createVertex(vertexPos[v]);
        vtx->atomIndex = vertexAtom[v];
        vtx->isCIS = vertexCis[v];
        vtx->color = vertexCis[v] ? Vector3(1.0, 1.0, 0.0) : Vector3(0.0, 0.0, 0.0);
    }
    outMesh.reserveFaces(faces.size());
    for(const auto& tri : faces) {
        outMesh.createFace({
            outMesh.vertex(tri[0]),
            outMesh.vertex(tri[1]),
            outMesh.vertex(tri[2])
        });
    }

    if(!outMesh.connectOppositeHalfedges()) {
        spdlog::warn("Unable to connect opposite half edges of the pruned surface mesh.");
    }

    return true;
}

}
