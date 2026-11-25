#include <volt/ilda_mesh.h>

#include <volt/core/frame_adapter.h>
#include <volt/core/progress.h>
#include <volt/ilda_log.h>
#include <volt/pipeline/delaunay_tessellation.h>
#include <volt/helpers/manifold_construction_helper.h>

#include <array>
#include <vector>

namespace Volt {

bool buildCisSurfaceMesh(
    const SimulationCell& cell,
    const LammpsParser::Frame& frame,
    const std::vector<std::uint8_t>& cis,
    const std::vector<int>& grain,
    const std::vector<int>& structureType,
    int grainB,
    int g2s,
    double Rsphere,
    IldaMesh& outMesh,
    IldaLog& log
) {
    outMesh.clear();

    const std::size_t n = static_cast<std::size_t>(frame.natoms);
    if(n == 0 || Rsphere <= 0.0) {
        log.dump(IldaIssue::SurfaceInput);
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
        log.dump(IldaIssue::SurfacePositions);
        return false;
    }

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

    const double alpha = Rsphere * Rsphere;
    ManifoldConstructionHelper<IldaMesh> helper{tessellation, outMesh, alpha, positions.get()};

    auto determineCellRegion = [](DelaunayTessellation::CellHandle) -> unsigned { return 1; };

    auto prepareMeshFace = [&](
        IldaMesh::Face* face,
        const std::array<int, 3>& vIdx,
        const std::array<DelaunayTessellation::VertexHandle, 3>& /*vH*/,
        DelaunayTessellation::CellHandle
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
        log.dump(IldaIssue::SurfaceConstruction);
        return false;
    }

    outMesh.makeManifold();
    log.dump(CisSurface{outMesh});

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
        std::vector<std::size_t> selectedFaces;
        selectedFaces.reserve(faces.size());
        for(std::size_t fi = 0; fi < faces.size(); ++fi) {
            if(selOf(faces[fi][0]) == 1 && selOf(faces[fi][1]) == 1) {
                selectedFaces.push_back(fi);
            }
        }

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
            break;
        }
    }
    log.dump(SurfacePruning{pruneCount});

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

    const bool halfEdgesConnected = outMesh.connectOppositeHalfedges();
    log.dump(CisSurfaceMesh{outMesh, halfEdgesConnected});
    return true;
}

}
