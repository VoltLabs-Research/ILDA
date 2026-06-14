#include <volt/ilda_plane.h>

#include <cmath>
#include <limits>

namespace Volt {

std::pair<Vector3, double> findTerracePlane2(
    IldaMesh& /*mesh*/,
    IldaMesh::Vertex* startVertex,
    IldaMesh::Vertex* endVertex,
    const SimulationCell& cell,
    double btol
) {
    Vector3 areaVec(0.0, 0.0, 0.0);

    // Parity: start vertex only (Python `for i in range(1)`) — do not also walk end_vertex.
    for(IldaMesh::Edge* edge = startVertex->edges(); edge; edge = edge->nextVertexEdge()) {
        IldaMesh::Edge* edge2 = edge->nextFaceEdge();
        IldaMesh::Edge* edge3 = edge->prevFaceEdge();

        IldaMesh::Vertex* vertex1 = edge->vertex1();
        IldaMesh::Vertex* vertex2 = edge2->vertex1();
        IldaMesh::Vertex* vertex3 = edge3->vertex1();

        const double bnorm1 = edge->bnorm;
        const double bnorm2 = edge2->bnorm;
        const double bnorm3 = edge3->bnorm;

        if(bnorm1 < btol && bnorm2 < btol && bnorm3 < btol) {
            const Point3& posv1 = vertex1->pos();
            const Vector3 v12 = cell.wrapVector(vertex2->pos() - posv1);
            const Vector3 v13 = cell.wrapVector(vertex3->pos() - posv1);

            areaVec += v12.cross(v13);
        }
    }

    Vector3 n(0.0, 0.0, 0.0);
    if(areaVec.length() != 0.0) {
        n = areaVec / areaVec.length();
    }

    double h;
    if(n.length() == 0.0) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        h = nan;
        n = Vector3(nan, nan, nan);  // PARITY: zero area -> NaN normal + NaN height
    } else {
        const Vector3 disp = cell.wrapVector(endVertex->pos() - startVertex->pos());
        h = std::abs(n.dot(disp));
    }

    return { n, h };
}

}
