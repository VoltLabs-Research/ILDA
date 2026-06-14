#include <volt/ilda_defects.h>
#include <volt/ilda_plane.h>
#include <volt/ilda_circuit.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

#include <spdlog/spdlog.h>

namespace Volt {

const Vector3 ILDA_COLORLIST[20] = {
    Vector3(0.0,    0.0,    0.6),
    Vector3(0.0,    0.0,    0.8),
    Vector3(0.0,    0.0,    1.0),
    Vector3(0.0,    0.2,    1.0),
    Vector3(0.0,    0.4,    1.0),
    Vector3(0.0,    0.6,    1.0),
    Vector3(0.0,    0.8,    1.0),
    Vector3(0.0,    1.0,    1.0),
    Vector3(0.2,    1.0,    0.8),
    Vector3(0.4,    1.0,    0.6),
    Vector3(0.6,    1.0,    0.4),
    Vector3(0.8,    1.0,    0.2),
    Vector3(1.0,    1.0,    0.0),
    Vector3(1.0,    0.8,    0.0),
    Vector3(1.0,    0.6,    0.0),
    Vector3(1.0,    0.4,    0.0),
    Vector3(1.0,    0.2,    0.0),
    Vector3(1.0,    0.0,    0.0),
    Vector3(0.8,    0.0,    0.0),
    Vector3(0.6,    0.0,    0.0),
};

namespace {

constexpr double ILDA_EPS = 1e-6;

double vecAngle(const Vector3& v1, const Vector3& v2) {
    const double l1 = std::sqrt(v1.dot(v1));
    const double l2 = std::sqrt(v2.dot(v2));
    double c = v1.dot(v2) / (l1 * l2);
    if(c > 1.0) c = 1.0;
    if(c < -1.0) c = -1.0;
    return std::acos(c);
}

std::pair<bool, double> parallelTest(const Vector3& v1, const Vector3& v2, double angtol) {
    double ang = std::abs(vecAngle(v1, v2));
    double sign = 1.0;
    if(ang > PI / 2.0) {
        sign = -1.0;  // PARITY: fold antiparallel onto parallel (sign flip)
        ang = PI - ang;
    }
    return { ang < angtol, sign };
}

std::vector<double> arange(double start, double stop, double step) {
    std::vector<double> out;
    if(step <= 0.0) return out;
    for(int i = 0; ; ++i) {
        const double v = start + static_cast<double>(i) * step;
        if(v >= stop) break;
        out.push_back(v);
    }
    return out;
}

std::vector<double> binMids(const std::vector<double>& edges) {
    std::vector<double> mids;
    if(edges.size() < 2) return mids;
    mids.reserve(edges.size() - 1);
    for(std::size_t i = 0; i + 1 < edges.size(); ++i) {
        mids.push_back((edges[i] + edges[i + 1]) / 2.0);
    }
    return mids;
}

int findBin(double v, const std::vector<double>& edges) {
    const int L = static_cast<int>(edges.size());
    if(L < 2) return -1;
    if(v < edges[0]) return -1;
    if(v > edges[static_cast<std::size_t>(L - 1)]) return -1;
    for(int b = 0; b < L - 1; ++b) {
        if(v < edges[static_cast<std::size_t>(b + 1)]) return b;
    }
    return L - 2;  // PARITY: last bin is closed [e[L-2], e[L-1]] (numpy.histogramdd)
}

int argminAbsDiff(const std::vector<double>& mids, double val) {
    if(mids.empty()) return 0;
    if(std::isnan(val)) return 0;  // PARITY: numpy.argmin returns 0 for all-NaN
    int best = 0;
    double bd = std::abs(mids[0] - val);
    for(std::size_t i = 1; i < mids.size(); ++i) {
        const double d = std::abs(mids[i] - val);
        if(d < bd) {  // PARITY: strict < keeps ties at the first index
            bd = d;
            best = static_cast<int>(i);
        }
    }
    return best;
}

double interp01(double x, const std::vector<double>& fp) {
    const int n = static_cast<int>(fp.size());
    if(n == 0) return 0.0;
    if(n == 1) return fp[0];
    if(x <= 0.0) return fp[0];
    if(x >= 1.0) return fp[static_cast<std::size_t>(n - 1)];
    const double scaled = x * static_cast<double>(n - 1);
    int lo = static_cast<int>(std::floor(scaled));
    if(lo >= n - 1) lo = n - 2;
    const double frac = scaled - static_cast<double>(lo);
    return fp[static_cast<std::size_t>(lo)] +
           frac * (fp[static_cast<std::size_t>(lo + 1)] - fp[static_cast<std::size_t>(lo)]);
}

int getId2(double bnorm, double hA, double hB, double btol,
           const std::vector<double>& bbinmids,
           const std::vector<double>& hbinmidsA,
           const std::vector<double>& hbinmidsB,
           const std::vector<int>& bhhistIds,
           int nHA, int nHB) {
    if(bnorm > btol) {
        const int bmatch = argminAbsDiff(bbinmids, bnorm);
        const int hmatchA = argminAbsDiff(hbinmidsA, hA);
        const int hmatchB = argminAbsDiff(hbinmidsB, hB);
        const int idx = (bmatch * nHA + hmatchA) * nHB + hmatchB;
        if(idx < 0 || idx >= static_cast<int>(bhhistIds.size())) return -1;
        return bhhistIds[static_cast<std::size_t>(idx)];
    }
    return -1;
}

enum class Edge1Kind { NoAtom, UsePcount, Real };

void setMidedge(IldaMesh::Edge* e, int value) {
    e->midedgeAtom = value;
    if(e->oppositeEdge()) {
        e->oppositeEdge()->midedgeAtom = value;
    }
}

std::pair<int, int> segmentAtom(Edge1Kind kind, IldaMesh::Edge* edge1,
                                IldaMesh::Edge* edge2, int pcount) {
    int atom1;
    if(kind == Edge1Kind::NoAtom) {
        atom1 = -1;
    } else if(kind == Edge1Kind::Real) {
        if(edge1->midedgeAtom > 0) {
            atom1 = edge1->midedgeAtom;
        } else {
            atom1 = pcount;
            setMidedge(edge1, pcount);
        }
    } else {
        atom1 = pcount;
    }

    int atom2;
    if(edge2->midedgeAtom > 0) {
        atom2 = edge2->midedgeAtom;
    } else if(atom1 == pcount) {
        atom2 = pcount + 1;
        setMidedge(edge2, pcount + 1);
    } else {
        atom2 = pcount;
        setMidedge(edge2, pcount);
    }
    return { atom1, atom2 };
}

}

void findDefects(
    IldaMesh& mesh,
    const CircuitContext& ctx,
    double btol,
    double htol,
    double angtol,
    double bmagmax,
    std::vector<IldaSegment>& outSegments,
    std::vector<DisconnectionMode>& outModes,
    double& outTotalLength
) {
    outSegments.clear();
    outModes.clear();
    outTotalLength = 0.0;

    const SimulationCell& cell = ctx.cell;
    const CisResult& cis = ctx.cis;
    const std::vector<std::uint8_t>& selection = ctx.selection;
    const std::size_t particleCount = ctx.positions.size();

    for(IldaMesh::Face* f : mesh.faces()) {
        IldaMesh::Edge* e = f->edges();
        if(!e) continue;
        do {
            e->burgers = Vector3(0.0, 0.0, 0.0);
            e->bnorm = 0.0;
            e->plane = Vector3(0.0, 0.0, 0.0);
            e->step = 0.0;
            e->stepIdealA = 0.0;
            e->stepIdealB = 0.0;
            e->bid = -1;
            e->midedgeAtom = 0;
            e = e->nextFaceEdge();
        } while(e != f->edges());
    }

    const std::size_t vertexCount = mesh.vertexCount();
    std::vector<std::uint8_t> pbcFailVertices(vertexCount, 0);
    const auto& pbc = cell.pbcFlags();
    bool anyPbcFail = false;
    for(IldaMesh::Face* f : mesh.faces()) {
        IldaMesh::Edge* e1 = f->edges();
        if(!e1) continue;
        IldaMesh::Edge* e2 = e1->nextFaceEdge();
        IldaMesh::Edge* e3 = e1->prevFaceEdge();
        IldaMesh::Vertex* fv[3] = { e1->vertex1(), e2->vertex1(), e3->vertex1() };

        const Point3 rv0 = cell.absoluteToReduced(fv[0]->pos());
        const Point3 rv1 = cell.absoluteToReduced(fv[1]->pos());
        const Point3 rv2 = cell.absoluteToReduced(fv[2]->pos());
        const Vector3 fe[3] = {
            Vector3(rv1.x() - rv0.x(), rv1.y() - rv0.y(), rv1.z() - rv0.z()),
            Vector3(rv2.x() - rv1.x(), rv2.y() - rv1.y(), rv2.z() - rv1.z()),
            Vector3(rv0.x() - rv2.x(), rv0.y() - rv2.y(), rv0.z() - rv2.z())
        };
        for(int dim = 0; dim < 3; ++dim) {
            if(!pbc[static_cast<std::size_t>(dim)]) continue;
            int longCount = 0;
            for(int i = 0; i < 3; ++i) {
                if(std::abs(fe[i][static_cast<std::size_t>(dim)]) >= 0.5) ++longCount;
            }
            if(longCount != 0 && longCount != 2) {
                for(int i = 0; i < 3; ++i) {
                    const int vi = fv[i]->index();
                    if(vi >= 0 && vi < static_cast<int>(vertexCount)) {
                        pbcFailVertices[static_cast<std::size_t>(vi)] = 1;
                    }
                }
                anyPbcFail = true;
            }
        }
    }
    if(anyPbcFail) {
        spdlog::warn("Some mesh faces are too large to use periodic boundary "
                     "conditions, these have been skipped by ILDA");
    }

    auto faceFail = [&](IldaMesh::Vertex* v) -> bool {
        const int vi = v->index();
        return vi >= 0 && vi < static_cast<int>(vertexCount) &&
               pbcFailVertices[static_cast<std::size_t>(vi)] == 1;
    };
    auto selOf = [&](int atom) -> bool {
        return atom >= 0 && static_cast<std::size_t>(atom) < selection.size() &&
               selection[static_cast<std::size_t>(atom)] != 0;
    };
    auto cisOf = [&](int atom) -> bool {
        return atom >= 0 && static_cast<std::size_t>(atom) < cis.cis.size() &&
               cis.cis[static_cast<std::size_t>(atom)] == 1;
    };

    bool printedMeshWarning = false;
    for(IldaMesh::Face* f : mesh.faces()) {
        IldaMesh::Edge* e1 = f->edges();
        if(!e1) continue;
        IldaMesh::Edge* e2 = e1->nextFaceEdge();
        IldaMesh::Edge* e3 = e1->prevFaceEdge();
        IldaMesh::Edge* faceEdges[3] = { e1, e2, e3 };
        for(IldaMesh::Edge* edge : faceEdges) {
            IldaMesh::Edge* opp = edge->oppositeEdge();
            if(opp && opp->bnorm > 0.0) {
                edge->burgers = -opp->burgers;  // PARITY: opposite solved -> flip Burgers sign
                edge->bnorm = opp->bnorm;
                edge->stepIdealA = opp->stepIdealA;
                edge->stepIdealB = opp->stepIdealB;
            } else {
                IldaMesh::Vertex* v1 = edge->vertex1();
                IldaMesh::Vertex* v2 = edge->vertex2();
                if(faceFail(v1) || faceFail(v2)) continue;
                const int atom1 = v1->atomIndex;
                const int atom2 = v2->atomIndex;
                if(!selOf(atom1) || !selOf(atom2)) continue;
                if(cisOf(atom1) && cisOf(atom2)) {
                    std::optional<CircuitResult> out = burgersCircuit(atom1, atom2, ctx);
                    if(out) {
                        edge->burgers = out->latticeSum;
                        edge->stepIdealA = out->hA;
                        edge->stepIdealB = out->hB;
                        edge->bnorm = edge->burgers.length();
                    }
                } else if(cisOf(atom1) != cisOf(atom2) && !printedMeshWarning) {
                    spdlog::warn("Mesh edge connects CIS atom to non-CIS atom, "
                                 "Rsphere probably needs to be increased");
                    printedMeshWarning = true;
                }
            }
        }
    }

    for(IldaMesh::Face* f : mesh.faces()) {
        IldaMesh::Edge* e1 = f->edges();
        if(!e1) continue;
        IldaMesh::Edge* e2 = e1->nextFaceEdge();
        IldaMesh::Edge* e3 = e1->prevFaceEdge();
        IldaMesh::Edge* faceEdges[3] = { e1, e2, e3 };
        for(IldaMesh::Edge* edge : faceEdges) {
            if(edge->bnorm > btol) {
                IldaMesh::Edge* opp = edge->oppositeEdge();
                if(opp && opp->step > 0.0) {
                    edge->plane = opp->plane;
                    edge->step = opp->step;
                } else {
                    IldaMesh::Vertex* v1 = edge->vertex1();
                    IldaMesh::Vertex* v2 = edge->vertex2();
                    if(faceFail(v1) || faceFail(v2)) continue;
                    auto [normal, h] = findTerracePlane2(mesh, v1, v2, cell, btol);
                    edge->plane = normal;
                    edge->step = h;
                }
            }
        }
    }

    double bMax = 0.0;
    double hMaxA = 0.0;
    double hMaxB = 0.0;
    for(IldaMesh::Face* f : mesh.faces()) {
        IldaMesh::Edge* e = f->edges();
        if(!e) continue;
        do {
            bMax = std::max(bMax, e->bnorm);
            hMaxA = std::max(hMaxA, e->stepIdealA);
            hMaxB = std::max(hMaxB, e->stepIdealB);
            e = e->nextFaceEdge();
        } while(e != f->edges());
    }

    // PARITY: pad stop with +eps so the intended final edge lands just inside arange.
    const std::vector<double> bbins = arange(btol, bMax + btol + ILDA_EPS, btol);
    const std::vector<double> bbinmids = binMids(bbins);
    const std::vector<double> hbinsA = arange(0.0, hMaxA + htol + ILDA_EPS, htol);
    const std::vector<double> hbinmidsA = binMids(hbinsA);
    const std::vector<double> hbinsB = arange(0.0, hMaxB + htol + ILDA_EPS, htol);
    const std::vector<double> hbinmidsB = binMids(hbinsB);

    const int nB = static_cast<int>(bbinmids.size());
    const int nHA = static_cast<int>(hbinmidsA.size());
    const int nHB = static_cast<int>(hbinmidsB.size());

    std::vector<std::int64_t> hist;
    if(nB > 0 && nHA > 0 && nHB > 0) {
        hist.assign(static_cast<std::size_t>(nB) * nHA * nHB, 0);
        for(IldaMesh::Face* f : mesh.faces()) {
            IldaMesh::Edge* e = f->edges();
            if(!e) continue;
            do {
                const int bi = findBin(e->bnorm, bbins);
                const int ai = findBin(e->stepIdealA, hbinsA);
                const int ci = findBin(e->stepIdealB, hbinsB);
                if(bi >= 0 && ai >= 0 && ci >= 0) {
                    ++hist[static_cast<std::size_t>((bi * nHA + ai) * nHB + ci)];
                }
                e = e->nextFaceEdge();
            } while(e != f->edges());
        }
    }

    std::vector<int> bhhistIds(hist.size(), -1);
    std::vector<double> bmaglistById;
    int idCount = -1;
    for(int i = 0; i < nB; ++i) {
        for(int j = 0; j < nHA; ++j) {
            for(int k = 0; k < nHB; ++k) {
                const std::size_t idx = static_cast<std::size_t>((i * nHA + j) * nHB + k);
                if(hist[idx] > 0) {
                    ++idCount;
                    bhhistIds[idx] = idCount;
                    bmaglistById.push_back(bbinmids[static_cast<std::size_t>(i)]);
                }
            }
        }
    }
    const int Nb = idCount + 1;

    const bool noTerraceNormal = (ctx.n.length() < ILDA_EPS);
    if(noTerraceNormal) {
        spdlog::warn("Terrace plane normal in coherent reference frame, n, not "
                     "specified, cannot determine ideal step heights");
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for(IldaMesh::Face* f : mesh.faces()) {
            IldaMesh::Edge* e = f->edges();
            if(!e) continue;
            do {
                e->stepIdealA = nan;
                e->stepIdealB = nan;
                e = e->nextFaceEdge();
            } while(e != f->edges());
        }
    }

    std::vector<Vector3> bcolors;
    bcolors.reserve(static_cast<std::size_t>(std::max(0, Nb)));
    std::vector<double> cmapR(20), cmapG(20), cmapB(20);
    for(int c = 0; c < 20; ++c) {
        cmapR[static_cast<std::size_t>(c)] = ILDA_COLORLIST[c].x();
        cmapG[static_cast<std::size_t>(c)] = ILDA_COLORLIST[c].y();
        cmapB[static_cast<std::size_t>(c)] = ILDA_COLORLIST[c].z();
    }
    for(int m = 0; m < Nb; ++m) {
        double bmagnorm = (bmagmax != 0.0) ? bmaglistById[static_cast<std::size_t>(m)] / bmagmax : 0.0;
        if(bmagnorm > 1.0) bmagnorm = 1.0;
        bcolors.emplace_back(interp01(bmagnorm, cmapR),
                             interp01(bmagnorm, cmapG),
                             interp01(bmagnorm, cmapB));
    }
    auto colorFor = [&](int bid) -> Vector3 {
        if(bid >= 0 && bid < Nb) return bcolors[static_cast<std::size_t>(bid)];
        if(Nb > 0) return bcolors[static_cast<std::size_t>(Nb - 1)];  // PARITY: bcolors[-1]
        return Vector3(0.0, 0.0, 0.0);
    };

    struct VariantAccum {
        std::vector<Vector3> protos;
        std::vector<Vector3> sums;
        std::vector<Vector3> sum2s;
        std::vector<double> hsums;
        std::vector<double> lsums;
        std::vector<double> hidealAsums;
        std::vector<double> hidealBsums;
    };
    std::vector<VariantAccum> acc(static_cast<std::size_t>(std::max(0, Nb)));

    int pcount = static_cast<int>(particleCount);
    std::int64_t bondId = 0;
    auto addSegment = [&](int atom1, int atom2, const Point3& pos1, const Point3& pos2,
                          const Vector3& b, double bnorm, double h,
                          double hidealA, double hidealB, const Vector3& plane,
                          const Vector3& bcolor, int bid) {
        if(atom1 >= pcount) pcount = std::max(pcount, atom1 + 1);
        if(atom2 >= pcount) pcount = std::max(pcount, atom2 + 1);
        IldaSegment seg;
        seg.atom1 = atom1;
        seg.atom2 = atom2;
        seg.pos1 = pos1;
        seg.pos2 = pos2;
        seg.burgers = b;
        seg.bnorm = bnorm;
        seg.step = h;
        seg.stepIdealA = std::abs(hidealA);
        seg.stepIdealB = std::abs(hidealB);
        seg.plane = plane;
        seg.color = bcolor;
        seg.bid = bid + 1;  // PARITY: 1-based ids
        (void)bondId;
        outSegments.push_back(seg);
    };

    for(IldaMesh::Face* f : mesh.faces()) {
        IldaMesh::Edge* e1 = f->edges();
        if(!e1) continue;
        IldaMesh::Edge* e2 = e1->nextFaceEdge();
        IldaMesh::Edge* e3 = e1->prevFaceEdge();

        IldaMesh::Vertex* v1 = e1->vertex1();
        IldaMesh::Vertex* v2 = e2->vertex1();
        IldaMesh::Vertex* v3 = e3->vertex1();

        if(faceFail(v1) || faceFail(v2) || faceFail(v3)) continue;

        const int atom1 = v1->atomIndex;
        const int atom2 = v2->atomIndex;
        const int atom3 = v3->atomIndex;

        if(!cisOf(atom1) || !cisOf(atom2) || !cisOf(atom3)) continue;
        if(!selOf(atom1) || !selOf(atom2) || !selOf(atom3)) continue;

        const Vector3 b1 = e1->burgers;
        const Vector3 b2 = e2->burgers;
        const Vector3 b3 = e3->burgers;
        const double h1 = e1->step;
        const double h2 = e2->step;
        const double h3 = e3->step;
        const double hidealA1 = e1->stepIdealA;
        const double hidealA2 = e2->stepIdealA;
        const double hidealA3 = e3->stepIdealA;
        const double hidealB1 = e1->stepIdealB;
        const double hidealB2 = e2->stepIdealB;
        const double hidealB3 = e3->stepIdealB;
        const Vector3 plane1 = e1->plane;
        const Vector3 plane2 = e2->plane;
        const Vector3 plane3 = e3->plane;

        const double bSumNorm = (b1 + b2 + b3).length();
        if(bSumNorm > btol) {
            v1->color = Vector3(1.0, 0.0, 0.0);
            v2->color = Vector3(1.0, 0.0, 0.0);
            v3->color = Vector3(1.0, 0.0, 0.0);
        }

        const double bnorm1 = e1->bnorm;
        const double bnorm2 = e2->bnorm;
        const double bnorm3 = e3->bnorm;

        // PARITY: NaN ideal heights never satisfy `> htol` (matches numpy).
        const bool edge1test = (bnorm1 > btol) || (hidealA1 > htol);
        const bool edge2test = (bnorm2 > btol) || (hidealA2 > htol);
        const bool edge3test = (bnorm3 > btol) || (hidealA3 > htol);

        if(edge1test || edge2test || edge3test) {
            const Point3 posv1 = v1->pos();
            const Point3 posv2 = posv1 + cell.wrapVector(v2->pos() - posv1);
            const Point3 posv3 = posv1 + cell.wrapVector(v3->pos() - posv1);
            const Point3 pos1 = Point3((posv1.x() + posv2.x()) / 2.0,
                                       (posv1.y() + posv2.y()) / 2.0,
                                       (posv1.z() + posv2.z()) / 2.0);
            const Point3 pos2 = Point3((posv2.x() + posv3.x()) / 2.0,
                                       (posv2.y() + posv3.y()) / 2.0,
                                       (posv2.z() + posv3.z()) / 2.0);
            const Point3 pos3 = Point3((posv1.x() + posv3.x()) / 2.0,
                                       (posv1.y() + posv3.y()) / 2.0,
                                       (posv1.z() + posv3.z()) / 2.0);

            const int bid1 = getId2(bnorm1, hidealA1, hidealB1, btol, bbinmids, hbinmidsA,
                                    hbinmidsB, bhhistIds, nHA, nHB);
            const int bid2 = getId2(bnorm2, hidealA2, hidealB2, btol, bbinmids, hbinmidsA,
                                    hbinmidsB, bhhistIds, nHA, nHB);
            const int bid3 = getId2(bnorm3, hidealA3, hidealB3, btol, bbinmids, hbinmidsA,
                                    hbinmidsB, bhhistIds, nHA, nHB);

            double L1 = 0.0, L2 = 0.0, L3 = 0.0;

            if(edge1test && edge2test && edge3test) {
                const Point3 mid = Point3((pos1.x() + pos2.x() + pos3.x()) / 3.0,
                                          (pos1.y() + pos2.y() + pos3.y()) / 3.0,
                                          (pos1.z() + pos2.z() + pos3.z()) / 3.0);
                L1 = (pos1 - mid).length();
                auto [a1, a2] = segmentAtom(Edge1Kind::UsePcount, nullptr, e1, pcount);
                const int center = a1;
                addSegment(a1, a2, mid, pos1, -b1, bnorm1, h1, hidealA1, hidealB1, plane1,
                           colorFor(bid1), bid1);

                L2 = (pos2 - mid).length();
                auto [ignore2, b2node] = segmentAtom(Edge1Kind::NoAtom, nullptr, e2, pcount);
                (void)ignore2;
                addSegment(center, b2node, mid, pos2, -b2, bnorm2, h2, hidealA2, hidealB2,
                           plane2, colorFor(bid2), bid2);

                L3 = (pos3 - mid).length();
                auto [ignore3, b3node] = segmentAtom(Edge1Kind::NoAtom, nullptr, e3, pcount);
                (void)ignore3;
                addSegment(center, b3node, mid, pos3, -b3, bnorm3, h3, hidealA3, hidealB3,
                           plane3, colorFor(bid3), bid3);
            } else if(edge1test && edge2test) {
                L1 = (pos1 - pos2).length();
                auto [a1, a2] = segmentAtom(Edge1Kind::Real, e1, e2, pcount);
                addSegment(a1, a2, pos1, pos2, b1, bnorm1, h1, hidealA1, hidealB1, plane1,
                           colorFor(bid1), bid1);
            } else if(edge2test && edge3test) {
                L2 = (pos2 - pos3).length();
                auto [a1, a2] = segmentAtom(Edge1Kind::Real, e2, e3, pcount);
                addSegment(a1, a2, pos2, pos3, b2, bnorm2, h2, hidealA2, hidealB2, plane2,
                           colorFor(bid2), bid2);
            } else if(edge1test && edge3test) {
                L3 = (pos1 - pos3).length();
                auto [a1, a2] = segmentAtom(Edge1Kind::Real, e3, e1, pcount);
                addSegment(a1, a2, pos3, pos1, b3, bnorm3, h3, hidealA3, hidealB3, plane3,
                           colorFor(bid3), bid3);
            }

            // PARITY: stats ordered/signed to the first prototype per bid.
            const Vector3 bs[3] = { b1, b2, b3 };
            const int bids[3] = { bid1, bid2, bid3 };
            const double Ls[3] = { L1, L2, L3 };
            const double hs[3] = { h1, h2, h3 };
            const double hAs[3] = { hidealA1, hidealA2, hidealA3 };
            const double hBs[3] = { hidealB1, hidealB2, hidealB3 };
            for(int t = 0; t < 3; ++t) {
                const int bidi = bids[t];
                const double Li = Ls[t];
                if(bidi >= 0 && bidi < Nb && Li > 0.0) {
                    VariantAccum& va = acc[static_cast<std::size_t>(bidi)];
                    const Vector3 bi = bs[t];
                    bool addNew = false;
                    double sign = 1.0;
                    std::size_t j = 0;
                    if(va.lsums.empty()) {
                        addNew = true;
                        sign = 1.0;
                        j = 0;
                    } else {
                        bool isParallel = false;
                        for(j = 0; j < va.protos.size(); ++j) {
                            auto pr = parallelTest(bi, va.protos[j], angtol);
                            isParallel = pr.first;
                            sign = pr.second;
                            if(isParallel) break;
                        }
                        if(!isParallel) {
                            addNew = true;
                            sign = 1.0;
                            j = va.protos.size();
                        }
                    }
                    if(addNew) {
                        va.protos.push_back(bi);
                        va.sums.emplace_back(0.0, 0.0, 0.0);
                        va.sum2s.emplace_back(0.0, 0.0, 0.0);
                        va.hsums.push_back(0.0);
                        va.lsums.push_back(0.0);
                        va.hidealAsums.push_back(0.0);
                        va.hidealBsums.push_back(0.0);
                    }
                    const Vector3 bstat = sign * bi;
                    va.sums[j] += Li * bstat;
                    va.sum2s[j] += Vector3(Li * bstat.x() * bstat.x(),
                                           Li * bstat.y() * bstat.y(),
                                           Li * bstat.z() * bstat.z());
                    if(!std::isnan(hs[t])) {
                        va.hsums[j] += Li * hs[t];
                    }
                    va.lsums[j] += Li;
                    va.hidealAsums[j] += Li * hAs[t];
                    va.hidealBsums[j] += Li * hBs[t];
                }
            }
        }
    }

    for(int i = 0; i < Nb; ++i) {
        VariantAccum& va = acc[static_cast<std::size_t>(i)];
        std::size_t j = 0;
        while(j < va.protos.size()) {
            const Vector3 bavgj = va.sums[j] / va.lsums[j];
            std::size_t k = j + 1;
            while(k < va.protos.size()) {
                const Vector3 bavgk = va.sums[k] / va.lsums[k];
                auto pr = parallelTest(bavgj, bavgk, angtol);
                if(pr.first) {
                    const double sign = pr.second;
                    va.sums[j] += sign * va.sums[k];
                    va.sum2s[j] += va.sum2s[k];
                    va.hsums[j] += va.hsums[k];
                    va.lsums[j] += va.lsums[k];
                    va.hidealAsums[j] += va.hidealAsums[k];
                    va.hidealBsums[j] += va.hidealBsums[k];
                    va.sums.erase(va.sums.begin() + static_cast<std::ptrdiff_t>(k));
                    va.sum2s.erase(va.sum2s.begin() + static_cast<std::ptrdiff_t>(k));
                    va.hsums.erase(va.hsums.begin() + static_cast<std::ptrdiff_t>(k));
                    va.lsums.erase(va.lsums.begin() + static_cast<std::ptrdiff_t>(k));
                    va.protos.erase(va.protos.begin() + static_cast<std::ptrdiff_t>(k));
                    va.hidealAsums.erase(va.hidealAsums.begin() + static_cast<std::ptrdiff_t>(k));
                    va.hidealBsums.erase(va.hidealBsums.begin() + static_cast<std::ptrdiff_t>(k));
                } else {
                    ++k;
                }
            }
            ++j;
        }
    }

    double Ltotal = 0.0;
    outModes.reserve(static_cast<std::size_t>(std::max(0, Nb)));
    for(int i = 0; i < Nb; ++i) {
        VariantAccum& va = acc[static_cast<std::size_t>(i)];
        DisconnectionMode mode;
        mode.id = i + 1;
        for(std::size_t j = 0; j < va.protos.size(); ++j) {
            const double Lj = va.lsums[j];
            const Vector3 bavg = va.sums[j] / Lj;
            IldaVariant variant;
            variant.bAvg = bavg;
            variant.bAvgNorm = bavg.length();
            variant.bStdX = std::sqrt(std::abs(Lj * va.sum2s[j].x() - va.sums[j].x() * va.sums[j].x())) / Lj;
            variant.bStdY = std::sqrt(std::abs(Lj * va.sum2s[j].y() - va.sums[j].y() * va.sums[j].y())) / Lj;
            variant.bStdZ = std::sqrt(std::abs(Lj * va.sum2s[j].z() - va.sums[j].z() * va.sums[j].z())) / Lj;
            variant.length = Lj;
            variant.hIdealA = va.hidealAsums[j] / Lj;
            variant.hIdealB = va.hidealBsums[j] / Lj;
            variant.hAvg = va.hsums[j] / Lj;
            mode.variants.push_back(variant);
            Ltotal += Lj;
        }
        outModes.push_back(std::move(mode));
    }
    outTotalLength = Ltotal;

    spdlog::info("Found {} unique disconnection modes (btol={}, htol={}); "
                 "total line length = {:.3f}", Nb, btol, htol, Ltotal);
}

}
