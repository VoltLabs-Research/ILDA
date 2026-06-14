#include <volt/ilda_circuit.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

#include <spdlog/spdlog.h>

namespace Volt {

namespace {

struct Neigh {
    double dist = 0.0;
    int index = -1;
    Vector3 delta{0.0, 0.0, 0.0};
    Vector3 a1{0.0, 0.0, 0.0};
    Matrix3 Ft2new{Matrix3::Identity()};
    Quaternion refOri2new{0.0, 0.0, 0.0, 0.0};
};

struct Frame {
    int atomIndex = -1;
    Vector3 spatialSum{0.0, 0.0, 0.0};
    Vector3 latticeSum{0.0, 0.0, 0.0};
    Matrix3 Ft2{Matrix3::Identity()};
    Quaternion refOri2{0.0, 0.0, 0.0, 0.0};
    int ntry = 0;                        // PARITY: depth value RECEIVED by the call (incremented on recurse)
    bool entered = false;

    Matrix3 Ft{Matrix3::Identity()};
    int grainId = 0;
    StructureType gs = OTHER;
    Quaternion refOri{0.0, 0.0, 0.0, 1.0};
    double ptmScale = 1.0;

    std::vector<Neigh> neighs;
    std::size_t nextNeighbor = 0;
};

enum class Ret { DEADEND, ABORT, SUCCESS };

struct PathResult {
    Vector3 spatialSum{0.0, 0.0, 0.0};
    Vector3 latticeSum{0.0, 0.0, 0.0};
};

}  // namespace

// PARITY: getHcpFccTrans unverified against the OVITO oracle (only fires at coherent FCC-HCP interfaces, never exercised by the parity test)
std::optional<std::pair<Matrix3, Quaternion>> getHcpFccTrans(
    const IldaPtmNeighborQuery::Neighbor& fccNeigh,
    int fccIndex,
    IldaPtmNeighborQuery& fccQuery,
    int hcpIndex,
    const CircuitContext& ctx
) {
    (void)ctx;
    const double tol = 1e-5;

    // PARITY: snapshot FCC neighbors before the HCP query clobbers the shared query (Python uses two distinct Query objects)
    std::vector<IldaPtmNeighborQuery::Neighbor> fcc;
    fcc.reserve(static_cast<std::size_t>(fccQuery.count()));
    for(int i = 0; i < fccQuery.count(); ++i) {
        fcc.push_back(fccQuery[i]);
    }

    fccQuery.findNeighbors(static_cast<std::size_t>(hcpIndex), nullptr);
    const Quaternion hcpRefOrientation = fccQuery.orientation();
    std::vector<IldaPtmNeighborQuery::Neighbor> hcp;
    hcp.reserve(static_cast<std::size_t>(fccQuery.count()));
    for(int j = 0; j < fccQuery.count(); ++j) {
        hcp.push_back(fccQuery[j]);
    }

    const Vector3 fcc1 = fccNeigh.idealVector.normalized();
    Vector3 hcp1{0.0, 0.0, 0.0};
    Vector3 hcp1_0{0.0, 0.0, 0.0};
    bool foundHcp1 = false;
    for(const auto& hn : hcp) {
        if(hn.index == fccIndex) {
            hcp1_0 = hn.idealVector;
            hcp1 = -(hn.idealVector.normalized());
            foundHcp1 = true;
            break;
        }
    }
    if(!foundHcp1) {
        return std::nullopt;
    }

    int vecCount = 1;
    Vector3 fcc2{0.0, 0.0, 0.0}, fcc3{0.0, 0.0, 0.0};
    Vector3 hcp2{0.0, 0.0, 0.0}, hcp3{0.0, 0.0, 0.0};
    Vector3 fcccross{0.0, 0.0, 0.0}, hcpcross{0.0, 0.0, 0.0};
    for(const auto& fn : fcc) {
        for(const auto& hn : hcp) {
            if(fn.index != hn.index) {
                continue;
            }
            if(vecCount == 1) {
                const Vector3 f2 = fn.idealVector.normalized();
                const Vector3 h2 = (hn.idealVector - hcp1_0).normalized();
                if(std::abs(fcc1.dot(f2)) > tol && std::abs(hcp1.dot(h2)) > tol) {
                    fcc2 = f2;
                    hcp2 = h2;
                    fcccross = fcc1.cross(fcc2).normalized();
                    hcpcross = hcp1.cross(hcp2).normalized();
                    vecCount = 2;
                }
            } else {
                const Vector3 f3 = fn.idealVector.normalized();
                const Vector3 h3 = (hn.idealVector - hcp1_0).normalized();
                if(std::abs(f3.dot(fcccross)) > tol && std::abs(h3.dot(hcpcross)) > tol) {
                    fcc3 = f3;
                    hcp3 = h3;
                    vecCount = 3;
                    break;
                }
            }
        }
        if(vecCount == 3) {
            break;
        }
    }
    if(vecCount < 3) {
        return std::nullopt;
    }

    const Matrix3 H(
        hcp1.x(), hcp1.y(), hcp1.z(),
        hcp2.x(), hcp2.y(), hcp2.z(),
        hcp3.x(), hcp3.y(), hcp3.z()
    );
    Matrix3 Hinv{Matrix3::Zero()};
    if(!H.inverse(Hinv)) {
        return std::nullopt;
    }
    Matrix3 Ft{Matrix3::Zero()};
    for(int r = 0; r < 3; ++r) {
        const Vector3 rhs(fcc1[static_cast<std::size_t>(r)],
                          fcc2[static_cast<std::size_t>(r)],
                          fcc3[static_cast<std::size_t>(r)]);
        const Vector3 g = Hinv * rhs;
        Ft(static_cast<Matrix3::size_type>(r), 0) = g.x();
        Ft(static_cast<Matrix3::size_type>(r), 1) = g.y();
        Ft(static_cast<Matrix3::size_type>(r), 2) = g.z();
    }

    return std::make_pair(Ft, hcpRefOrientation);
}

namespace {

bool computeGrainProps(Frame& f, const CircuitContext& ctx) {
    f.grainId = ctx.grain[static_cast<std::size_t>(f.atomIndex)];
    if(f.grainId == ctx.gA.grainId) {
        f.refOri = ctx.gA.orientation;
        f.gs = ctx.gA.structureType;
        f.Ft = ctx.gA.F;
        f.ptmScale = scalePTMVec(ctx.gA.structureType, ctx.gA.a);
        return true;
    }
    if(f.grainId == ctx.gB.grainId) {
        f.refOri = ctx.gB.orientation;
        f.gs = ctx.gB.structureType;
        f.Ft = ctx.gB.F;
        f.ptmScale = scalePTMVec(ctx.gB.structureType, ctx.gB.a);
        return true;
    }
    return false;
}

void buildNeighs(Frame& f, int endIndex, const Vector3& target,
                 const std::vector<std::uint8_t>& visited,
                 const CircuitContext& ctx) {
    IldaPtmNeighborQuery& query = ctx.query;
    const auto& grain = ctx.grain;
    const auto& structureType = ctx.structureType;
    const auto& selection = ctx.selection;

    // PARITY: on a stacking fault use the inherited ref_orientation2, else the grain reference
    const StructureType at = static_cast<StructureType>(structureType[static_cast<std::size_t>(f.atomIndex)]);
    const Quaternion* qref = &f.refOri;
    if((f.gs == FCC && at == HCP) || (f.gs == HCP && at == FCC)) {
        qref = &f.refOri2;
    }
    query.findNeighbors(static_cast<std::size_t>(f.atomIndex), qref);

    // PARITY: snapshot neighbors so a getHcpFccTrans call (re-targets the shared query) cannot corrupt this iteration
    std::vector<IldaPtmNeighborQuery::Neighbor> snap;
    snap.reserve(static_cast<std::size_t>(query.count()));
    for(int j = 0; j < query.count(); ++j) {
        snap.push_back(query[j]);
    }

    std::vector<Neigh> neighs;
    neighs.reserve(snap.size());

    for(const auto& neigh : snap) {
        // PARITY: if neighbor is the end atom, take it (dist 0) with inherited Ft2/ref_orientation2 and stop scanning
        if(neigh.index == endIndex) {
            Neigh nn;
            nn.dist = 0.0;
            nn.index = neigh.index;
            nn.delta = neigh.delta;
            nn.a1 = neigh.idealVector * f.ptmScale;
            nn.Ft2new = f.Ft2;
            nn.refOri2new = f.refOri2;
            neighs.push_back(nn);
            break;
        }

        if(neigh.index < 0 || static_cast<std::size_t>(neigh.index) >= grain.size()) {
            continue;
        }
        const std::size_t ni = static_cast<std::size_t>(neigh.index);
        if(selection[ni] == 0 || f.grainId != grain[ni]) {
            continue;
        }
        if(visited[ni]) {
            continue;
        }

        const StructureType nt = static_cast<StructureType>(structureType[ni]);
        Matrix3 Ft2new{Matrix3::Identity()};
        Quaternion refOri2new{0.0, 0.0, 0.0, 0.0};

        if(f.gs == FCC && nt == HCP && at == FCC) {
            query.findNeighbors(static_cast<std::size_t>(f.atomIndex), &f.refOri);
            auto tr = getHcpFccTrans(neigh, f.atomIndex, query, neigh.index, ctx);
            if(!tr) {
                continue;
            }
            Ft2new = tr->first;
            refOri2new = tr->second;
        } else if(f.gs == HCP && nt == FCC && at == HCP) {
            query.findNeighbors(static_cast<std::size_t>(f.atomIndex), &f.refOri);
            auto tr = getHcpFccTrans(neigh, f.atomIndex, query, neigh.index, ctx);
            if(!tr) {
                continue;
            }
            Ft2new = tr->first;
            refOri2new = tr->second;
        } else if(nt != f.gs && nt != at) {
            continue;
        } else if(f.gs == FCC && nt == FCC && at == HCP) {
            Ft2new = Matrix3::Identity();
            refOri2new = Quaternion(0.0, 0.0, 0.0, 0.0);
        } else if(f.gs == HCP && nt == HCP && at == FCC) {
            Ft2new = Matrix3::Identity();
            refOri2new = Quaternion(0.0, 0.0, 0.0, 0.0);
        } else {
            Ft2new = f.Ft2;
            refOri2new = f.refOri2;
        }

        Neigh nn;
        nn.dist = (f.spatialSum + neigh.delta - target).length();
        nn.index = neigh.index;
        nn.delta = neigh.delta;
        nn.a1 = neigh.idealVector * f.ptmScale;
        nn.Ft2new = Ft2new;
        nn.refOri2new = refOri2new;
        neighs.push_back(nn);
    }

    // PARITY: tie-break by atom index on equal distance (Python tuple sort fallback; deterministic)
    std::sort(neighs.begin(), neighs.end(), [](const Neigh& a, const Neigh& b) {
        if(a.dist != b.dist) {
            return a.dist < b.dist;
        }
        return a.index < b.index;
    });

    f.neighs = std::move(neighs);
    f.nextNeighbor = 0;
}

// PARITY: prev_atoms is a monotonic per-search visited array (Python list shared by reference, only appended, never popped); -1 and None both map to nullopt
std::optional<PathResult> recursiveCircuitSearch(
    int startNbr, int startIndex, int endIndex, const Vector3& target,
    const CircuitContext& ctx
) {
    const std::size_t n = ctx.positions.size();
    std::vector<std::uint8_t> visited(n, 0);

    std::vector<Frame> stack;
    {
        Frame root;
        root.atomIndex = startNbr;
        root.spatialSum = Vector3(0.0, 0.0, 0.0);
        root.latticeSum = Vector3(0.0, 0.0, 0.0);
        root.Ft2 = Matrix3::Identity();
        root.refOri2 = Quaternion(0.0, 0.0, 0.0, 0.0);
        root.ntry = 0;
        stack.push_back(std::move(root));
    }

    Ret childRet = Ret::DEADEND;
    PathResult success;

    while(!stack.empty()) {
        Frame& f = stack.back();

        if(!f.entered) {
            f.entered = true;

            if(f.atomIndex == endIndex) {
                const bool closed = (f.spatialSum - target).length() <= ILDA_CIRCUIT_TOL;
                if(closed) {
                    success.spatialSum = f.spatialSum;
                    success.latticeSum = f.latticeSum;
                    childRet = Ret::SUCCESS;
                } else {
                    childRet = Ret::DEADEND;
                }
                stack.pop_back();
                continue;
            }

            if(f.atomIndex >= 0 && static_cast<std::size_t>(f.atomIndex) < n) {
                visited[static_cast<std::size_t>(f.atomIndex)] = 1;
            }

            // PARITY: call increments ntry, so the checked value is ntry+1
            if(f.ntry + 1 > ILDA_NTRY_MAX) {
                childRet = Ret::ABORT;
                stack.pop_back();
                continue;
            }

            if(!computeGrainProps(f, ctx)) {
                childRet = Ret::DEADEND;
                stack.pop_back();
                continue;
            }

            buildNeighs(f, endIndex, target, visited, ctx);
        } else {
            if(childRet == Ret::DEADEND) {
                ++f.nextNeighbor;
            } else if(childRet == Ret::SUCCESS) {
                stack.pop_back();
                continue;
            } else {
                stack.pop_back();
                continue;
            }
        }

        if(f.nextNeighbor < f.neighs.size()) {
            const Neigh nb = f.neighs[f.nextNeighbor];

            // PARITY: aa = Ft * Ft2 * a1, Ft is THIS atom's grain F
            const Vector3 aa = f.Ft * (f.Ft2 * nb.a1);
            if(aa.length() < 0.001) {
                childRet = Ret::DEADEND;
                stack.pop_back();
                continue;
            }

            Frame child;
            child.atomIndex = nb.index;
            child.spatialSum = f.spatialSum + nb.delta;
            child.latticeSum = f.latticeSum + aa;
            child.Ft2 = nb.Ft2new;
            child.refOri2 = nb.refOri2new;
            child.ntry = f.ntry + 1;
            stack.push_back(std::move(child));
            continue;
        }

        const bool isStart = (f.atomIndex == startIndex);
        stack.pop_back();
        childRet = isStart ? Ret::ABORT : Ret::DEADEND;
    }

    if(childRet == Ret::SUCCESS) {
        return success;
    }
    return std::nullopt;
}

}  // namespace

std::optional<CircuitResult> burgersCircuit(
    int startIndex,
    int endIndex,
    const CircuitContext& ctx
) {
    const CisResult& cis = ctx.cis;

    if(startIndex < 0 || endIndex < 0 ||
       static_cast<std::size_t>(startIndex) >= cis.cis.size() ||
       static_cast<std::size_t>(endIndex) >= cis.cis.size()) {
        spdlog::warn("burgersCircuit: invalid indices {} {}", startIndex, endIndex);
        return std::nullopt;
    }
    if(!cis.cis[static_cast<std::size_t>(startIndex)] || !cis.cis[static_cast<std::size_t>(endIndex)]) {
        spdlog::warn("Wrong start or end indices: not a co-incident site {} {}", startIndex, endIndex);
        return std::nullopt;
    }
    if(startIndex == endIndex) {
        spdlog::warn("Wrong start or end indices, both are same: {}", startIndex);
        return std::nullopt;
    }

    const double ptmScales[2] = {
        scalePTMVec(ctx.gA.structureType, ctx.gA.a),
        scalePTMVec(ctx.gB.structureType, ctx.gB.a)
    };
    const Matrix3 Fs[2] = {ctx.gA.F, ctx.gB.F};

    const Point3 startPoint = ctx.positions[static_cast<std::size_t>(startIndex)];
    const Vector3 startToEnd = ctx.cell.wrapVector(
        ctx.positions[static_cast<std::size_t>(endIndex)] - ctx.positions[static_cast<std::size_t>(startIndex)]);
    const Point3 endPoint = startPoint + startToEnd;

    Vector3 halfLatticeSum{0.0, 0.0, 0.0};
    Vector3 latticeSum{0.0, 0.0, 0.0};
    Vector3 spatialSum{0.0, 0.0, 0.0};
    double hA = 0.0;
    double hB = 0.0;

    for(int i = 0; i < 2; ++i) {
        const std::size_t si = static_cast<std::size_t>(startIndex);
        const std::size_t ei = static_cast<std::size_t>(endIndex);
        const int startNbr = cis.cisNbrs[si][static_cast<std::size_t>(i)];
        const int endNbr = cis.cisNbrs[ei][static_cast<std::size_t>(i)];
        const Vector3 startDelta = cis.cisSpatialVecs[si][static_cast<std::size_t>(i)];
        const Vector3 endDelta = cis.cisSpatialVecs[ei][static_cast<std::size_t>(i)];

        const Vector3 target = (endPoint + endDelta) - (startPoint + startDelta);

        const std::optional<PathResult> circuit =
            recursiveCircuitSearch(startNbr, startNbr, endNbr, target, ctx);
        if(!circuit) {
            return std::nullopt;
        }

        const Vector3 spatialSumi = circuit->spatialSum;
        const Vector3 latticeSumi = circuit->latticeSum;

        // PARITY: signi = 1 - 2*i flips the sign for the return leg (i=1)
        const double signi = 1.0 - 2.0 * static_cast<double>(i);
        latticeSum += signi * latticeSumi;
        spatialSum += signi * spatialSumi;
        halfLatticeSum += signi * latticeSumi;

        const Vector3 startIdeal =
            Fs[i] * (cis.cisIdealVecs[si][static_cast<std::size_t>(i)] * ptmScales[i]);
        const Vector3 endIdeal =
            Fs[i] * (cis.cisIdealVecs[ei][static_cast<std::size_t>(i)] * ptmScales[i]);
        latticeSum += signi * (startIdeal - endIdeal);
        spatialSum += signi * (startDelta - endDelta);
        halfLatticeSum += signi * (startIdeal - endIdeal);

        if(i == 0) {
            hA = std::abs(halfLatticeSum.dot(ctx.n));
        } else {
            hB = std::abs(halfLatticeSum.dot(ctx.n));
        }
        halfLatticeSum = Vector3(0.0, 0.0, 0.0);
    }

    if(spatialSum.length() < ILDA_CIRCUIT_TOL) {
        CircuitResult res;
        res.latticeSum = Vector3(
            (std::abs(latticeSum.x()) < ILDA_CIRCUIT_TOL) ? 0.0 : latticeSum.x(),
            (std::abs(latticeSum.y()) < ILDA_CIRCUIT_TOL) ? 0.0 : latticeSum.y(),
            (std::abs(latticeSum.z()) < ILDA_CIRCUIT_TOL) ? 0.0 : latticeSum.z()
        );
        res.hA = (std::abs(hA) < ILDA_CIRCUIT_TOL) ? 0.0 : hA;
        res.hB = (std::abs(hB) < ILDA_CIRCUIT_TOL) ? 0.0 : hB;
        return res;
    }

    spdlog::warn("Warning: Not a closed circuit! {} {} spatialSum=({}, {}, {})",
                 startIndex, endIndex, spatialSum.x(), spatialSum.y(), spatialSum.z());
    return std::nullopt;
}

}  // namespace Volt
