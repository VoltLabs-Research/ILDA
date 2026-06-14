#include <volt/ilda_cis.h>

#include <volt/core/frame_adapter.h>
#include <volt/analysis/nearest_neighbor_finder.h>

#include <array>
#include <cmath>
#include <vector>

#include <spdlog/spdlog.h>

namespace Volt {

namespace {

// PARITY: Volt {x,y,z,w} == ILDA (qx,qy,qz,qw) quaternion order
Matrix3 quatToRotMat(const Quaternion& q) {
    const double qx = q.x();
    const double qy = q.y();
    const double qz = q.z();
    const double qw = q.w();

    Matrix3 R{Matrix3::Zero()};
    R(0, 0) = 1.0 - 2.0 * (qy * qy + qz * qz);
    R(0, 1) = 2.0 * (qx * qy - qw * qz);
    R(0, 2) = 2.0 * (qx * qz + qw * qy);
    R(1, 0) = 2.0 * (qx * qy + qw * qz);
    R(1, 1) = 1.0 - 2.0 * (qx * qx + qz * qz);
    R(1, 2) = 2.0 * (qy * qz - qw * qx);
    R(2, 0) = 2.0 * (qx * qz - qw * qy);
    R(2, 1) = 2.0 * (qy * qz + qw * qx);
    R(2, 2) = 1.0 - 2.0 * (qx * qx + qy * qy);
    return R;
}

}

CisResult findCIS(
    const LammpsParser::Frame& frame,
    const std::vector<int>& grain,
    const std::vector<int>& structureType,
    const GrainRef& gA,
    const GrainRef& gB,
    const std::vector<std::uint8_t>& selection,
    IldaPtmNeighborQuery& query,
    double cisTol
) {
    const std::size_t n = static_cast<std::size_t>(frame.natoms);

    CisResult result;
    result.cis.assign(n, 0);
    result.cisNbrs.assign(n, std::array<int, 2>{-1, -1});
    result.cisIdealVecs.assign(n, std::array<Vector3, 2>{Vector3(0.0, 0.0, 0.0), Vector3(0.0, 0.0, 0.0)});
    result.cisSpatialVecs.assign(n, std::array<Vector3, 2>{Vector3(0.0, 0.0, 0.0), Vector3(0.0, 0.0, 0.0)});
    result.color.assign(n, Vector3(0.0, 0.0, 0.0));

    if(n == 0) {
        return result;
    }

    const std::array<const GrainRef*, 2> grains{&gA, &gB};

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions) {
        spdlog::error("findCIS: failed to create position property");
        return result;
    }
    NearestNeighborFinder finder(26);
    if(!finder.prepare(positions.get(), frame.simulationCell)) {
        spdlog::error("findCIS: NearestNeighborFinder::prepare failed");
        return result;
    }

    std::vector<std::size_t> g1g2atoms;    g1g2atoms.reserve(n);
    for(std::size_t i = 0; i < n; ++i) {
        if(!selection[i]) {
            continue;
        }
        const bool inA = (structureType[i] == static_cast<int>(gA.structureType)) && (grain[i] == gA.grainId);
        const bool inB = (structureType[i] == static_cast<int>(gB.structureType)) && (grain[i] == gB.grainId);
        if(inA || inB) {
            g1g2atoms.push_back(i);
        }
    }

    std::vector<std::size_t> g1g2candidates;
    g1g2candidates.reserve(g1g2atoms.size());
    {
        NearestNeighborFinder::Query<32> nnq(finder);
        for(const std::size_t m : g1g2atoms) {
            nnq.findNeighbors(m);
            const auto& res = nnq.results();
            int otherCount = 0;
            for(int j = 0; j < res.size(); ++j) {
                const std::size_t ni = res[j].index;
                if(ni < n && structureType[ni] == OTHER) {
                    ++otherCount;
                }
            }
            if(otherCount > 0) {
                g1g2candidates.push_back(m);
            }
        }
    }

    std::array<std::vector<Vector3>, 2> idealVecs;
    std::array<std::vector<Vector3>, 2> rotIdealVecs;
    std::array<Matrix3, 2> grainRs{
        quatToRotMat(gA.tableOrientation),
        quatToRotMat(gB.tableOrientation)
    };
    std::array<double, 2> grainScales{
        scalePTMInteratomicDistance(gA.structureType, 1.0) * scalePTMVec(gA.structureType, 1.0),
        scalePTMInteratomicDistance(gB.structureType, 1.0) * scalePTMVec(gB.structureType, 1.0)
    };
    for(int grainInd = 0; grainInd < 2; ++grainInd) {
        const int targetGrainId = grains[grainInd]->grainId;
        const Quaternion tableOri = grains[grainInd]->tableOrientation;
        const Quaternion refOri = grains[grainInd]->orientation;
        for(std::size_t m = 0; m < n; ++m) {
            if(grain[m] != targetGrainId) {
                continue;
            }
            query.findNeighbors(m, &tableOri);
            const int spatialCount = query.count();
            std::vector<Vector3> spatialIdeal(static_cast<std::size_t>(spatialCount));
            for(int jj = 0; jj < spatialCount; ++jj) {
                spatialIdeal[static_cast<std::size_t>(jj)] = query[jj].idealVector;
            }
            query.findNeighbors(m, &refOri);
            const int refCount = query.count();
            const int sharedCount = std::min(spatialCount, refCount);
            for(int jj = 0; jj < sharedCount; ++jj) {
                idealVecs[grainInd].push_back(query[jj].idealVector);
                rotIdealVecs[grainInd].push_back(
                    grainRs[grainInd] * (spatialIdeal[static_cast<std::size_t>(jj)] * grainScales[grainInd])
                );
            }
            break;  // PARITY: only the first atom of the grain
        }
    }

    for(std::size_t i = 0; i < g1g2candidates.size(); ++i) {
        const std::size_t m = g1g2candidates[i];

        Quaternion refOrientation;
        int grainInd = -1;
        if(grain[m] == gA.grainId) {
            refOrientation = gA.orientation;
            grainInd = 0;
        } else if(grain[m] == gB.grainId) {
            refOrientation = gB.orientation;
            grainInd = 1;
        } else {
            continue;
        }

        query.findNeighbors(m, &refOrientation);

        const double atomicDistM = (m < query.atomStates().size())
            ? query.atomStates()[m].interatomicDistance : 0.0;

        for(int jj = 0; jj < query.count(); ++jj) {
            const IldaPtmNeighborQuery::Neighbor neigh = query[jj];

            if(neigh.index < 0 || static_cast<std::size_t>(neigh.index) >= n) {
                continue;
            }
            if(structureType[static_cast<std::size_t>(neigh.index)] != OTHER) {
                continue;
            }
            // PARITY: skip if already explored for this grain (test is >0, not >=0)
            if(result.cisNbrs[static_cast<std::size_t>(neigh.index)][static_cast<std::size_t>(grainInd)] > 0) {
                continue;
            }

            result.cisNbrs[static_cast<std::size_t>(neigh.index)][static_cast<std::size_t>(grainInd)] = static_cast<int>(m);
            result.cisIdealVecs[static_cast<std::size_t>(neigh.index)][static_cast<std::size_t>(grainInd)] = -neigh.idealVector;
            result.cisSpatialVecs[static_cast<std::size_t>(neigh.index)][static_cast<std::size_t>(grainInd)] = -neigh.delta;

            if(cisTol > 0.0) {
                NearestNeighborFinder::Query<32> nnq2(finder);
                nnq2.findNeighbors(static_cast<std::size_t>(neigh.index));
                const auto& res2 = nnq2.results();
                for(int a = 0; a < res2.size(); ++a) {
                    const std::size_t neigh2Index = res2[a].index;
                    const Vector3 neigh2Delta = res2[a].delta;
                    if(neigh2Index >= n || structureType[neigh2Index] != OTHER) {
                        continue;
                    }
                    if(result.cisNbrs[neigh2Index][static_cast<std::size_t>(grainInd)] > 0) {
                        continue;
                    }

                    for(std::size_t k = 0; k < idealVecs[grainInd].size(); ++k) {
                        const Vector3 nbrVec = atomicDistM * rotIdealVecs[grainInd][k];
                        if((neigh2Delta - nbrVec).length() < cisTol) {
                            result.cisNbrs[neigh2Index][static_cast<std::size_t>(grainInd)] = static_cast<int>(m);
                            result.cisIdealVecs[neigh2Index][static_cast<std::size_t>(grainInd)] =
                                -neigh.idealVector - idealVecs[grainInd][k];
                            result.cisSpatialVecs[neigh2Index][static_cast<std::size_t>(grainInd)] =
                                -neigh.delta - neigh2Delta;

                            // PARITY: only the first eligible neigh3 (unconditional break)
                            NearestNeighborFinder::Query<32> nnq3(finder);
                            nnq3.findNeighbors(neigh2Index);
                            const auto& res3 = nnq3.results();
                            for(int b = 0; b < res3.size(); ++b) {
                                const std::size_t neigh3Index = res3[b].index;
                                const Vector3 neigh3Delta = res3[b].delta;
                                if(neigh3Index >= n || structureType[neigh3Index] != OTHER) {
                                    continue;
                                }
                                if(result.cisNbrs[neigh3Index][static_cast<std::size_t>(grainInd)] > 0) {
                                    continue;
                                }
                                for(std::size_t l = 0; l < idealVecs[grainInd].size(); ++l) {
                                    const Vector3 nbrVec3 = atomicDistM * rotIdealVecs[grainInd][l];
                                    if((neigh3Delta - nbrVec3).length() < cisTol) {
                                        result.cisNbrs[neigh3Index][static_cast<std::size_t>(grainInd)] = static_cast<int>(m);
                                        result.cisIdealVecs[neigh3Index][static_cast<std::size_t>(grainInd)] =
                                            -neigh.idealVector - idealVecs[grainInd][k] - idealVecs[grainInd][l];
                                        result.cisSpatialVecs[neigh3Index][static_cast<std::size_t>(grainInd)] =
                                            -neigh.delta - neigh2Delta - neigh3Delta;
                                        break;  // PARITY: break the l loop
                                    }
                                }
                                break;  // PARITY: break the neigh3 loop
                            }
                        }
                    }
                }
            }
        }
    }

    // PARITY: CIS atom requires entries >=0 for BOTH grains
    for(std::size_t i = 0; i < n; ++i) {
        if(result.cisNbrs[i][0] >= 0 && result.cisNbrs[i][1] >= 0) {
            result.cis[i] = 1;
            result.color[i] = Vector3(1.0, 1.0, 0.0);
        }
    }

    long cisCount = 0;
    for(const auto v : result.cis) {
        cisCount += v;
    }
    spdlog::info("Found {} CIS atoms", cisCount);

    return result;
}

}
