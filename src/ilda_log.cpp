#include <volt/ilda_log.h>

#include <volt/analysis/ptm.h>
#include <volt/core/type_counts.h>
#include <volt/ilda_crf.h>
#include <volt/ilda_mesh.h>
#include <volt/ilda_options.h>
#include <volt/ilda_ptm_neighbor_query.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace Volt {

void IldaLog::dump(const GrainLabels& labels) {
    spdlog::info("  Aligned {}/{} atoms to grain labels from {} in {}",
                 labels.aligned, labels.atoms, labels.path, phaseElapsed());
    if(labels.aligned != labels.atoms) {
        spdlog::warn("{} atoms in the frame were not present in the grain-atoms parquet "
                     "and stay unassigned (grain 0, OTHER)", labels.atoms - labels.aligned);
    }
}

void IldaLog::dump(const GrainTable& table) {
    spdlog::info("  Loaded {} grains from {} in {}", table.grains, table.path, phaseElapsed());
}

void IldaLog::dump(const ResolvedGrains& grains) {
    spdlog::info("Grain A: id={} type={} | grain B: id={} type={}",
                 grains.grainA.grainId, structureTypeName(grains.grainA.structureType),
                 grains.grainB.grainId, structureTypeName(grains.grainB.structureType));
    phaseMilliseconds();
}

void IldaLog::dump(const IldaOptions& options) {
    if(!options.single_circuit && !options.extract_lines) {
        spdlog::warn("Neither extract_lines nor single_circuit is enabled: only CIS atoms are detected, "
                     "no line defects are written");
    }
}

void IldaLog::dump(const ReferenceFrameAxes& axes) {
    if(std::abs(axes.xA.dot(axes.yA)) > 1e-6) {
        spdlog::error("computeF: xA and yA are not orthogonal!");
    }
    if(std::abs(axes.xB.dot(axes.yB)) > 1e-6) {
        spdlog::error("computeF: xB and yB are not orthogonal!");
    }
}

void IldaLog::dump(const UserReferenceFrame&) {
    spdlog::info("User specified coherent reference frame");
    phaseMilliseconds();
}

void IldaLog::dump(const PTM&) {
    _ptmNeighbors = phaseDuration();
}

void IldaLog::dump(const IldaPtmNeighborQuery& query) {
    const Stopwatch::Clock::duration matching = phaseDuration();
    const auto& states = query.atomStates();
    if(!states.empty()) {
        std::array<std::size_t, NUM_STRUCTURE_TYPES> totals{};
        for(const auto& state : states) {
            if(static_cast<std::size_t>(state.structureType) < totals.size()) {
                ++totals[static_cast<std::size_t>(state.structureType)];
            }
        }
        TypeCounts counts{"  Matched PTM templates in " + formatDuration(matching), {}, states.size()};
        for(std::size_t type = 1; type < totals.size(); ++type) {
            if(totals[type] > 0) {
                counts.counts.emplace_back(structureTypeName(static_cast<int>(type)), totals[type]);
            }
        }
        counts.counts.emplace_back("Other", totals[OTHER]);
        dump(counts);
    }
    if(_ptmNeighbors) {
        spdlog::info("Polyhedral template matching done in {}", formatDuration(*_ptmNeighbors + matching));
        _ptmNeighbors.reset();
    }
}

void IldaLog::dump(const InterfaceCandidates& candidates) {
    _pass = phaseDuration();
    spdlog::info("  Found {} interface candidates among {} grain A and {} grain B atoms in {}",
                 candidates.candidates, candidates.grainAAtoms, candidates.grainBAtoms, formatDuration(_pass));
}

void IldaLog::dump(const ReferenceLatticeVectors& vectors) {
    if(vectors.cisTol > 0.0 && !vectors.found) {
        spdlog::warn("The first atom of grain {} has no PTM neighbors to take reference lattice vectors "
                     "from: cis_tol has no effect for this grain", vectors.grainId);
    }
}

void IldaLog::dump(const CisResult& cis) {
    _pass += phaseDuration();
    const auto cisCount = std::count(cis.cis.begin(), cis.cis.end(), static_cast<std::uint8_t>(1));
    spdlog::info("Found {} CIS atoms in {}", cisCount, formatDuration(_pass));
}

void IldaLog::dump(const InterfaceSkin& skin) {
    _pass = phaseDuration();
    spdlog::info("  Collected {} atoms within distF of the interface in {}", skin.atoms, formatDuration(_pass));
}

void IldaLog::dump(const DeformationGradientFit& fit) {
    _pass += phaseDuration();
    if(fit.atomsA <= 0) {
        spdlog::warn("estimateF: no atoms found to estimate FA");
    }
    if(fit.atomsB <= 0) {
        spdlog::warn("estimateF: no atoms found to estimate FB");
    }
    spdlog::info("  Used {} and {} atoms to estimate FA and FB, respectively.", fit.atomsA, fit.atomsB);

    Matrix3 VA, RA, VB, RB;
    polarDecomp(fit.grainA.F, VA, RA);
    polarDecomp(fit.grainB.F, VB, RB);
    spdlog::info("Estimated coherent reference frame (with polar decomposition) in {}: "
                 "trace(VA)={:.6f}, trace(VB)={:.6f}",
                 formatDuration(_pass), VA.trace(), VB.trace());
}

void IldaLog::dump(const CircuitFailure& failure) {
    switch(failure.reason) {
        case CircuitFailure::Reason::InvalidIndices:
            spdlog::warn("burgersCircuit: invalid indices {} {}", failure.startIndex, failure.endIndex);
            break;
        case CircuitFailure::Reason::NotCoincidentSites:
            spdlog::warn("Wrong start or end indices: not a co-incident site {} {}", failure.startIndex, failure.endIndex);
            break;
        case CircuitFailure::Reason::SameAtom:
            spdlog::warn("Wrong start or end indices, both are same: {}", failure.startIndex);
            break;
        case CircuitFailure::Reason::NotClosed:
            spdlog::warn("Warning: Not a closed circuit! {} {} spatialSum=({}, {}, {})",
                         failure.startIndex, failure.endIndex,
                         failure.spatialSum.x(), failure.spatialSum.y(), failure.spatialSum.z());
            break;
    }
}

void IldaLog::dump(const SingleCircuit& circuit) {
    switch(circuit.outcome) {
        case SingleCircuit::Outcome::InvalidAtom:
            spdlog::error("Error: invalid single circuit atom identifier");
            break;
        case SingleCircuit::Outcome::NotCoincidentSites:
            spdlog::error("Error: Cannot construct specified Burgers circuit because "
                          "one or both atoms is not a CIS atom");
            break;
        case SingleCircuit::Outcome::Traced:
            if(circuit.circuit) {
                const CircuitResult& result = *circuit.circuit;
                spdlog::info("Info for Burgers circuit using atoms {} and {}:", circuit.atom1, circuit.atom2);
                spdlog::info("b = ({}, {}, {})", result.latticeSum.x(), result.latticeSum.y(), result.latticeSum.z());
                spdlog::info("Ideal step height A = {:.3f}", result.hA);
                spdlog::info("Ideal step height B = {:.3f}", result.hB);
            } else {
                spdlog::warn("single_circuit produced no closed circuit");
            }
            break;
    }
    phaseMilliseconds();
}

void IldaLog::dump(const CisSurface& surface) {
    _pass = phaseDuration();
    spdlog::info("  Surface constructed in {}: {} faces, {} vertices",
                 formatDuration(_pass), surface.mesh.faceCount(), surface.mesh.vertexCount());
}

void IldaLog::dump(const SurfacePruning& pruning) {
    const Stopwatch::Clock::duration pruned = phaseDuration();
    _pass += pruned;
    if(pruning.prunedAtoms > ILDA_PRUNE_MAX) {
        spdlog::warn("Number of surface mesh pruning events exceeded {} - pruning was aborted", ILDA_PRUNE_MAX);
    }
    spdlog::info("  Pruned {} non-CIS atoms from the surface mesh in {}", pruning.prunedAtoms, formatDuration(pruned));
}

void IldaLog::dump(const CisSurfaceMesh& surfaceMesh) {
    _pass += phaseDuration();
    if(!surfaceMesh.halfEdgesConnected) {
        spdlog::warn("Unable to connect opposite half edges of the pruned surface mesh.");
    }
    spdlog::info("CIS surface mesh ready in {}: {} faces", formatDuration(_pass), surfaceMesh.mesh.faceCount());
}

void IldaLog::dump(const PeriodicFaces& faces) {
    if(faces.skipped > 0) {
        spdlog::warn("{} mesh faces are too large to use periodic boundary "
                     "conditions, these have been skipped by ILDA", faces.skipped);
    }
}

void IldaLog::dump(const CircuitTally& tally) {
    _pass = phaseDuration();
    spdlog::info("  Traced {} Burgers circuits in {} ({} with |b| > btol, {} failed)",
                 tally.traced, formatDuration(_pass), tally.aboveTolerance, tally.failed);
}

void IldaLog::dump(const DefectLines& defects) {
    _pass += phaseDuration();
    if(defects.stepHeightsUnknown) {
        spdlog::warn("Terrace plane normal in coherent reference frame, n, not "
                     "specified, cannot determine ideal step heights");
    }
    spdlog::info("Found {} unique disconnection modes and {} line segments in {} (total line length {:.3f})",
                 defects.modes.size(), defects.segments.size(), formatDuration(_pass), defects.totalLength);
}

void IldaLog::dump(const WrittenLineSegments& table) {
    spdlog::info("Wrote ILDA line segments -> {}", table.path);
    phaseMilliseconds();
}

void IldaLog::dump(const WrittenSummary& table) {
    std::size_t variants = 0;
    for(const DisconnectionMode& mode : table.modes) {
        variants += mode.variants.size();
    }
    spdlog::info("Wrote ILDA summary ({} modes, {} variants, L={:.3f}) -> {}",
                 table.modes.size(), variants, table.totalLength, table.path);
    phaseMilliseconds();
}

void IldaLog::dump(const PtmTemplateRemapFailure& failure) {
    spdlog::warn("IldaPtmNeighborQuery: ptm_remap_template failed for atom {}", failure.atom);
}

void IldaLog::dump(const UntestedPtmScale& scale) {
    spdlog::warn("scalePTMVec: untested structure type {}", static_cast<int>(scale.type));
}

void IldaLog::dump(IldaIssue issue) {
    switch(issue) {
        case IldaIssue::CisPositions:
            spdlog::error("findCIS: failed to create position property");
            break;
        case IldaIssue::CisNeighborLists:
            spdlog::error("findCIS: NearestNeighborFinder::prepare failed");
            break;
        case IldaIssue::InterfaceSkinPositions:
            spdlog::error("estimateF: failed to create position property");
            break;
        case IldaIssue::InterfaceSkinNeighborLists:
            spdlog::error("estimateF: CutoffNeighborFinder::prepare failed");
            break;
        case IldaIssue::SurfaceInput:
            spdlog::warn("buildCisSurfaceMesh: empty frame or non-positive Rsphere");
            break;
        case IldaIssue::SurfacePositions:
            spdlog::error("buildCisSurfaceMesh: failed to create position property");
            break;
        case IldaIssue::SurfaceConstruction:
            spdlog::error("buildCisSurfaceMesh: ManifoldConstructionHelper::construct failed");
            break;
        case IldaIssue::MixedCisMeshEdge:
            spdlog::warn("Mesh edge connects CIS atom to non-CIS atom, "
                         "Rsphere probably needs to be increased");
            break;
    }
}

}
