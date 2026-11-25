#pragma once

#include <volt/core/run_log.h>
#include <volt/ilda_types.h>

#include <optional>

namespace Volt {

class PTM;
class IldaPtmNeighborQuery;
struct IldaOptions;

// The execution log of one ILDA run, and of its self-tests. Each pass hands it what it built, and the log writes what
// that holds and how long the pass took; the warnings a pass raises along the way leave its time alone. Scripts parse
// these lines, so their wording stays as it is.
class IldaLog : public RunLog {
public:
    using RunLog::dump;

    void dump(const GrainLabels& labels);
    void dump(const GrainTable& table);
    void dump(const ResolvedGrains& grains);
    // Warns when the options ask for neither a circuit nor lines: the run then only finds the CIS atoms.
    void dump(const IldaOptions& options);
    void dump(const ReferenceFrameAxes& axes);
    void dump(const UserReferenceFrame& frame);
    // PTM's neighbor lists: their time counts in the line of the template matching that follows.
    void dump(const PTM& ptm);
    void dump(const IldaPtmNeighborQuery& query);
    void dump(const InterfaceCandidates& candidates);
    void dump(const ReferenceLatticeVectors& vectors);
    void dump(const CisResult& cis);
    void dump(const InterfaceSkin& skin);
    void dump(const DeformationGradientFit& fit);
    void dump(const CircuitFailure& failure);
    void dump(const SingleCircuit& circuit);
    void dump(const CisSurface& surface);
    void dump(const SurfacePruning& pruning);
    void dump(const CisSurfaceMesh& surfaceMesh);
    void dump(const PeriodicFaces& faces);
    void dump(const CircuitTally& tally);
    void dump(const DefectLines& defects);
    void dump(const WrittenLineSegments& table);
    void dump(const WrittenSummary& table);
    void dump(const PtmTemplateRemapFailure& failure);
    void dump(const UntestedPtmScale& scale);
    void dump(IldaIssue issue);

private:
    // The PTM neighbor lists' time, until the template matching line reports it.
    std::optional<Stopwatch::Clock::duration> _ptmNeighbors;
    // A pass that logs lines along the way closes with its whole time: the sum of the times its lines took.
    Stopwatch::Clock::duration _pass{};
};

}
