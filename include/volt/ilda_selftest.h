#pragma once

#include <volt/ilda_types.h>
#include <volt/core/lammps_parser.h>
#include <volt/plugin/plugin_main.h>

#include <nlohmann/json.hpp>

#include <cstddef>
#include <optional>
#include <vector>

namespace Volt {

class IldaLog;
class IldaPtmNeighborQuery;
class PTM;

// The self-tests, which run one stage of the analysis on a dump with the inputs an oracle file gives, and write what
// the stage found as JSON for comparison with the oracle's own results.
namespace IldaSelftest {

using json = nlohmann::json;
using Plugin::OptsMap;

// Runs the self-test the arguments ask for, and gives its exit code; none when they ask for none.
std::optional<int> run(int argc, char* argv[]);

int runPtm(const OptsMap& opts);
int runCis(const OptsMap& opts);
int runCrfCompute(const OptsMap& opts);
int runCrfEstimate(const OptsMap& opts);
int runMesh(const OptsMap& opts);
int runCircuit(const OptsMap& opts);
int runDefects(const OptsMap& opts);

// The grain and structure the oracle gives each atom of the frame, matched by particle identifier.
struct OracleLabels {
    std::vector<int> grain;
    std::vector<int> structureType;
    // Whether the oracle lists the atom.
    std::vector<std::uint8_t> listed;
    std::size_t aligned = 0;
};

StructureType ovitoTypeToVolt(int type);
Matrix3 jsonToMatrix3RowMajor(const json& values);
json matrix3ToJsonRowMajor(const Matrix3& matrix);
int particleId(const LammpsParser::Frame& frame, std::size_t atom);

// The oracle file of a self-test; false, once it said why, when it cannot be read.
bool readOracle(const std::string& path, json& oracle);
OracleLabels alignToOracle(const LammpsParser::Frame& frame, const json& oracle);
// Grain "g1" or "g2" of the oracle, in the slot it takes in the analysis.
GrainRef oracleGrain(const json& oracle, const char* name, int slot, double latticeConstant);
// Matches the PTM templates of every atom of the frame; false, once it said why, when PTM cannot be set up.
bool matchPtm(PTM& ptm, IldaPtmNeighborQuery& query, const LammpsParser::Frame& frame, double rmsd, IldaLog& log);
// Prints the result, or writes it to the --out file and says so with the summary.
void writeResult(const json& out, const std::string& outPath, const std::string& summary);

}

}
