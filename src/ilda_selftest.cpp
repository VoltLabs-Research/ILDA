#include <volt/ilda_selftest.h>

#include <volt/analysis/ptm.h>
#include <volt/cli/common.h>
#include <volt/core/progress.h>
#include <volt/ilda_circuit.h>
#include <volt/ilda_cis.h>
#include <volt/ilda_crf.h>
#include <volt/ilda_defects.h>
#include <volt/ilda_json_exporter.h>
#include <volt/ilda_log.h>
#include <volt/ilda_mesh.h>
#include <volt/ilda_ptm_neighbor_query.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>

namespace Volt::IldaSelftest {

using namespace Volt::CLI;

std::optional<int> run(int argc, char* argv[]) {
    std::string filename;
    std::string outputBase;
    const auto opts = parseArgs(argc, argv, filename, outputBase);
    static const std::pair<const char*, int (*)(const OptsMap&)> selftests[] = {
        {"--selftest-ptm", runPtm},
        {"--selftest-cis", runCis},
        {"--selftest-crf-compute", runCrfCompute},
        {"--selftest-crf-estimate", runCrfEstimate},
        {"--selftest-mesh", runMesh},
        {"--selftest-circuit", runCircuit},
        {"--selftest-defects", runDefects},
    };
    for(const auto& [flag, selftest] : selftests) {
        if(hasOption(opts, flag)) {
            initLogging("ilda-selftest");
            return selftest(opts);
        }
    }
    return std::nullopt;
}

StructureType ovitoTypeToVolt(int type) {
    switch(type) {
        case 1: return FCC;
        case 2: return HCP;
        case 3: return BCC;
        case 4: return ICO;
        case 5: return SC;
        case 6: return CUBIC_DIAMOND;
        case 7: return HEX_DIAMOND;
        case 8: return GRAPHENE;
        default: return OTHER;
    }
}

Matrix3 jsonToMatrix3RowMajor(const json& values) {
    Matrix3 m{Matrix3::Zero()};
    for(int r = 0; r < 3; ++r)
        for(int c = 0; c < 3; ++c)
            m(static_cast<Matrix3::size_type>(r), static_cast<Matrix3::size_type>(c)) =
                values[r * 3 + c].get<double>();
    return m;
}

json matrix3ToJsonRowMajor(const Matrix3& matrix) {
    json values = json::array();
    for(int r = 0; r < 3; ++r)
        for(int c = 0; c < 3; ++c)
            values.push_back(matrix(static_cast<Matrix3::size_type>(r), static_cast<Matrix3::size_type>(c)));
    return values;
}

int particleId(const LammpsParser::Frame& frame, std::size_t atom) {
    return (atom < frame.ids.size()) ? frame.ids[atom] : static_cast<int>(atom);
}

bool readOracle(const std::string& path, json& oracle) {
    std::ifstream oracleStream(path);
    if(!oracleStream.is_open()) {
        std::cerr << "Cannot open oracle: " << path << "\n";
        return false;
    }
    oracle = json::parse(oracleStream);
    return true;
}

OracleLabels alignToOracle(const LammpsParser::Frame& frame, const json& oracle) {
    const std::size_t n = static_cast<std::size_t>(frame.natoms);
    const auto& oraclePids = oracle.at("particle_id");
    const auto& oracleGrain = oracle.at("grain");
    const auto& oracleStype = oracle.at("structure_type");

    std::unordered_map<int, std::size_t> oraclePidToRow;
    oraclePidToRow.reserve(oraclePids.size());
    for(std::size_t r = 0; r < oraclePids.size(); ++r) {
        oraclePidToRow.emplace(oraclePids[r].get<int>(), r);
    }

    OracleLabels labels;
    labels.grain.assign(n, 0);
    labels.structureType.assign(n, static_cast<int>(OTHER));
    labels.listed.assign(n, 0);
    for(std::size_t i = 0; i < n; ++i) {
        auto it = oraclePidToRow.find(particleId(frame, i));
        if(it == oraclePidToRow.end()) continue;
        const std::size_t r = it->second;
        labels.grain[i] = oracleGrain[r].get<int>();
        labels.structureType[i] = static_cast<int>(ovitoTypeToVolt(oracleStype[r].get<int>()));
        labels.listed[i] = 1;
        ++labels.aligned;
    }
    std::cerr << "Aligned " << labels.aligned << "/" << n << " atoms to the oracle by particle id\n";
    return labels;
}

GrainRef oracleGrain(const json& oracle, const char* name, int slot, double latticeConstant) {
    const std::string prefix(name);
    GrainRef grain;
    grain.grainId = oracle.at(prefix + "_id").get<int>();
    grain.grainIndex = slot;
    grain.a = latticeConstant;
    grain.structureType = ovitoTypeToVolt(oracle.at(prefix + "_s").get<int>());
    const auto& q = oracle.at(prefix + "_orientation_xyzw");
    grain.tableOrientation = Quaternion(q[0].get<double>(), q[1].get<double>(),
                                        q[2].get<double>(), q[3].get<double>());
    grain.orientation = grain.tableOrientation;
    return grain;
}

bool matchPtm(PTM& ptm, IldaPtmNeighborQuery& query, const LammpsParser::Frame& frame, double rmsd, IldaLog& log) {
    ptm.setRmsdCutoff(rmsd);
    ptm.setCalculateDefGradient(false);
    if(!ptm.prepare(frame.positions.data(), frame.positions.size(), frame.simulationCell)) {
        std::cerr << "PTM prepare failed\n";
        return false;
    }
    query.computeAtomStates();
    log.dump(query);
    return true;
}

void writeResult(const json& out, const std::string& outPath, const std::string& summary) {
    if(outPath.empty()) {
        std::cout << out.dump(2) << "\n";
    } else {
        std::ofstream os(outPath);
        os << out.dump(2) << "\n";
        std::cerr << "Wrote " << outPath << " (" << summary << ")\n";
    }
}

int runPtm(const OptsMap& opts) {
    const Stopwatch timer;
    const std::string dump = getString(opts, "--selftest-ptm");
    const std::string oraclePath = getString(opts, "--oracle");
    const std::string outPath = getString(opts, "--out", "");
    if(dump.empty() || oraclePath.empty()) {
        std::cerr << "Usage: ilda --selftest-ptm <dump> --oracle <oracle.json> [--out <cpp.json>]\n";
        return 1;
    }
    json oracle;
    if(!readOracle(oraclePath, oracle)) return 1;
    const double rmsd = oracle.value("rmsd", 0.077);
    const auto qextArr = oracle.at("qext_xyzw");
    const Quaternion qext(qextArr[0].get<double>(), qextArr[1].get<double>(),
                          qextArr[2].get<double>(), qextArr[3].get<double>());
    std::cerr << "PTM self-test: rmsd=" << rmsd << "\n";

    IldaLog log;
    LammpsParser::Frame frame;
    if(!parseFrame(dump, frame)) return 1;

    PTM ptm;
    IldaPtmNeighborQuery query(ptm, log);
    if(!matchPtm(ptm, query, frame, rmsd, log)) return 1;

    auto dumpCase = [&](std::size_t idx, const Quaternion* qref) {
        query.findNeighbors(idx, qref);
        json rec;
        rec["index"] = static_cast<int>(idx);
        rec["count"] = query.count();
        const Quaternion o = query.orientation();
        rec["orientation"] = {o.x(), o.y(), o.z(), o.w()};
        json neighbors = json::array();
        for(int j = 0; j < query.count(); ++j) {
            const auto& n = query[j];
            neighbors.push_back({
                {"slot", j},
                {"index", n.index},
                {"delta", {n.delta.x(), n.delta.y(), n.delta.z()}},
                {"ideal_vector", {n.idealVector.x(), n.idealVector.y(), n.idealVector.z()}}
            });
        }
        rec["neighbors"] = neighbors;
        return rec;
    };

    json out;
    out["dump"] = dump;
    out["rmsd"] = rmsd;
    json atoms = json::array();
    for(const auto& a : oracle.at("atoms")) {
        const std::size_t idx = a.at("none_case").at("index").get<std::size_t>();
        json entry;
        entry["type_name"] = a.value("type_name", "");
        entry["structure_type"] = static_cast<int>(query.atomStates()[idx].structureType);
        entry["none_case"] = dumpCase(idx, nullptr);
        entry["qref_external_case"] = dumpCase(idx, &qext);
        const Quaternion qself = query.atomStates()[idx].orientation;
        entry["qref_self_case"] = dumpCase(idx, &qself);
        atoms.push_back(entry);
    }
    out["atoms"] = atoms;

    std::ostringstream summary;
    summary << atoms.size() << " atoms";
    writeResult(out, outPath, summary.str());
    std::cerr << "PTM self-test completed in " << timer.elapsed() << "\n";
    return 0;
}

int runCis(const OptsMap& opts) {
    const Stopwatch timer;
    const std::string dump = getString(opts, "--selftest-cis");
    const std::string oraclePath = getString(opts, "--oracle");
    const std::string outPath = getString(opts, "--out", "");
    if(dump.empty() || oraclePath.empty()) {
        std::cerr << "Usage: ilda --selftest-cis <dump> --oracle <cis_oracle.json> [--out <cis_cpp.json>]\n";
        return 1;
    }
    json oracle;
    if(!readOracle(oraclePath, oracle)) return 1;
    const double rmsd = oracle.value("rmsd", 0.077);
    const double cisTol = oracle.value("cis_tol", 0.0);
    const double aA = oracle.value("aA", 1.0);
    const double aB = oracle.value("aB", 1.0);
    std::cerr << "CIS self-test: rmsd=" << rmsd << ", cis_tol=" << cisTol << ", aA=" << aA << ", aB=" << aB << "\n";

    IldaLog log;
    LammpsParser::Frame frame;
    if(!parseFrame(dump, frame)) return 1;
    const std::size_t n = static_cast<std::size_t>(frame.natoms);

    const OracleLabels labels = alignToOracle(frame, oracle);
    const std::vector<std::uint8_t> selection(n, 1);
    const GrainRef gA = oracleGrain(oracle, "g1", 0, aA);
    const GrainRef gB = oracleGrain(oracle, "g2", 1, aB);

    PTM ptm;
    IldaPtmNeighborQuery query(ptm, log);
    if(!matchPtm(ptm, query, frame, rmsd, log)) return 1;

    CisResult cis = findCIS(frame, labels.grain, labels.structureType, gA, gB, selection, query, cisTol, log);

    long cisCount = 0;
    std::vector<int> cisPids;
    for(std::size_t i = 0; i < n; ++i) {
        if(cis.cis[i]) {
            ++cisCount;
            cisPids.push_back(particleId(frame, i));
        }
    }
    std::sort(cisPids.begin(), cisPids.end());

    json out;
    out["dump"] = dump;
    out["rmsd"] = rmsd;
    out["cis_tol"] = cisTol;
    out["natoms"] = static_cast<int>(n);
    out["sum_cis"] = static_cast<int>(cisCount);
    out["aligned"] = static_cast<int>(labels.aligned);
    out["cis_particle_ids"] = cisPids;

    std::ostringstream summary;
    summary << "sum_cis=" << cisCount;
    writeResult(out, outPath, summary.str());
    std::cerr << "CIS self-test completed in " << timer.elapsed() << "\n";
    return 0;
}

int runCrfCompute(const OptsMap& opts) {
    const Stopwatch timer;
    const std::string oraclePath = getString(opts, "--oracle");
    const std::string outPath = getString(opts, "--out", "");
    if(oraclePath.empty()) {
        std::cerr << "Usage: ilda --selftest-crf-compute --oracle <oracle.json> [--out <cpp.json>]\n";
        return 1;
    }
    json oracle;
    if(!readOracle(oraclePath, oracle)) return 1;

    auto readVec = [&](const char* key) {
        const auto& a = oracle.at(key);
        return Vector3(a[0].get<double>(), a[1].get<double>(), a[2].get<double>());
    };
    const Vector3 xA = readVec("xA");
    const Vector3 yA = readVec("yA");
    const Vector3 xB = readVec("xB");
    const Vector3 yB = readVec("yB");
    const Matrix3 EcohA = jsonToMatrix3RowMajor(oracle.at("EcohA"));
    const Matrix3 EcohB = jsonToMatrix3RowMajor(oracle.at("EcohB"));
    const double coveraA = oracle.value("coveraA", 0.0);
    const double coveraB = oracle.value("coveraB", 0.0);
    std::cerr << "CRF compute self-test: coveraA=" << coveraA << ", coveraB=" << coveraB << "\n";

    IldaLog log;
    GrainRef gA, gB;
    gA.grainId = 1;
    gB.grainId = 2;
    log.dump(ReferenceFrameAxes{xA, yA, xB, yB});
    computeF(xA, yA, xB, yB, EcohA, EcohB, coveraA, coveraB, gA, gB);

    json out;
    out["mode"] = "compute";
    out["FA"] = matrix3ToJsonRowMajor(gA.F);
    out["FB"] = matrix3ToJsonRowMajor(gB.F);
    out["g1_orientation_xyzw"] = {gA.orientation.x(), gA.orientation.y(), gA.orientation.z(), gA.orientation.w()};
    out["g2_orientation_xyzw"] = {gB.orientation.x(), gB.orientation.y(), gB.orientation.z(), gB.orientation.w()};

    writeResult(out, outPath, "crf-compute");
    std::cerr << "CRF compute self-test completed in " << timer.elapsed() << "\n";
    return 0;
}

int runCrfEstimate(const OptsMap& opts) {
    const Stopwatch timer;
    const std::string dump = getString(opts, "--selftest-crf-estimate");
    const std::string oraclePath = getString(opts, "--oracle");
    const std::string outPath = getString(opts, "--out", "");
    if(dump.empty() || oraclePath.empty()) {
        std::cerr << "Usage: ilda --selftest-crf-estimate <dump> --oracle <oracle.json> [--out <cpp.json>]\n";
        return 1;
    }
    json oracle;
    if(!readOracle(oraclePath, oracle)) return 1;
    const double rmsd = oracle.value("rmsd", 0.077);
    const double cisTol = oracle.value("cis_tol", 0.0);
    const double aA = oracle.value("aA", 1.0);
    const double aB = oracle.value("aB", 1.0);
    const double distF = oracle.value("distF", 10.0);
    std::cerr << "CRF estimate self-test: rmsd=" << rmsd << ", cis_tol=" << cisTol << ", aA=" << aA
              << ", aB=" << aB << ", distF=" << distF << "\n";

    IldaLog log;
    LammpsParser::Frame frame;
    if(!parseFrame(dump, frame)) return 1;
    const std::size_t n = static_cast<std::size_t>(frame.natoms);

    const OracleLabels labels = alignToOracle(frame, oracle);
    const std::vector<std::uint8_t> selection(n, 1);
    GrainRef gA = oracleGrain(oracle, "g1", 0, aA);
    GrainRef gB = oracleGrain(oracle, "g2", 1, aB);

    PTM ptm;
    IldaPtmNeighborQuery query(ptm, log);
    if(!matchPtm(ptm, query, frame, rmsd, log)) return 1;

    CisResult cis = findCIS(frame, labels.grain, labels.structureType, gA, gB, selection, query, cisTol, log);
    estimateF(frame, cis, labels.grain, labels.structureType, selection, distF, query, gA, gB, log);

    json out;
    out["mode"] = "estimate";
    out["dump"] = dump;
    out["aligned"] = static_cast<int>(labels.aligned);
    out["FA"] = matrix3ToJsonRowMajor(gA.F);
    out["FB"] = matrix3ToJsonRowMajor(gB.F);

    writeResult(out, outPath, "crf-estimate");
    std::cerr << "CRF estimate self-test completed in " << timer.elapsed() << "\n";
    return 0;
}

int runMesh(const OptsMap& opts) {
    const Stopwatch timer;
    const std::string dump = getString(opts, "--selftest-mesh");
    const std::string oraclePath = getString(opts, "--oracle");
    const std::string outPath = getString(opts, "--out", "");
    if(dump.empty() || oraclePath.empty()) {
        std::cerr << "Usage: ilda --selftest-mesh <dump> --oracle <mesh_oracle.json> [--out <cpp.json>]\n";
        return 1;
    }
    json oracle;
    if(!readOracle(oraclePath, oracle)) return 1;
    const double Rsphere = oracle.value("Rsphere", 15.0);
    const int grainB = oracle.at("grain_b").get<int>();
    const int g2sVolt = static_cast<int>(ovitoTypeToVolt(oracle.at("g2_s").get<int>()));
    std::cerr << "Mesh self-test: Rsphere=" << Rsphere << ", grain B=" << grainB
              << " (" << structureTypeName(g2sVolt) << ")\n";

    IldaLog log;
    LammpsParser::Frame frame;
    if(!parseFrame(dump, frame)) return 1;
    const std::size_t n = static_cast<std::size_t>(frame.natoms);

    std::set<int> cisPidSet;
    for(const auto& pid : oracle.at("cis_particle_ids")) {
        cisPidSet.insert(pid.get<int>());
    }
    const OracleLabels labels = alignToOracle(frame, oracle);
    std::vector<std::uint8_t> cis(n, 0);
    for(std::size_t i = 0; i < n; ++i) {
        cis[i] = (labels.listed[i] && cisPidSet.count(particleId(frame, i))) ? 1 : 0;
    }

    long sumCis = 0;
    for(const auto v : cis) sumCis += v;

    IldaMesh mesh;
    const bool ok = buildCisSurfaceMesh(frame.simulationCell, frame, cis, labels.grain,
                                        labels.structureType, grainB, g2sVolt, Rsphere, mesh, log);

    const std::size_t faceCount = mesh.faceCount();
    std::set<int> usedVerts;
    for(IldaMesh::Face* f : mesh.faces()) {
        IldaMesh::Edge* e = f->edges();
        if(!e) continue;
        usedVerts.insert(e->vertex1()->index());
        usedVerts.insert(e->vertex2()->index());
        usedVerts.insert(e->nextFaceEdge()->vertex2()->index());
    }
    int usedCis = 0;
    int usedNonCis = 0;
    std::set<int> usedAtoms;
    for(const int vidx : usedVerts) {
        IldaMesh::Vertex* v = mesh.vertex(vidx);
        if(v->isCIS) ++usedCis; else ++usedNonCis;
        if(v->atomIndex >= 0) usedAtoms.insert(v->atomIndex);
    }
    std::vector<int> usedPids;
    for(const int atom : usedAtoms) {
        usedPids.push_back(particleId(frame, static_cast<std::size_t>(atom)));
    }
    std::sort(usedPids.begin(), usedPids.end());

    json out;
    out["dump"] = dump;
    out["Rsphere"] = Rsphere;
    out["ok"] = ok;
    out["aligned"] = static_cast<int>(labels.aligned);
    out["sum_cis"] = static_cast<int>(sumCis);
    out["final_faces"] = static_cast<int>(faceCount);
    out["final_verts"] = static_cast<int>(mesh.vertexCount());
    out["final_used_verts"] = static_cast<int>(usedVerts.size());
    out["final_used_cis"] = usedCis;
    out["final_used_noncis"] = usedNonCis;
    out["final_used_pids"] = usedPids;

    std::ostringstream summary;
    summary << "final_faces=" << faceCount << " used_cis=" << usedCis << " used_noncis=" << usedNonCis;
    writeResult(out, outPath, summary.str());
    std::cerr << "Mesh self-test completed in " << timer.elapsed() << "\n";
    return 0;
}

int runCircuit(const OptsMap& opts) {
    const Stopwatch timer;
    const std::string dump = getString(opts, "--selftest-circuit");
    const std::string oraclePath = getString(opts, "--oracle");
    const std::string outPath = getString(opts, "--out", "");
    if(dump.empty() || oraclePath.empty()) {
        std::cerr << "Usage: ilda --selftest-circuit <dump> --oracle <oracle.json> [--out <cpp.json>]\n";
        return 1;
    }
    json oracle;
    if(!readOracle(oraclePath, oracle)) return 1;
    const double rmsd = oracle.value("rmsd", 0.077);
    const double cisTol = oracle.value("cis_tol", 0.0);
    const double aA = oracle.value("aA", 1.0);
    const double aB = oracle.value("aB", 1.0);
    std::cerr << "Circuit self-test: rmsd=" << rmsd << ", cis_tol=" << cisTol << ", aA=" << aA << ", aB=" << aB << "\n";

    IldaLog log;
    LammpsParser::Frame frame;
    if(!parseFrame(dump, frame)) return 1;
    const std::size_t n = static_cast<std::size_t>(frame.natoms);

    const OracleLabels labels = alignToOracle(frame, oracle);
    const std::vector<std::uint8_t> selection(n, 1);

    std::unordered_map<int, int> pidToIndex;
    pidToIndex.reserve(n);
    for(std::size_t i = 0; i < n; ++i) {
        pidToIndex.emplace(particleId(frame, i), static_cast<int>(i));
    }

    GrainRef gA = oracleGrain(oracle, "g1", 0, aA);
    gA.F = jsonToMatrix3RowMajor(oracle.at("FA"));
    GrainRef gB = oracleGrain(oracle, "g2", 1, aB);
    gB.F = jsonToMatrix3RowMajor(oracle.at("FB"));

    const auto& nArr = oracle.at("n");
    const Vector3 nNormal(nArr[0].get<double>(), nArr[1].get<double>(), nArr[2].get<double>());

    PTM ptm;
    IldaPtmNeighborQuery query(ptm, log);
    if(!matchPtm(ptm, query, frame, rmsd, log)) return 1;

    CisResult cis = findCIS(frame, labels.grain, labels.structureType, gA, gB, selection, query, cisTol, log);

    CircuitContext ctx{
        frame, frame.simulationCell, frame.positions,
        labels.grain, labels.structureType, selection, cis, gA, gB,
        nNormal, query
    };

    json circuitsOut = json::array();
    for(const auto& c : oracle.at("circuits")) {
        const int pid1 = c.at("circuitAtom1").get<int>();
        const int pid2 = c.at("circuitAtom2").get<int>();
        auto i1 = pidToIndex.find(pid1);
        auto i2 = pidToIndex.find(pid2);
        json rec;
        rec["circuitAtom1"] = pid1;
        rec["circuitAtom2"] = pid2;
        if(i1 == pidToIndex.end() || i2 == pidToIndex.end()) {
            rec["ok"] = false;
            rec["error"] = "particle id not found in frame";
            circuitsOut.push_back(rec);
            continue;
        }
        auto result = burgersCircuit(i1->second, i2->second, ctx, log);
        if(result) {
            rec["ok"] = true;
            rec["b"] = {result->latticeSum.x(), result->latticeSum.y(), result->latticeSum.z()};
            rec["hA"] = result->hA;
            rec["hB"] = result->hB;
        } else {
            rec["ok"] = false;
            rec["error"] = "no closed circuit";
        }
        circuitsOut.push_back(rec);
    }

    long sumCis = 0;
    for(const auto v : cis.cis) sumCis += v;

    json out;
    out["dump"] = dump;
    out["aligned"] = static_cast<int>(labels.aligned);
    out["sum_cis"] = static_cast<int>(sumCis);
    out["circuits"] = circuitsOut;

    std::ostringstream summary;
    summary << "circuit; sum_cis=" << sumCis << ", " << circuitsOut.size() << " circuits";
    writeResult(out, outPath, summary.str());
    std::cerr << "Circuit self-test completed in " << timer.elapsed() << "\n";
    return 0;
}

int runDefects(const OptsMap& opts) {
    const Stopwatch timer;
    const std::string dump = getString(opts, "--selftest-defects");
    const std::string oraclePath = getString(opts, "--oracle");
    const std::string outPath = getString(opts, "--out", "");
    if(dump.empty() || oraclePath.empty()) {
        std::cerr << "Usage: ilda --selftest-defects <dump> --oracle <circuit_oracle.json> [--out <cpp.json>]\n";
        return 1;
    }
    json oracle;
    if(!readOracle(oraclePath, oracle)) return 1;
    const double rmsd = oracle.value("rmsd", 0.077);
    const double cisTol = oracle.value("cis_tol", 0.0);
    const double aA = oracle.value("aA", 1.0);
    const double aB = oracle.value("aB", 1.0);
    const double btol = oracle.value("btol", 0.01);
    const double htol = oracle.value("htol", 0.5);
    const double angtolDeg = oracle.value("angtol", 5.0);
    const double Rsphere = oracle.value("Rsphere", 10.0);
    std::cerr << "Defects self-test: rmsd=" << rmsd << ", cis_tol=" << cisTol << ", aA=" << aA << ", aB=" << aB
              << ", btol=" << btol << ", htol=" << htol << ", angtol=" << angtolDeg << ", Rsphere=" << Rsphere << "\n";

    IldaLog log;
    LammpsParser::Frame frame;
    if(!parseFrame(dump, frame)) return 1;
    const std::size_t n = static_cast<std::size_t>(frame.natoms);

    const OracleLabels labels = alignToOracle(frame, oracle);
    const std::vector<std::uint8_t> selection(n, 1);

    GrainRef gA = oracleGrain(oracle, "g1", 0, aA);
    gA.F = jsonToMatrix3RowMajor(oracle.at("FA"));
    GrainRef gB = oracleGrain(oracle, "g2", 1, aB);
    gB.F = jsonToMatrix3RowMajor(oracle.at("FB"));

    const auto& nArr = oracle.at("n");
    const Vector3 nNormal(nArr[0].get<double>(), nArr[1].get<double>(), nArr[2].get<double>());

    PTM ptm;
    IldaPtmNeighborQuery query(ptm, log);
    if(!matchPtm(ptm, query, frame, rmsd, log)) return 1;

    CisResult cis = findCIS(frame, labels.grain, labels.structureType, gA, gB, selection, query, cisTol, log);
    long sumCis = 0;
    for(const auto v : cis.cis) sumCis += v;

    IldaMesh mesh;
    const bool meshOk = buildCisSurfaceMesh(frame.simulationCell, frame, cis.cis, labels.grain,
                                            labels.structureType, gB.grainId, gB.structureType,
                                            Rsphere, mesh, log);

    Matrix3 VA, RA, VB, RB;
    polarDecomp(gA.F, VA, RA);
    polarDecomp(gB.F, VB, RB);
    const double bmagmax = std::max(gA.a * VA.trace() / 3.0, gB.a * VB.trace() / 3.0);

    CircuitContext ctx{
        frame, frame.simulationCell, frame.positions,
        labels.grain, labels.structureType, selection, cis, gA, gB,
        nNormal, query
    };

    std::vector<IldaSegment> segments;
    std::vector<DisconnectionMode> modes;
    double totalLength = 0.0;
    const double angtolRad = angtolDeg * PI / 180.0;
    findDefects(mesh, ctx, btol, htol, angtolRad, bmagmax, segments, modes, totalLength, log);

    IldaJsonExporter exporter;
    const std::string base = outPath.empty() ? std::string("/tmp/ilda_defects_selftest") : outPath;
    const std::string bondsPath = base + "_bonds.parquet";
    exporter.writeBondsParquet(segments, frame.simulationCell, bondsPath);
    log.dump(WrittenLineSegments{bondsPath});
    const std::string summaryPath = base + "_disconnection_summary.parquet";
    exporter.writeSummaryParquet(modes, totalLength, static_cast<int>(modes.size()), summaryPath);
    log.dump(WrittenSummary{modes, totalLength, summaryPath});

    int variantTotal = 0;
    for(const auto& m : modes) variantTotal += static_cast<int>(m.variants.size());

    json out;
    out["dump"] = dump;
    out["aligned"] = static_cast<int>(labels.aligned);
    out["sum_cis"] = static_cast<int>(sumCis);
    out["mesh_ok"] = meshOk;
    out["mesh_faces"] = static_cast<int>(mesh.faceCount());
    out["bmagmax"] = bmagmax;
    out["num_segments"] = static_cast<int>(segments.size());
    out["num_modes"] = static_cast<int>(modes.size());
    out["num_variants"] = variantTotal;
    out["total_line_length"] = totalLength;

    std::ostringstream summary;
    summary << "segments=" << segments.size() << " modes=" << modes.size() << " L=" << totalLength;
    writeResult(out, outPath, summary.str());
    std::cerr << "Defects self-test completed in " << timer.elapsed() << "\n";
    return 0;
}

}
