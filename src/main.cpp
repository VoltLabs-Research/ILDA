#include <volt/cli/common.h>
#include <volt/ilda_service.h>
#include <volt/ilda_options.h>
#include <volt/ilda_ptm_neighbor_query.h>
#include <volt/ilda_cis.h>
#include <volt/ilda_crf.h>
#include <volt/ilda_circuit.h>
#include <volt/ilda_mesh.h>
#include <volt/ilda_defects.h>
#include <volt/ilda_json_exporter.h>
#include <volt/ilda_types.h>
#include <volt/analysis/ptm.h>

#include <oneapi/tbb/global_control.h>
#include <tbb/info.h>

#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include <nlohmann/json.hpp>

using namespace Volt;
using namespace Volt::CLI;

namespace {

using json = nlohmann::json;

int runPtmSelftest(const std::map<std::string, std::string>& opts) {
    const std::string dump = getString(opts, "--selftest-ptm");
    const std::string oraclePath = getString(opts, "--oracle");
    const std::string outPath = getString(opts, "--out", "");
    if(dump.empty() || oraclePath.empty()) {
        std::cerr << "Usage: ilda --selftest-ptm <dump> --oracle <oracle.json> [--out <cpp.json>]\n";
        return 1;
    }

    std::ifstream oracleStream(oraclePath);
    if(!oracleStream.is_open()) {
        std::cerr << "Cannot open oracle: " << oraclePath << "\n";
        return 1;
    }
    json oracle = json::parse(oracleStream);
    const double rmsd = oracle.value("rmsd", 0.077);
    const auto qextArr = oracle.at("qext_xyzw");
    const Quaternion qext(qextArr[0].get<double>(), qextArr[1].get<double>(),
                          qextArr[2].get<double>(), qextArr[3].get<double>());

    LammpsParser::Frame frame;
    if(!parseFrame(dump, frame)) return 1;

    PTM ptm;
    ptm.setRmsdCutoff(rmsd);
    ptm.setCalculateDefGradient(false);
    if(!ptm.prepare(frame.positions.data(), frame.positions.size(), frame.simulationCell)) {
        std::cerr << "PTM prepare failed\n";
        return 1;
    }

    IldaPtmNeighborQuery query(ptm);
    query.computeAtomStates();

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

    if(outPath.empty()) {
        std::cout << out.dump(2) << "\n";
    } else {
        std::ofstream os(outPath);
        os << out.dump(2) << "\n";
        std::cerr << "Wrote " << outPath << " (" << atoms.size() << " atoms)\n";
    }
    return 0;
}

StructureType ovitoTypeToVolt(int t) {
    switch(t) {
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

int runCisSelftest(const std::map<std::string, std::string>& opts) {
    const std::string dump = getString(opts, "--selftest-cis");
    const std::string oraclePath = getString(opts, "--oracle");
    const std::string outPath = getString(opts, "--out", "");
    if(dump.empty() || oraclePath.empty()) {
        std::cerr << "Usage: ilda --selftest-cis <dump> --oracle <cis_oracle.json> [--out <cis_cpp.json>]\n";
        return 1;
    }

    std::ifstream oracleStream(oraclePath);
    if(!oracleStream.is_open()) {
        std::cerr << "Cannot open oracle: " << oraclePath << "\n";
        return 1;
    }
    json oracle = json::parse(oracleStream);
    const double rmsd = oracle.value("rmsd", 0.077);
    const double cisTol = oracle.value("cis_tol", 0.0);
    const double aA = oracle.value("aA", 1.0);
    const double aB = oracle.value("aB", 1.0);

    LammpsParser::Frame frame;
    if(!parseFrame(dump, frame)) return 1;
    const std::size_t n = static_cast<std::size_t>(frame.natoms);

    const auto& oraclePids = oracle.at("particle_id");
    const auto& oracleGrain = oracle.at("grain");
    const auto& oracleStype = oracle.at("structure_type");

    std::unordered_map<int, std::size_t> oraclePidToRow;
    oraclePidToRow.reserve(oraclePids.size());
    for(std::size_t r = 0; r < oraclePids.size(); ++r) {
        oraclePidToRow.emplace(oraclePids[r].get<int>(), r);
    }

    std::vector<int> grain(n, 0);
    std::vector<int> structureType(n, static_cast<int>(OTHER));
    std::vector<std::uint8_t> selection(n, 1);
    std::size_t aligned = 0;
    for(std::size_t i = 0; i < n; ++i) {
        const int pid = (i < frame.ids.size()) ? frame.ids[i] : static_cast<int>(i);
        auto it = oraclePidToRow.find(pid);
        if(it == oraclePidToRow.end()) continue;
        const std::size_t r = it->second;
        grain[i] = oracleGrain[r].get<int>();
        structureType[i] = static_cast<int>(ovitoTypeToVolt(oracleStype[r].get<int>()));
        ++aligned;
    }
    std::cerr << "Aligned " << aligned << "/" << n << " atoms to the oracle by particle id\n";

    GrainRef gA;
    gA.grainId = oracle.at("g1_id").get<int>();
    gA.grainIndex = 0;
    gA.a = aA;
    gA.structureType = ovitoTypeToVolt(oracle.at("g1_s").get<int>());
    {
        const auto& q = oracle.at("g1_orientation_xyzw");
        gA.tableOrientation = Quaternion(q[0].get<double>(), q[1].get<double>(),
                                         q[2].get<double>(), q[3].get<double>());
        gA.orientation = gA.tableOrientation;
    }

    GrainRef gB;
    gB.grainId = oracle.at("g2_id").get<int>();
    gB.grainIndex = 1;
    gB.a = aB;
    gB.structureType = ovitoTypeToVolt(oracle.at("g2_s").get<int>());
    {
        const auto& q = oracle.at("g2_orientation_xyzw");
        gB.tableOrientation = Quaternion(q[0].get<double>(), q[1].get<double>(),
                                         q[2].get<double>(), q[3].get<double>());
        gB.orientation = gB.tableOrientation;
    }

    PTM ptm;
    ptm.setRmsdCutoff(rmsd);
    ptm.setCalculateDefGradient(false);
    if(!ptm.prepare(frame.positions.data(), frame.positions.size(), frame.simulationCell)) {
        std::cerr << "PTM prepare failed\n";
        return 1;
    }
    IldaPtmNeighborQuery query(ptm);
    query.computeAtomStates();

    CisResult cis = findCIS(frame, grain, structureType, gA, gB, selection, query, cisTol);

    long cisCount = 0;
    std::vector<int> cisPids;
    for(std::size_t i = 0; i < n; ++i) {
        if(cis.cis[i]) {
            ++cisCount;
            cisPids.push_back((i < frame.ids.size()) ? frame.ids[i] : static_cast<int>(i));
        }
    }
    std::sort(cisPids.begin(), cisPids.end());

    json out;
    out["dump"] = dump;
    out["rmsd"] = rmsd;
    out["cis_tol"] = cisTol;
    out["natoms"] = static_cast<int>(n);
    out["sum_cis"] = static_cast<int>(cisCount);
    out["aligned"] = static_cast<int>(aligned);
    out["cis_particle_ids"] = cisPids;

    if(outPath.empty()) {
        std::cout << out.dump(2) << "\n";
    } else {
        std::ofstream os(outPath);
        os << out.dump(2) << "\n";
        std::cerr << "Wrote " << outPath << " (sum_cis=" << cisCount << ")\n";
    }
    return 0;
}

Matrix3 jsonToMatrix3RowMajor(const json& arr) {
    Matrix3 m{Matrix3::Zero()};
    for(int r = 0; r < 3; ++r)
        for(int c = 0; c < 3; ++c)
            m(static_cast<Matrix3::size_type>(r), static_cast<Matrix3::size_type>(c)) =
                arr[r * 3 + c].get<double>();
    return m;
}

json matrix3ToJsonRowMajor(const Matrix3& m) {
    json arr = json::array();
    for(int r = 0; r < 3; ++r)
        for(int c = 0; c < 3; ++c)
            arr.push_back(m(static_cast<Matrix3::size_type>(r), static_cast<Matrix3::size_type>(c)));
    return arr;
}

int runCrfComputeSelftest(const std::map<std::string, std::string>& opts) {
    const std::string oraclePath = getString(opts, "--oracle");
    const std::string outPath = getString(opts, "--out", "");
    if(oraclePath.empty()) {
        std::cerr << "Usage: ilda --selftest-crf-compute --oracle <oracle.json> [--out <cpp.json>]\n";
        return 1;
    }
    std::ifstream oracleStream(oraclePath);
    if(!oracleStream.is_open()) {
        std::cerr << "Cannot open oracle: " << oraclePath << "\n";
        return 1;
    }
    json oracle = json::parse(oracleStream);

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

    GrainRef gA, gB;
    gA.grainId = 1;
    gB.grainId = 2;
    computeF(xA, yA, xB, yB, EcohA, EcohB, coveraA, coveraB, gA, gB);

    json out;
    out["mode"] = "compute";
    out["FA"] = matrix3ToJsonRowMajor(gA.F);
    out["FB"] = matrix3ToJsonRowMajor(gB.F);
    out["g1_orientation_xyzw"] = {gA.orientation.x(), gA.orientation.y(), gA.orientation.z(), gA.orientation.w()};
    out["g2_orientation_xyzw"] = {gB.orientation.x(), gB.orientation.y(), gB.orientation.z(), gB.orientation.w()};

    if(outPath.empty()) {
        std::cout << out.dump(2) << "\n";
    } else {
        std::ofstream os(outPath);
        os << out.dump(2) << "\n";
        std::cerr << "Wrote " << outPath << " (crf-compute)\n";
    }
    return 0;
}

int runCrfEstimateSelftest(const std::map<std::string, std::string>& opts) {
    const std::string dump = getString(opts, "--selftest-crf-estimate");
    const std::string oraclePath = getString(opts, "--oracle");
    const std::string outPath = getString(opts, "--out", "");
    if(dump.empty() || oraclePath.empty()) {
        std::cerr << "Usage: ilda --selftest-crf-estimate <dump> --oracle <oracle.json> [--out <cpp.json>]\n";
        return 1;
    }
    std::ifstream oracleStream(oraclePath);
    if(!oracleStream.is_open()) {
        std::cerr << "Cannot open oracle: " << oraclePath << "\n";
        return 1;
    }
    json oracle = json::parse(oracleStream);
    const double rmsd = oracle.value("rmsd", 0.077);
    const double cisTol = oracle.value("cis_tol", 0.0);
    const double aA = oracle.value("aA", 1.0);
    const double aB = oracle.value("aB", 1.0);
    const double distF = oracle.value("distF", 10.0);

    LammpsParser::Frame frame;
    if(!parseFrame(dump, frame)) return 1;
    const std::size_t n = static_cast<std::size_t>(frame.natoms);

    const auto& oraclePids = oracle.at("particle_id");
    const auto& oracleGrain = oracle.at("grain");
    const auto& oracleStype = oracle.at("structure_type");
    std::unordered_map<int, std::size_t> oraclePidToRow;
    oraclePidToRow.reserve(oraclePids.size());
    for(std::size_t r = 0; r < oraclePids.size(); ++r) {
        oraclePidToRow.emplace(oraclePids[r].get<int>(), r);
    }

    std::vector<int> grain(n, 0);
    std::vector<int> structureType(n, static_cast<int>(OTHER));
    std::vector<std::uint8_t> selection(n, 1);
    std::size_t aligned = 0;
    for(std::size_t i = 0; i < n; ++i) {
        const int pid = (i < frame.ids.size()) ? frame.ids[i] : static_cast<int>(i);
        auto it = oraclePidToRow.find(pid);
        if(it == oraclePidToRow.end()) continue;
        const std::size_t r = it->second;
        grain[i] = oracleGrain[r].get<int>();
        structureType[i] = static_cast<int>(ovitoTypeToVolt(oracleStype[r].get<int>()));
        ++aligned;
    }
    std::cerr << "Aligned " << aligned << "/" << n << " atoms to the oracle by particle id\n";

    GrainRef gA;
    gA.grainId = oracle.at("g1_id").get<int>();
    gA.grainIndex = 0;
    gA.a = aA;
    gA.structureType = ovitoTypeToVolt(oracle.at("g1_s").get<int>());
    {
        const auto& q = oracle.at("g1_orientation_xyzw");
        gA.tableOrientation = Quaternion(q[0].get<double>(), q[1].get<double>(),
                                         q[2].get<double>(), q[3].get<double>());
        gA.orientation = gA.tableOrientation;
    }
    GrainRef gB;
    gB.grainId = oracle.at("g2_id").get<int>();
    gB.grainIndex = 1;
    gB.a = aB;
    gB.structureType = ovitoTypeToVolt(oracle.at("g2_s").get<int>());
    {
        const auto& q = oracle.at("g2_orientation_xyzw");
        gB.tableOrientation = Quaternion(q[0].get<double>(), q[1].get<double>(),
                                         q[2].get<double>(), q[3].get<double>());
        gB.orientation = gB.tableOrientation;
    }

    PTM ptm;
    ptm.setRmsdCutoff(rmsd);
    ptm.setCalculateDefGradient(false);
    if(!ptm.prepare(frame.positions.data(), frame.positions.size(), frame.simulationCell)) {
        std::cerr << "PTM prepare failed\n";
        return 1;
    }
    IldaPtmNeighborQuery query(ptm);
    query.computeAtomStates();

    CisResult cis = findCIS(frame, grain, structureType, gA, gB, selection, query, cisTol);
    estimateF(frame, cis, grain, structureType, selection, distF, query, gA, gB);

    json out;
    out["mode"] = "estimate";
    out["dump"] = dump;
    out["aligned"] = static_cast<int>(aligned);
    out["FA"] = matrix3ToJsonRowMajor(gA.F);
    out["FB"] = matrix3ToJsonRowMajor(gB.F);

    if(outPath.empty()) {
        std::cout << out.dump(2) << "\n";
    } else {
        std::ofstream os(outPath);
        os << out.dump(2) << "\n";
        std::cerr << "Wrote " << outPath << " (crf-estimate)\n";
    }
    return 0;
}


int runMeshSelftest(const std::map<std::string, std::string>& opts) {
    const std::string dump = getString(opts, "--selftest-mesh");
    const std::string oraclePath = getString(opts, "--oracle");
    const std::string outPath = getString(opts, "--out", "");
    if(dump.empty() || oraclePath.empty()) {
        std::cerr << "Usage: ilda --selftest-mesh <dump> --oracle <mesh_oracle.json> [--out <cpp.json>]\n";
        return 1;
    }

    std::ifstream oracleStream(oraclePath);
    if(!oracleStream.is_open()) {
        std::cerr << "Cannot open oracle: " << oraclePath << "\n";
        return 1;
    }
    json oracle = json::parse(oracleStream);
    const double Rsphere = oracle.value("Rsphere", 15.0);
    const int grainB = oracle.at("grain_b").get<int>();
    const int g2sVolt = static_cast<int>(ovitoTypeToVolt(oracle.at("g2_s").get<int>()));

    LammpsParser::Frame frame;
    if(!parseFrame(dump, frame)) return 1;
    const std::size_t n = static_cast<std::size_t>(frame.natoms);

    const auto& oraclePids = oracle.at("particle_id");
    const auto& oracleGrain = oracle.at("grain");
    const auto& oracleStype = oracle.at("structure_type");
    std::unordered_map<int, std::size_t> oraclePidToRow;
    oraclePidToRow.reserve(oraclePids.size());
    for(std::size_t r = 0; r < oraclePids.size(); ++r) {
        oraclePidToRow.emplace(oraclePids[r].get<int>(), r);
    }

    std::set<int> cisPidSet;
    for(const auto& pid : oracle.at("cis_particle_ids")) {
        cisPidSet.insert(pid.get<int>());
    }

    std::vector<int> grain(n, 0);
    std::vector<int> structureType(n, static_cast<int>(OTHER));
    std::vector<std::uint8_t> cis(n, 0);
    std::size_t aligned = 0;
    for(std::size_t i = 0; i < n; ++i) {
        const int pid = (i < frame.ids.size()) ? frame.ids[i] : static_cast<int>(i);
        auto it = oraclePidToRow.find(pid);
        if(it == oraclePidToRow.end()) continue;
        const std::size_t r = it->second;
        grain[i] = oracleGrain[r].get<int>();
        structureType[i] = static_cast<int>(ovitoTypeToVolt(oracleStype[r].get<int>()));
        cis[i] = cisPidSet.count(pid) ? 1 : 0;
        ++aligned;
    }
    std::cerr << "Aligned " << aligned << "/" << n << " atoms to the oracle by particle id\n";

    long sumCis = 0;
    for(const auto v : cis) sumCis += v;

    IldaMesh mesh;
    const bool ok = buildCisSurfaceMesh(frame.simulationCell, frame, cis, grain,
                                        structureType, grainB, g2sVolt, Rsphere, mesh);

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
        usedPids.push_back((static_cast<std::size_t>(atom) < frame.ids.size())
                               ? frame.ids[static_cast<std::size_t>(atom)] : atom);
    }
    std::sort(usedPids.begin(), usedPids.end());

    json out;
    out["dump"] = dump;
    out["Rsphere"] = Rsphere;
    out["ok"] = ok;
    out["aligned"] = static_cast<int>(aligned);
    out["sum_cis"] = static_cast<int>(sumCis);
    out["final_faces"] = static_cast<int>(faceCount);
    out["final_verts"] = static_cast<int>(mesh.vertexCount());
    out["final_used_verts"] = static_cast<int>(usedVerts.size());
    out["final_used_cis"] = usedCis;
    out["final_used_noncis"] = usedNonCis;
    out["final_used_pids"] = usedPids;

    if(outPath.empty()) {
        std::cout << out.dump(2) << "\n";
    } else {
        std::ofstream os(outPath);
        os << out.dump(2) << "\n";
        std::cerr << "Wrote " << outPath << " (final_faces=" << faceCount
                  << " used_cis=" << usedCis << " used_noncis=" << usedNonCis << ")\n";
    }
    return 0;
}


int runCircuitSelftest(const std::map<std::string, std::string>& opts) {
    const std::string dump = getString(opts, "--selftest-circuit");
    const std::string oraclePath = getString(opts, "--oracle");
    const std::string outPath = getString(opts, "--out", "");
    if(dump.empty() || oraclePath.empty()) {
        std::cerr << "Usage: ilda --selftest-circuit <dump> --oracle <oracle.json> [--out <cpp.json>]\n";
        return 1;
    }
    std::ifstream oracleStream(oraclePath);
    if(!oracleStream.is_open()) {
        std::cerr << "Cannot open oracle: " << oraclePath << "\n";
        return 1;
    }
    json oracle = json::parse(oracleStream);
    const double rmsd = oracle.value("rmsd", 0.077);
    const double cisTol = oracle.value("cis_tol", 0.0);
    const double aA = oracle.value("aA", 1.0);
    const double aB = oracle.value("aB", 1.0);

    LammpsParser::Frame frame;
    if(!parseFrame(dump, frame)) return 1;
    const std::size_t n = static_cast<std::size_t>(frame.natoms);

    const auto& oraclePids = oracle.at("particle_id");
    const auto& oracleGrain = oracle.at("grain");
    const auto& oracleStype = oracle.at("structure_type");
    std::unordered_map<int, std::size_t> oraclePidToRow;
    oraclePidToRow.reserve(oraclePids.size());
    for(std::size_t r = 0; r < oraclePids.size(); ++r) {
        oraclePidToRow.emplace(oraclePids[r].get<int>(), r);
    }

    std::vector<int> grain(n, 0);
    std::vector<int> structureType(n, static_cast<int>(OTHER));
    std::vector<std::uint8_t> selection(n, 1);
    std::size_t aligned = 0;
    for(std::size_t i = 0; i < n; ++i) {
        const int pid = (i < frame.ids.size()) ? frame.ids[i] : static_cast<int>(i);
        auto it = oraclePidToRow.find(pid);
        if(it == oraclePidToRow.end()) continue;
        const std::size_t r = it->second;
        grain[i] = oracleGrain[r].get<int>();
        structureType[i] = static_cast<int>(ovitoTypeToVolt(oracleStype[r].get<int>()));
        ++aligned;
    }
    std::cerr << "Aligned " << aligned << "/" << n << " atoms to the oracle by particle id\n";

    std::unordered_map<int, int> pidToIndex;
    pidToIndex.reserve(n);
    for(std::size_t i = 0; i < n; ++i) {
        const int pid = (i < frame.ids.size()) ? frame.ids[i] : static_cast<int>(i);
        pidToIndex.emplace(pid, static_cast<int>(i));
    }

    GrainRef gA;
    gA.grainId = oracle.at("g1_id").get<int>();
    gA.grainIndex = 0;
    gA.a = aA;
    gA.structureType = ovitoTypeToVolt(oracle.at("g1_s").get<int>());
    gA.F = jsonToMatrix3RowMajor(oracle.at("FA"));
    {
        const auto& q = oracle.at("g1_orientation_xyzw");
        gA.tableOrientation = Quaternion(q[0].get<double>(), q[1].get<double>(),
                                         q[2].get<double>(), q[3].get<double>());
        gA.orientation = gA.tableOrientation;
    }
    GrainRef gB;
    gB.grainId = oracle.at("g2_id").get<int>();
    gB.grainIndex = 1;
    gB.a = aB;
    gB.structureType = ovitoTypeToVolt(oracle.at("g2_s").get<int>());
    gB.F = jsonToMatrix3RowMajor(oracle.at("FB"));
    {
        const auto& q = oracle.at("g2_orientation_xyzw");
        gB.tableOrientation = Quaternion(q[0].get<double>(), q[1].get<double>(),
                                         q[2].get<double>(), q[3].get<double>());
        gB.orientation = gB.tableOrientation;
    }

    const auto& nArr = oracle.at("n");
    const Vector3 nNormal(nArr[0].get<double>(), nArr[1].get<double>(), nArr[2].get<double>());

    PTM ptm;
    ptm.setRmsdCutoff(rmsd);
    ptm.setCalculateDefGradient(false);
    if(!ptm.prepare(frame.positions.data(), frame.positions.size(), frame.simulationCell)) {
        std::cerr << "PTM prepare failed\n";
        return 1;
    }
    IldaPtmNeighborQuery query(ptm);
    query.computeAtomStates();

    CisResult cis = findCIS(frame, grain, structureType, gA, gB, selection, query, cisTol);

    CircuitContext ctx{
        frame, frame.simulationCell, frame.positions,
        grain, structureType, selection, cis, gA, gB,
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
        auto result = burgersCircuit(i1->second, i2->second, ctx);
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
    out["aligned"] = static_cast<int>(aligned);
    out["sum_cis"] = static_cast<int>(sumCis);
    out["circuits"] = circuitsOut;

    if(outPath.empty()) {
        std::cout << out.dump(2) << "\n";
    } else {
        std::ofstream os(outPath);
        os << out.dump(2) << "\n";
        std::cerr << "Wrote " << outPath << " (circuit; sum_cis=" << sumCis
                  << ", " << circuitsOut.size() << " circuits)\n";
    }
    return 0;
}


int runDefectsSelftest(const std::map<std::string, std::string>& opts) {
    const std::string dump = getString(opts, "--selftest-defects");
    const std::string oraclePath = getString(opts, "--oracle");
    const std::string outPath = getString(opts, "--out", "");
    if(dump.empty() || oraclePath.empty()) {
        std::cerr << "Usage: ilda --selftest-defects <dump> --oracle <circuit_oracle.json> [--out <cpp.json>]\n";
        return 1;
    }
    std::ifstream oracleStream(oraclePath);
    if(!oracleStream.is_open()) {
        std::cerr << "Cannot open oracle: " << oraclePath << "\n";
        return 1;
    }
    json oracle = json::parse(oracleStream);
    const double rmsd = oracle.value("rmsd", 0.077);
    const double cisTol = oracle.value("cis_tol", 0.0);
    const double aA = oracle.value("aA", 1.0);
    const double aB = oracle.value("aB", 1.0);
    const double btol = oracle.value("btol", 0.01);
    const double htol = oracle.value("htol", 0.5);
    const double angtolDeg = oracle.value("angtol", 5.0);
    const double Rsphere = oracle.value("Rsphere", 10.0);

    LammpsParser::Frame frame;
    if(!parseFrame(dump, frame)) return 1;
    const std::size_t n = static_cast<std::size_t>(frame.natoms);

    const auto& oraclePids = oracle.at("particle_id");
    const auto& oracleGrain = oracle.at("grain");
    const auto& oracleStype = oracle.at("structure_type");
    std::unordered_map<int, std::size_t> oraclePidToRow;
    oraclePidToRow.reserve(oraclePids.size());
    for(std::size_t r = 0; r < oraclePids.size(); ++r) {
        oraclePidToRow.emplace(oraclePids[r].get<int>(), r);
    }

    std::vector<int> grain(n, 0);
    std::vector<int> structureType(n, static_cast<int>(OTHER));
    std::vector<std::uint8_t> selection(n, 1);
    std::size_t aligned = 0;
    for(std::size_t i = 0; i < n; ++i) {
        const int pid = (i < frame.ids.size()) ? frame.ids[i] : static_cast<int>(i);
        auto it = oraclePidToRow.find(pid);
        if(it == oraclePidToRow.end()) continue;
        const std::size_t r = it->second;
        grain[i] = oracleGrain[r].get<int>();
        structureType[i] = static_cast<int>(ovitoTypeToVolt(oracleStype[r].get<int>()));
        ++aligned;
    }
    std::cerr << "Aligned " << aligned << "/" << n << " atoms to the oracle by particle id\n";

    GrainRef gA;
    gA.grainId = oracle.at("g1_id").get<int>();
    gA.grainIndex = 0;
    gA.a = aA;
    gA.structureType = ovitoTypeToVolt(oracle.at("g1_s").get<int>());
    gA.F = jsonToMatrix3RowMajor(oracle.at("FA"));
    {
        const auto& q = oracle.at("g1_orientation_xyzw");
        gA.tableOrientation = Quaternion(q[0].get<double>(), q[1].get<double>(),
                                         q[2].get<double>(), q[3].get<double>());
        gA.orientation = gA.tableOrientation;
    }
    GrainRef gB;
    gB.grainId = oracle.at("g2_id").get<int>();
    gB.grainIndex = 1;
    gB.a = aB;
    gB.structureType = ovitoTypeToVolt(oracle.at("g2_s").get<int>());
    gB.F = jsonToMatrix3RowMajor(oracle.at("FB"));
    {
        const auto& q = oracle.at("g2_orientation_xyzw");
        gB.tableOrientation = Quaternion(q[0].get<double>(), q[1].get<double>(),
                                         q[2].get<double>(), q[3].get<double>());
        gB.orientation = gB.tableOrientation;
    }

    const auto& nArr = oracle.at("n");
    const Vector3 nNormal(nArr[0].get<double>(), nArr[1].get<double>(), nArr[2].get<double>());

    PTM ptm;
    ptm.setRmsdCutoff(rmsd);
    ptm.setCalculateDefGradient(false);
    if(!ptm.prepare(frame.positions.data(), frame.positions.size(), frame.simulationCell)) {
        std::cerr << "PTM prepare failed\n";
        return 1;
    }
    IldaPtmNeighborQuery query(ptm);
    query.computeAtomStates();

    CisResult cis = findCIS(frame, grain, structureType, gA, gB, selection, query, cisTol);
    long sumCis = 0;
    for(const auto v : cis.cis) sumCis += v;

    IldaMesh mesh;
    const bool meshOk = buildCisSurfaceMesh(frame.simulationCell, frame, cis.cis, grain,
                                            structureType, gB.grainId, gB.structureType,
                                            Rsphere, mesh);

    Matrix3 VA, RA, VB, RB;
    polarDecomp(gA.F, VA, RA);
    polarDecomp(gB.F, VB, RB);
    const double trA = VA[0][0] + VA[1][1] + VA[2][2];
    const double trB = VB[0][0] + VB[1][1] + VB[2][2];
    const double bmagmax = std::max(gA.a * trA / 3.0, gB.a * trB / 3.0);

    CircuitContext ctx{
        frame, frame.simulationCell, frame.positions,
        grain, structureType, selection, cis, gA, gB,
        nNormal, query
    };

    std::vector<IldaSegment> segments;
    std::vector<DisconnectionMode> modes;
    double totalLength = 0.0;
    const double angtolRad = angtolDeg * PI / 180.0;
    findDefects(mesh, ctx, btol, htol, angtolRad, bmagmax, segments, modes, totalLength);

    IldaJsonExporter exporter;
    const std::string base = outPath.empty() ? std::string("/tmp/ilda_defects_selftest") : outPath;
    exporter.writeBondsParquet(segments, frame.simulationCell, base + "_bonds.parquet");
    exporter.writeSummaryParquet(modes, totalLength, static_cast<int>(modes.size()),
                                 base + "_disconnection_summary.parquet");

    int variantTotal = 0;
    for(const auto& m : modes) variantTotal += static_cast<int>(m.variants.size());

    json out;
    out["dump"] = dump;
    out["aligned"] = static_cast<int>(aligned);
    out["sum_cis"] = static_cast<int>(sumCis);
    out["mesh_ok"] = meshOk;
    out["mesh_faces"] = static_cast<int>(mesh.faceCount());
    out["bmagmax"] = bmagmax;
    out["num_segments"] = static_cast<int>(segments.size());
    out["num_modes"] = static_cast<int>(modes.size());
    out["num_variants"] = variantTotal;
    out["total_line_length"] = totalLength;

    if(outPath.empty()) {
        std::cout << out.dump(2) << "\n";
    } else {
        std::ofstream os(outPath);
        os << out.dump(2) << "\n";
        std::cerr << "Wrote " << outPath << " (segments=" << segments.size()
                  << " modes=" << modes.size() << " L=" << totalLength << ")\n";
    }
    return 0;
}


Vector3 parseVector3Option(const std::map<std::string, std::string>& opts,
                           const std::string& key, const Vector3& fallback) {
    const std::string raw = getString(opts, key, "");
    if(raw.empty()) return fallback;
    try {
        json parsed = json::parse(raw);
        if(parsed.is_array() && parsed.size() == 3) {
            return Vector3(parsed[0].get<double>(), parsed[1].get<double>(), parsed[2].get<double>());
        }
        if(parsed.is_object()) {
            return Vector3(
                parsed.value("x", fallback.x()),
                parsed.value("y", fallback.y()),
                parsed.value("z", fallback.z())
            );
        }
    } catch(const std::exception& e) {
        spdlog::warn("Could not parse vector option {} = '{}': {}", key, raw, e.what());
    }
    return fallback;
}

Matrix3 parseMatrix3Option(const std::map<std::string, std::string>& opts,
                           const std::string& key, const Matrix3& fallback) {
    const std::string raw = getString(opts, key, "");
    if(raw.empty()) return fallback;
    try {
        json parsed = json::parse(raw);
        double m[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        if(parsed.is_array() && parsed.size() == 9) {
            for(int i = 0; i < 9; ++i) m[i] = parsed[i].get<double>();
        } else if(parsed.is_array() && parsed.size() == 3) {
            for(int r = 0; r < 3; ++r) {
                const auto& row = parsed[r];
                if(row.is_array() && row.size() == 3) {
                    for(int c = 0; c < 3; ++c) m[r * 3 + c] = row[c].get<double>();
                } else if(row.is_object()) {
                    m[r * 3 + 0] = row.value("x", 0.0);
                    m[r * 3 + 1] = row.value("y", 0.0);
                    m[r * 3 + 2] = row.value("z", 0.0);
                } else {
                    return fallback;
                }
            }
        } else {
            return fallback;
        }
        Matrix3 result{Matrix3::Zero()};
        for(int r = 0; r < 3; ++r)
            for(int c = 0; c < 3; ++c)
                result(static_cast<Matrix3::size_type>(r), static_cast<Matrix3::size_type>(c)) = m[r * 3 + c];
        return result;
    } catch(const std::exception& e) {
        spdlog::warn("Could not parse matrix option {} = '{}': {}", key, raw, e.what());
    }
    return fallback;
}

void showUsage(const std::string& name) {
    printUsageHeader(name, "Volt - Interfacial Line Defect Analysis (ILDA)");
    std::cerr
        << "  --grain-atoms <path>      Path to upstream *_atoms.parquet (grain-segmentation).\n"
        << "  --grains <path>           Path to upstream *_grains.parquet (grain-segmentation).\n"
        << "  --atomA <id>              Particle Identifier of fiducial atom A.\n"
        << "  --atomB <id>              Particle Identifier of fiducial atom B.\n"
        << "  --circuitAtom1 <id>       Particle Identifier for single-circuit start atom.\n"
        << "  --circuitAtom2 <id>       Particle Identifier for single-circuit end atom.\n"
        << "  --aA/--cA/--aB/--cB <f>   Lattice constants (a/c axes) for grains A and B.\n"
        << "  --typeA/--typeB <int>     Structure type override (-1 => derive). [default: -1]\n"
        << "  --Rsphere <float>         Probe sphere radius. [default: 10]\n"
        << "  --htol <float>            Step height tolerance. [default: 0.5]\n"
        << "  --btol <float>            Burgers vector length tolerance. [default: 0.01]\n"
        << "  --angtol <float>          Burgers vector angular tolerance (deg). [default: 5]\n"
        << "  --distF <float>           Interface skin distance. [default: 10]\n"
        << "  --cis_tol <float>         Co-incidence site tolerance. [default: 0]\n"
        << "  --n '{\"x\":..}'            Interface plane normal vector (JSON).\n"
        << "  --xA/--yA/--xB/--yB '{}'  Orientation vectors (JSON; used when estimateF=false).\n"
        << "  --EcohA/--EcohB '[[..]]'  Coherency strain matrices (JSON; used when estimateF=false).\n"
        << "  --estimateF [bool]        Estimate coherent reference frame. [default: true]\n"
        << "  --single_circuit          Run a single Burgers circuit only.\n"
        << "  --extract_lines           Run the full line-extraction pipeline.\n"
        << "  --selection_only          Use only selected particles.\n"
        << "  --print_results           Print a results summary.\n"
        << "  --crystalStructure <t>    Reference crystal structure for PTM. [default: FCC]\n"
        << "  --rmsd <float>            PTM RMSD cutoff. [default: 0.1]\n"
        << "  --threads <int>           Max worker threads (TBB). [default: physical cores]\n";
    printHelpOption();
}

}

int main(int argc, char* argv[]) {
    if(argc < 2) {
        showUsage(argv[0]);
        return 1;
    }

    std::string filename;
    std::string outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);

    if(hasOption(opts, "--selftest-ptm")) {
        initLogging("ilda-selftest");
        return runPtmSelftest(opts);
    }

    if(hasOption(opts, "--selftest-cis")) {
        initLogging("ilda-selftest");
        return runCisSelftest(opts);
    }

    if(hasOption(opts, "--selftest-crf-compute")) {
        initLogging("ilda-selftest");
        return runCrfComputeSelftest(opts);
    }

    if(hasOption(opts, "--selftest-crf-estimate")) {
        initLogging("ilda-selftest");
        return runCrfEstimateSelftest(opts);
    }

    if(hasOption(opts, "--selftest-mesh")) {
        initLogging("ilda-selftest");
        return runMeshSelftest(opts);
    }

    if(hasOption(opts, "--selftest-circuit")) {
        initLogging("ilda-selftest");
        return runCircuitSelftest(opts);
    }

    if(hasOption(opts, "--selftest-defects")) {
        initLogging("ilda-selftest");
        return runDefectsSelftest(opts);
    }

    if(hasOption(opts, "--help") || filename.empty()) {
        showUsage(argv[0]);
        return filename.empty() ? 1 : 0;
    }

    if(!hasOption(opts, "--threads")) {
        const int maxAvailableThreads = static_cast<int>(oneapi::tbb::info::default_concurrency());
        int physicalCores = 0;
        std::ifstream cpuinfo("/proc/cpuinfo");
        if(cpuinfo.is_open()) {
            std::set<std::pair<int, int>> physicalCoreIds;
            int fallbackCpuCores = 0;
            int physicalId = -1;
            int coreId = -1;
            std::string line;
            while(std::getline(cpuinfo, line)) {
                if(line.empty()) {
                    if(physicalId >= 0 && coreId >= 0) {
                        physicalCoreIds.emplace(physicalId, coreId);
                    }
                    physicalId = -1;
                    coreId = -1;
                    continue;
                }
                int parsedValue = 0;
                const auto separator = line.find(':');
                if(separator != std::string::npos) {
                    try {
                        parsedValue = std::stoi(line.substr(separator + 1));
                    } catch(const std::exception&) {
                        parsedValue = 0;
                    }
                }
                if(line.rfind("physical id", 0) == 0) {
                    physicalId = parsedValue;
                } else if(line.rfind("core id", 0) == 0) {
                    coreId = parsedValue;
                } else if(line.rfind("cpu cores", 0) == 0) {
                    fallbackCpuCores = std::max(fallbackCpuCores, parsedValue);
                }
            }
            if(physicalId >= 0 && coreId >= 0) {
                physicalCoreIds.emplace(physicalId, coreId);
            }
            physicalCores = physicalCoreIds.empty() ? fallbackCpuCores : static_cast<int>(physicalCoreIds.size());
        }
        opts["--threads"] = std::to_string(std::max(1, std::min(maxAvailableThreads, physicalCores > 0 ? physicalCores : maxAvailableThreads)));
    }

    const int requestedThreads = getInt(opts, "--threads");
    oneapi::tbb::global_control parallelControl(
        oneapi::tbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(std::max(1, requestedThreads))
    );
    initLogging("ilda");
    spdlog::info("Using {} threads (OneTBB)", requestedThreads);

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)) return 1;

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    IldaOptions options;
    options.grainAtomsPath = getString(opts, "--grain-atoms");
    options.grainsPath = getString(opts, "--grains");

    options.atomA = getInt(opts, "--atomA", options.atomA);
    options.atomB = getInt(opts, "--atomB", options.atomB);
    options.circuitAtom1 = getInt(opts, "--circuitAtom1", options.circuitAtom1);
    options.circuitAtom2 = getInt(opts, "--circuitAtom2", options.circuitAtom2);

    options.aA = getDouble(opts, "--aA", options.aA);
    options.cA = getDouble(opts, "--cA", options.cA);
    options.aB = getDouble(opts, "--aB", options.aB);
    options.cB = getDouble(opts, "--cB", options.cB);

    options.typeA = getInt(opts, "--typeA", options.typeA);
    options.typeB = getInt(opts, "--typeB", options.typeB);

    options.cis_tol = getDouble(opts, "--cis_tol", options.cis_tol);
    options.Rsphere = getDouble(opts, "--Rsphere", options.Rsphere);
    options.htol = getDouble(opts, "--htol", options.htol);
    options.btol = getDouble(opts, "--btol", options.btol);
    options.angtol = getDouble(opts, "--angtol", options.angtol);
    options.distF = getDouble(opts, "--distF", options.distF);

    options.estimateF = getBool(opts, "--estimateF", options.estimateF);
    options.single_circuit = getBool(opts, "--single_circuit", options.single_circuit);
    options.extract_lines = getBool(opts, "--extract_lines", options.extract_lines);
    options.selection_only = getBool(opts, "--selection_only", options.selection_only);
    options.print_results = getBool(opts, "--print_results", options.print_results);

    options.n = parseVector3Option(opts, "--n", options.n);
    options.xA = parseVector3Option(opts, "--xA", options.xA);
    options.yA = parseVector3Option(opts, "--yA", options.yA);
    options.xB = parseVector3Option(opts, "--xB", options.xB);
    options.yB = parseVector3Option(opts, "--yB", options.yB);
    options.EcohA = parseMatrix3Option(opts, "--EcohA", options.EcohA);
    options.EcohB = parseMatrix3Option(opts, "--EcohB", options.EcohB);

    options.crystalStructure = getString(opts, "--crystalStructure", options.crystalStructure);
    options.rmsdCutoff = getDouble(opts, "--rmsd", options.rmsdCutoff);

    spdlog::info("Starting ILDA analysis...");
    ILDAService analyzer;
    analyzer.setOptions(options);
    json result = analyzer.run(frame, outputBase);
    if(result.value("is_failed", false)) {
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }
    {
        const std::string resultPath = outputBase + "_result.json";
        std::ofstream resultStream(resultPath);
        if(resultStream.is_open()) {
            resultStream << result.dump(2) << "\n";
            spdlog::info("Wrote {}", resultPath);
        } else {
            spdlog::warn("Could not write result document to {}", resultPath);
        }
    }
    spdlog::info("Analysis completed successfully.");
    return 0;
}
