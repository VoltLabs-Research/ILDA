#include <volt/ilda_service.h>
#include <volt/ilda_cis.h>
#include <volt/ilda_crf.h>
#include <volt/ilda_mesh.h>
#include <volt/ilda_circuit.h>
#include <volt/ilda_defects.h>
#include <volt/ilda_ptm_neighbor_query.h>
#include <volt/analysis/ptm.h>

#include <duckdb.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <spdlog/spdlog.h>

namespace Volt {

namespace {

json readPayloadParquet(const std::string& path) {
    duckdb::DuckDB db(nullptr);
    duckdb::Connection con(db);

    std::string quoted = "'";
    for(char c : path) {
        if(c == '\'') quoted.push_back('\'');
        quoted.push_back(c);
    }
    quoted.push_back('\'');

    auto result = con.Query("SELECT payload FROM read_parquet(" + quoted + ")");
    if(result->HasError()) {
        throw std::runtime_error("DuckDB read of " + path + " failed: " + result->GetError());
    }
    if(result->RowCount() == 0) {
        throw std::runtime_error("Parquet " + path + " has no rows");
    }
    const std::string payload = result->GetValue(0, 0).ToString();
    return json::parse(payload);
}

}

ILDA::ILDA() = default;

int ILDA::particleIdToIndex(int particleId) const {
    auto it = _particleIdToIndex.find(particleId);
    if(it == _particleIdToIndex.end()) return -1;
    return it->second;
}

bool ILDA::loadUpstream(const LammpsParser::Frame& frame, std::string& error) {
    const std::size_t n = static_cast<std::size_t>(frame.natoms);

    _particleIds = frame.ids;
    if(_particleIds.size() != n) {
        _particleIds.resize(n);
        for(std::size_t i = 0; i < n; ++i) _particleIds[i] = static_cast<int>(i);
    }
    _particleIdToIndex.clear();
    _particleIdToIndex.reserve(n);
    for(std::size_t i = 0; i < n; ++i) {
        _particleIdToIndex.emplace(_particleIds[i], static_cast<int>(i));
    }

    if(_options.grainAtomsPath.empty() || _options.grainsPath.empty()) {
        error = "ILDA requires grain-segmentation upstream artifacts "
                "(--grain-atoms <_atoms.parquet> and --grains <_grains.parquet>)";
        return false;
    }

    json atomsDoc;
    try {
        atomsDoc = readPayloadParquet(_options.grainAtomsPath);
    } catch(const std::exception& e) {
        error = std::string("Could not read grain-atoms parquet: ") + e.what();
        return false;
    }
    if(!atomsDoc.contains("per-atom-properties") || !atomsDoc["per-atom-properties"].is_array()) {
        error = "grain-atoms parquet missing 'per-atom-properties' array";
        return false;
    }

    _grain.assign(n, 0);
    _structureType.assign(n, static_cast<int>(OTHER));
    _selection.assign(n, 1);

    std::size_t aligned = 0;
    for(const auto& a : atomsDoc["per-atom-properties"]) {
        const int pid = a.value("id", -1);
        auto it = _particleIdToIndex.find(pid);
        if(it == _particleIdToIndex.end()) continue;
        const std::size_t idx = static_cast<std::size_t>(it->second);
        _grain[idx] = a.value("grain_id", 0);
        _structureType[idx] = a.value("structure_id", static_cast<int>(OTHER));
        ++aligned;
    }
    spdlog::info("ILDA: aligned {}/{} atoms from {}", aligned, n, _options.grainAtomsPath);
    if(aligned != n) {
        spdlog::warn("ILDA: {} atoms in the frame were not present in the grain-atoms parquet",
                     n - aligned);
    }

    json grainsDoc;
    try {
        grainsDoc = readPayloadParquet(_options.grainsPath);
    } catch(const std::exception& e) {
        error = std::string("Could not read grains parquet: ") + e.what();
        return false;
    }
    if(!grainsDoc.contains("sub_listings") ||
       !grainsDoc["sub_listings"].contains("grains") ||
       !grainsDoc["sub_listings"]["grains"].is_array()) {
        error = "grains parquet missing 'sub_listings.grains' array";
        return false;
    }
    _grainTable.clear();
    for(const auto& g : grainsDoc["sub_listings"]["grains"]) {
        GrainTableEntry e;
        e.id = g.value("id", 0);
        e.structureType = g.value("structure_type", static_cast<int>(OTHER));
        if(g.contains("orientation") && g["orientation"].is_array() && g["orientation"].size() == 4) {
            const auto& q = g["orientation"];
            e.orientation = Quaternion(q[0].get<double>(), q[1].get<double>(),
                                       q[2].get<double>(), q[3].get<double>());
        }
        _grainTable[e.id] = e;
    }
    spdlog::info("ILDA: loaded {} grains from {}", _grainTable.size(), _options.grainsPath);

    return true;
}

bool ILDA::resolveGrains(std::string& error) {
    const int idxA = particleIdToIndex(_options.atomA);
    if(idxA < 0) {
        error = "Invalid Particle Identifier for atomA";
        return false;
    }
    const int idxB = particleIdToIndex(_options.atomB);
    if(idxB < 0) {
        error = "Invalid Particle Identifier for atomB";
        return false;
    }

    const int g1Id = _grain[static_cast<std::size_t>(idxA)];
    const int g2Id = _grain[static_cast<std::size_t>(idxB)];
    if(g1Id == g2Id) {
        error = "atomA and atomB are both part of the same grain";
        return false;
    }

    auto resolveOne = [&](GrainRef& g, int grainId, int slot, int typeOverride,
                          double a, double c) -> bool {
        g.grainId = grainId;
        g.grainIndex = slot;
        g.a = a;
        g.covera = (a > 0.0) ? c / a : 0.0;

        auto it = _grainTable.find(grainId);
        if(it == _grainTable.end()) {
            error = "Grain " + std::to_string(grainId) +
                    " (from fiducial atom) is not present in the grains table";
            return false;
        }
        g.structureType = (typeOverride > 0)
            ? static_cast<StructureType>(typeOverride)
            : static_cast<StructureType>(it->second.structureType);
        g.tableOrientation = it->second.orientation;
        g.orientation = it->second.orientation;
        return true;
    };

    if(!resolveOne(_gA, g1Id, 0, _options.typeA, _options.aA, _options.cA)) return false;
    if(!resolveOne(_gB, g2Id, 1, _options.typeB, _options.aB, _options.cB)) return false;

    spdlog::info("ILDA: grain A id={} type={} | grain B id={} type={}",
                 _gA.grainId, structureTypeName(_gA.structureType),
                 _gB.grainId, structureTypeName(_gB.structureType));
    return true;
}

json ILDA::run(const LammpsParser::Frame& frame, const std::string& outputBase) {
    spdlog::info("ILDA: starting analysis");

    std::string error;
    if(!loadUpstream(frame, error)) {
        return AnalysisResult::failure(error.empty() ? "Failed to load upstream artifacts" : error);
    }
    if(!resolveGrains(error)) {
        return AnalysisResult::failure(error.empty() ? "Failed to resolve grains" : error);
    }

    if(!_options.estimateF) {
        const double coveraA = (_options.aA > 0.0) ? _options.cA / _options.aA : 0.0;
        const double coveraB = (_options.aB > 0.0) ? _options.cB / _options.aB : 0.0;
        computeF(_options.xA, _options.yA, _options.xB, _options.yB,
                 _options.EcohA, _options.EcohB, coveraA, coveraB, _gA, _gB);
    }

    PTM ptm;
    ptm.setRmsdCutoff(_options.rmsdCutoff);
    ptm.setCalculateDefGradient(false);
    if(!ptm.prepare(frame.positions.data(), frame.positions.size(), frame.simulationCell)) {
        return AnalysisResult::failure("PTM prepare failed");
    }
    IldaPtmNeighborQuery query(ptm);
    query.computeAtomStates();

    CisResult cis = findCIS(frame, _grain, _structureType, _gA, _gB, _selection, query, _options.cis_tol);

    const long cisCount = std::count(cis.cis.begin(), cis.cis.end(), static_cast<std::uint8_t>(1));
    spdlog::info("Found {} CIS atoms", cisCount);

    json result = AnalysisResult::success();
    result["cis_count"] = static_cast<int>(cisCount);
    result["grain_a"] = _gA.grainId;
    result["grain_b"] = _gB.grainId;

    if(_options.single_circuit || _options.extract_lines) {
        if(cisCount < 3) {
            return AnalysisResult::failure("Not enough CIS atoms to proceed with line extraction");
        }

        if(_options.estimateF) {
            estimateF(frame, cis, _grain, _structureType, _selection, _options.distF, query, _gA, _gB);
        } else {
            spdlog::info("User specified coherent reference frame");
        }

        Matrix3 VA, RA, VB, RB;
        polarDecomp(_gA.F, VA, RA);
        polarDecomp(_gB.F, VB, RB);
        const double trA = VA[0][0] + VA[1][1] + VA[2][2];
        const double trB = VB[0][0] + VB[1][1] + VB[2][2];
        const double bmagmax = std::max(_gA.a * trA / 3.0, _gB.a * trB / 3.0);

        CircuitContext ctx{
            frame, frame.simulationCell, frame.positions,
            _grain, _structureType, _selection, cis, _gA, _gB,
            _options.n, query
        };

        if(_options.single_circuit) {
            const int idx1 = particleIdToIndex(_options.circuitAtom1);
            const int idx2 = particleIdToIndex(_options.circuitAtom2);
            if(idx1 < 0 || idx2 < 0) {
                spdlog::error("Error: invalid single circuit atom identifier");
            } else if(!cis.cis[static_cast<std::size_t>(idx1)] || !cis.cis[static_cast<std::size_t>(idx2)]) {
                spdlog::error("Error: Cannot construct specified Burgers circuit because "
                              "one or both atoms is not a CIS atom");
            } else {
                auto circuitResult = burgersCircuit(idx1, idx2, ctx);
                if(circuitResult) {
                    spdlog::info("Info for Burgers circuit using atoms {} and {}:",
                                 _options.circuitAtom1, _options.circuitAtom2);
                    spdlog::info("b = ({}, {}, {})",
                                 circuitResult->latticeSum.x(), circuitResult->latticeSum.y(),
                                 circuitResult->latticeSum.z());
                    spdlog::info("Ideal step height A = {:.3f}", circuitResult->hA);
                    spdlog::info("Ideal step height B = {:.3f}", circuitResult->hB);
                    result["single_circuit"] = {
                        {"b", {circuitResult->latticeSum.x(), circuitResult->latticeSum.y(),
                               circuitResult->latticeSum.z()}},
                        {"hA", circuitResult->hA},
                        {"hB", circuitResult->hB}
                    };
                } else {
                    spdlog::warn("single_circuit produced no closed circuit");
                }
            }
        }

        if(_options.extract_lines) {
            IldaMesh mesh;
            if(!buildCisSurfaceMesh(frame.simulationCell, frame, cis.cis, _grain,
                                    _structureType, _gB.grainId, _gB.structureType,
                                    _options.Rsphere, mesh)) {
                return AnalysisResult::failure("Failed to build the CIS surface mesh");
            }

            std::vector<IldaSegment> segments;
            std::vector<DisconnectionMode> modes;
            double totalLength = 0.0;
            const double angtolRad = _options.angtol * PI / 180.0;
            findDefects(mesh, ctx, _options.btol, _options.htol, angtolRad,
                        bmagmax, segments, modes, totalLength);

            _jsonExporter.writeBondsParquet(segments, frame.simulationCell, outputBase + "_bonds.parquet");
            _jsonExporter.writeSummaryParquet(modes, totalLength,
                                              static_cast<int>(modes.size()),
                                              outputBase + "_disconnection_summary.parquet");

            spdlog::info("ILDA: emitted {} line segments across {} disconnection modes "
                         "(total line length {:.3f})",
                         segments.size(), modes.size(), totalLength);
            result["num_segments"] = static_cast<int>(segments.size());
            result["num_modes"] = static_cast<int>(modes.size());
            result["total_line_length"] = totalLength;
        }
    }

    spdlog::info("ILDA: analysis completed");
    return result;
}

}
