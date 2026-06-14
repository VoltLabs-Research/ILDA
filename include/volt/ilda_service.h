#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/analysis_result.h>
#include <volt/ilda_options.h>
#include <volt/ilda_types.h>
#include <volt/ilda_json_exporter.h>
#include <nlohmann/json.hpp>

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace Volt {

using json = nlohmann::json;

class ILDA {
public:
    ILDA();

    void setOptions(const IldaOptions& options) { _options = options; }
    const IldaOptions& options() const { return _options; }

    json run(const LammpsParser::Frame& frame, const std::string& outputBase);

private:
    bool loadUpstream(const LammpsParser::Frame& frame, std::string& error);

    int particleIdToIndex(int particleId) const;

    bool resolveGrains(std::string& error);

    IldaOptions _options;
    IldaJsonExporter _jsonExporter;

    std::vector<int> _grain;
    std::vector<int> _structureType;
    std::vector<int> _particleIds;
    std::vector<std::uint8_t> _selection;
    std::unordered_map<int, int> _particleIdToIndex;

    struct GrainTableEntry {
        int id = 0;
        int structureType = 0;
        Quaternion orientation{0.0, 0.0, 0.0, 1.0};
    };
    std::unordered_map<int, GrainTableEntry> _grainTable;

    GrainRef _gA;
    GrainRef _gB;
};

using ILDAService = ILDA;

}
