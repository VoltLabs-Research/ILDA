#include <volt/ilda_plugin.h>
#include <volt/ilda_selftest.h>
#include <volt/ilda_service.h>
#include <volt/plugin/option_binding.h>

#include <fstream>
#include <string>
#include <vector>

using namespace Volt;
using namespace Volt::Plugin;
using namespace Volt::IldaPlugin;

static const std::vector<OptionBinding<IldaOptions>> bindings = {
    pathOption("--grain-atoms", "Per-atom table from grain-segmentation.", &IldaOptions::grainAtomsPath),
    pathOption("--grains", "Grains table from grain-segmentation.", &IldaOptions::grainsPath),
    opt("--atomA", "Particle identifier of fiducial atom A.", &IldaOptions::atomA),
    opt("--atomB", "Particle identifier of fiducial atom B.", &IldaOptions::atomB),
    opt("--circuitAtom1", "Particle identifier for single-circuit start atom.", &IldaOptions::circuitAtom1),
    opt("--circuitAtom2", "Particle identifier for single-circuit end atom.", &IldaOptions::circuitAtom2),
    opt("--aA", "Lattice constant a for grain A.", &IldaOptions::aA),
    opt("--cA", "Lattice constant c for grain A.", &IldaOptions::cA),
    opt("--aB", "Lattice constant a for grain B.", &IldaOptions::aB),
    opt("--cB", "Lattice constant c for grain B.", &IldaOptions::cB),
    opt("--typeA", "Structure type override for grain A (-1 derives it).", &IldaOptions::typeA),
    opt("--typeB", "Structure type override for grain B (-1 derives it).", &IldaOptions::typeB),
    opt("--Rsphere", "Probe sphere radius.", &IldaOptions::Rsphere),
    opt("--htol", "Step height tolerance.", &IldaOptions::htol),
    opt("--btol", "Burgers vector length tolerance.", &IldaOptions::btol),
    opt("--angtol", "Burgers vector angular tolerance (degrees).", &IldaOptions::angtol),
    opt("--distF", "Interface skin distance.", &IldaOptions::distF),
    opt("--cis_tol", "Co-incidence site tolerance.", &IldaOptions::cis_tol),
    vectorOption("--n", "Interface plane normal vector, as JSON {\"x\":..,\"y\":..,\"z\":..}.", &IldaOptions::n),
    vectorOption("--xA", "Grain A orientation vector x, as JSON (used when estimateF is false).", &IldaOptions::xA),
    vectorOption("--yA", "Grain A orientation vector y, as JSON (used when estimateF is false).", &IldaOptions::yA),
    vectorOption("--xB", "Grain B orientation vector x, as JSON (used when estimateF is false).", &IldaOptions::xB),
    vectorOption("--yB", "Grain B orientation vector y, as JSON (used when estimateF is false).", &IldaOptions::yB),
    matrixOption("--EcohA", "Grain A coherency strain matrix, as JSON (used when estimateF is false).", &IldaOptions::EcohA),
    matrixOption("--EcohB", "Grain B coherency strain matrix, as JSON (used when estimateF is false).", &IldaOptions::EcohB),
    opt("--estimateF", "Estimate the coherent reference frame.", &IldaOptions::estimateF),
    opt("--single_circuit", "Run a single Burgers circuit only.", &IldaOptions::single_circuit),
    opt("--extract_lines", "Run the full line-extraction pipeline.", &IldaOptions::extract_lines),
    opt("--selection_only", "Use only selected particles.", &IldaOptions::selection_only),
    opt("--print_results", "Print a results summary.", &IldaOptions::print_results),
    opt<IldaOptions>(CliOption{"--crystalStructure", "enum", "Reference crystal structure for PTM.",
            IldaOptions{}.crystalStructure, {"SC", "FCC", "HCP", "BCC", "CUBIC_DIAMOND", "HEX_DIAMOND"}},
        [](IldaOptions& options, const OptsMap& opts) {
            options.crystalStructure = CLI::getString(opts, "--crystalStructure", options.crystalStructure);
        }),
    opt("--rmsd", "PTM RMSD cutoff.", &IldaOptions::rmsdCutoff),
};

namespace Volt::IldaPlugin {

Binding pathOption(const char* name, const char* help, std::string IldaOptions::* field) {
    return {{name, "path", help, ""},
        [name, field](IldaOptions& options, const OptsMap& opts) {
            options.*field = CLI::getString(opts, name);
        }};
}

Binding vectorOption(const char* name, const char* help, Vector3 IldaOptions::* field) {
    return {{name, "string", help, ""},
        [name, field](IldaOptions& options, const OptsMap& opts) {
            options.*field = parseVector3Option(opts, name, options.*field);
        }};
}

Binding matrixOption(const char* name, const char* help, Matrix3 IldaOptions::* field) {
    return {{name, "string", help, ""},
        [name, field](IldaOptions& options, const OptsMap& opts) {
            options.*field = parseMatrix3Option(opts, name, options.*field);
        }};
}

Vector3 parseVector3Option(const OptsMap& opts, const std::string& key, const Vector3& fallback) {
    const std::string raw = CLI::getString(opts, key, "");
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

Matrix3 parseMatrix3Option(const OptsMap& opts, const std::string& key, const Matrix3& fallback) {
    const std::string raw = CLI::getString(opts, key, "");
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

void logResolvedParameters(const IldaOptions& options) {
    spdlog::info("Parameters: atomA={}, atomB={}, aA={}, cA={}, aB={}, cB={}, typeA={}, typeB={}, rmsd={}, "
                 "cis_tol={}, estimateF={}, distF={}, extract_lines={}, Rsphere={}, btol={}, htol={}, "
                 "angtol={}, n=({}, {}, {}), single_circuit={}, circuitAtom1={}, circuitAtom2={}",
                 options.atomA, options.atomB, options.aA, options.cA, options.aB, options.cB,
                 options.typeA, options.typeB, options.rmsdCutoff, options.cis_tol, options.estimateF,
                 options.distF, options.extract_lines, options.Rsphere, options.btol, options.htol,
                 options.angtol, options.n.x(), options.n.y(), options.n.z(), options.single_circuit,
                 options.circuitAtom1, options.circuitAtom2);
}

void writeResultDocument(const json& result, const std::string& path) {
    std::ofstream resultStream(path);
    if(resultStream.is_open()) {
        resultStream << result.dump(2) << "\n";
        spdlog::info("Wrote {}", path);
    } else {
        spdlog::warn("Could not write result document to {}", path);
    }
}

PluginDescriptor buildDescriptor() {
    PluginDescriptor descriptor{"ilda", "Interfacial Line Defect Analysis (ILDA)", optionsMeta(bindings)};
    descriptor.defaultThreads = DefaultThreads::PhysicalCores;
    descriptor.extraModes = IldaSelftest::run;
    return descriptor;
}

json run(const IldaOptions& options, const LammpsParser::Frame& frame, const LammpsParser::Frame*, const std::string& outputBase) {
    ILDAService analyzer;
    analyzer.setOptions(options);
    const json result = analyzer.run(frame, outputBase);
    if(!result.value("is_failed", false)) {
        writeResultDocument(result, outputBase + "_result.json");
    }
    return result;
}

}

VOLT_PLUGIN_MAIN(IldaPlugin::buildDescriptor(), bindings, IldaPlugin::logResolvedParameters, IldaPlugin::run)
