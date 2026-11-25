#pragma once

#include <volt/core/volt.h>
#include <volt/ilda_options.h>
#include <volt/plugin/option_binding.h>

#include <nlohmann/json.hpp>

#include <string>

namespace Volt::IldaPlugin {

using json = nlohmann::json;
using Binding = Plugin::OptionBinding<IldaOptions>;

Plugin::PluginDescriptor buildDescriptor();

json run(
    const IldaOptions& options,
    const LammpsParser::Frame& frame,
    const LammpsParser::Frame* reference,
    const std::string& outputBase
);

// An option that sets one field of the options, with the field's own default.
Binding pathOption(const char* name, const char* help, std::string IldaOptions::* field);
Binding vectorOption(const char* name, const char* help, Vector3 IldaOptions::* field);
Binding matrixOption(const char* name, const char* help, Matrix3 IldaOptions::* field);

// A vector option given as JSON, [x, y, z] or {"x": .., "y": .., "z": ..}; the fallback when it is missing or cannot
// be read.
Vector3 parseVector3Option(const Plugin::OptsMap& opts, const std::string& key, const Vector3& fallback);
// A matrix option given as JSON, nine numbers row by row or three rows as arrays or objects; the fallback when it is
// missing or cannot be read.
Matrix3 parseMatrix3Option(const Plugin::OptsMap& opts, const std::string& key, const Matrix3& fallback);

// The values the options resolved to, the run's "Parameters:" line. A negative number such as "--typeA -1" reaches
// parseArgs as "true", and the option then keeps its default, so the raw flags would not say what the run used.
void logResolvedParameters(const IldaOptions& options);
void writeResultDocument(const json& result, const std::string& path);

}
