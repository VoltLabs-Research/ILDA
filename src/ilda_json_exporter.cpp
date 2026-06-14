#include <volt/ilda_json_exporter.h>
#include <volt/utilities/parquet_bond_writer.h>
#include <volt/utilities/json_utils.h>

#include <vector>

#include <spdlog/spdlog.h>

namespace Volt {

void IldaJsonExporter::writeBondsParquet(
    const std::vector<IldaSegment>& segments,
    const SimulationCell& /*cell*/,
    const std::string& filePath
) {
    std::vector<Bond> bonds;
    bonds.reserve(segments.size());
    for(std::size_t i = 0; i < segments.size(); ++i) {
        const IldaSegment& s = segments[i];
        Bond b;
        b.id = static_cast<std::int64_t>(i);
        b.atomA = s.atom1;
        b.atomB = s.atom2;
        b.pbcShift = { {0, 0, 0} };  // PARITY: pos1/pos2 already carry the image shift
        b.posA = s.pos1;
        b.posB = s.pos2;
        b.distance = (s.pos2 - s.pos1).length();
        bonds.push_back(b);
    }

    streamBondsToParquet(
        filePath,
        bonds,
        [&](ColumnarBondWriter& writer, std::size_t i) {
            const IldaSegment& s = segments[i];
            writer.field("burgers_vector", std::vector<double>{ s.burgers.x(), s.burgers.y(), s.burgers.z() });
            writer.field("terrace_plane", std::vector<double>{ s.plane.x(), s.plane.y(), s.plane.z() });
            writer.field("color", std::vector<double>{ s.color.x(), s.color.y(), s.color.z() });
            writer.field("step_height", s.step);
            writer.field("ideal_step_height_a", s.stepIdealA);
            writer.field("ideal_step_height_b", s.stepIdealB);
            writer.field("burgers_id", static_cast<std::int64_t>(s.bid));
            writer.field("magnitude", s.bnorm);
            writer.field("length", (s.pos2 - s.pos1).length());
        }
    );

    spdlog::info("Wrote {} ILDA line segments -> {}", bonds.size(), filePath);
}

void IldaJsonExporter::writeSummaryParquet(
    const std::vector<DisconnectionMode>& modes,
    double totalLength,
    int Nb,
    const std::string& filePath
) {
    json rows = json::array();
    int variantTotal = 0;
    for(const DisconnectionMode& mode : modes) {
        for(std::size_t j = 0; j < mode.variants.size(); ++j) {
            const IldaVariant& v = mode.variants[j];
            // PARITY: column order must mirror the Python listing.
            rows.push_back({
                {"ID", mode.id},
                {"Variant", static_cast<int>(j) + 1},
                {"Length", v.length},
                {"Avg b", v.bAvgNorm},
                {"Avg bx", v.bAvg.x()},
                {"Avg by", v.bAvg.y()},
                {"Avg bz", v.bAvg.z()},
                {"Std bx", v.bStdX},
                {"Std by", v.bStdY},
                {"Std bz", v.bStdZ},
                {"Ideal hA", v.hIdealA},
                {"Ideal hB", v.hIdealB},
                {"Avg h", v.hAvg}
            });
            ++variantTotal;
        }
    }

    json payload;
    payload["main_listing"] = {
        {"disconnection_modes", Nb},
        {"variants", variantTotal},
        {"total_line_length", totalLength}
    };
    payload["sub_listings"] = {
        {"disconnection_modes", rows}
    };

    JsonUtils::writeJsonToParquet(payload, filePath);
    spdlog::info("Wrote ILDA summary ({} modes, {} variants, L={:.3f}) -> {}",
                 Nb, variantTotal, totalLength, filePath);
}

}
