#pragma once

#include <volt/ilda_types.h>
#include <volt/core/simulation_cell.h>

#include <string>
#include <vector>

namespace Volt {

class IldaJsonExporter {
public:
    explicit IldaJsonExporter() = default;

    void writeBondsParquet(
        const std::vector<IldaSegment>& segments,
        const SimulationCell& cell,
        const std::string& filePath
    );

    void writeSummaryParquet(
        const std::vector<DisconnectionMode>& modes,
        double totalLength,
        int Nb,
        const std::string& filePath
    );
};

}
