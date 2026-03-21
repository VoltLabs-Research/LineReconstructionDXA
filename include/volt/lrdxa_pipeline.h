#pragma once

#include <volt/core/lammps_parser.h>
#include <volt/lrdxa_options.h>
#include <nlohmann/json.hpp>
#include <functional>
#include <string>
#include <string_view>

namespace Volt {

using json = nlohmann::json;

class AnalysisContext;
class StructureAnalysis;
class LineReconstructionJsonExporter;

namespace DXA {

class LineReconstructionDXAPipeline {
public:
    static json run(
        const LammpsParser::Frame& frame,
        const std::string& outputFile,
        const LineReconstructionDXAOptions& options,
        StructureAnalysis& reconstructionStructureAnalysis,
        AnalysisContext& reconstructionContext,
        LineReconstructionJsonExporter& jsonExporter,
        const std::function<void(std::string_view)>& markStage
    );
};

}  // namespace DXA

}  // namespace Volt
