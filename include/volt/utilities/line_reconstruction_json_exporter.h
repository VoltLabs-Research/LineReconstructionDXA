#pragma once

#include <nlohmann/json.hpp>
#include <string>
#include <vector>
#include <volt/analysis/analysis_context.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/lammps_parser.h>
#include <volt/dxa/line_reconstruction_dxa/line_reconstruction_dxa_types.h>
#include <volt/math/lin_alg.h>

namespace Volt {

using json = nlohmann::json;

class LineReconstructionJsonExporter {
public:
    explicit LineReconstructionJsonExporter() = default;

    void exportForStructureIdentification(
        const LammpsParser::Frame& frame,
        const StructureAnalysis& structureAnalysis,
        const std::string& outputFilename
    );

    void exportPTMData(
        const AnalysisContext& context,
        const std::vector<int>& ids,
        const std::string& outputFilename
    );

    json getExtendedSimulationCellInfo(const SimulationCell& cell);

    void writeLineSegmentsMsgpackToFile(
        const std::vector<DXA::LineReconstructionSegment>& segments,
        const SimulationCell& simulationCell,
        const std::string& filePath
    );

    void writeDislocationLinesMsgpackToFile(
        const std::vector<DXA::LineReconstructionDislocationLine>& lines,
        const SimulationCell& simulationCell,
        const std::string& filePath
    );

    void writeUnassignedEdgesMsgpackToFile(
        const std::vector<DXA::LineReconstructionUnassignedEdge>& edges,
        const std::string& filePath
    );

    void writeInterfaceMeshMsgpackToFile(
        const DXA::LineReconstructionInterfaceMesh& mesh,
        const StructureAnalysis& structureAnalysis,
        const std::string& filePath
    );

private:
    json exportLineSegmentsToJson(
        const std::vector<DXA::LineReconstructionSegment>& segments,
        const SimulationCell* simulationCell
    );
    json exportDislocationLinesToJson(
        const std::vector<DXA::LineReconstructionDislocationLine>& lines,
        const SimulationCell* simulationCell
    );
    json exportUnassignedEdgesToJson(const std::vector<DXA::LineReconstructionUnassignedEdge>& edges);

    template <typename MeshType>
    json getMeshData(const MeshType& mesh, const StructureAnalysis& structureAnalysis);

    json vectorToJson(const Vector3& vector);
    json affineTransformationToJson(const AffineTransformation& transform);
    json simulationCellToJson(const SimulationCell& cell);
    double calculateAngle(const Vector3& a, const Vector3& b);
};

}  // namespace Volt
