#include <volt/lrdxa_service.h>
#include <volt/lrdxa_pipeline.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/analysis_context.h>
#include <spdlog/spdlog.h>
#include <chrono>
#include <string_view>

namespace Volt {

using namespace Volt::Particles;

LineReconstructionDXA::LineReconstructionDXA()
    : _inputCrystalStructure(LATTICE_FCC)
    , _identificationMode(StructureAnalysis::Mode::CNA)
    , _rmsd(0.12f)
    , _structureIdentificationOnly(false)
    , _onlyPerfectDislocations(false)
    , _minClusterSize(50)
    , _crystalPathSteps(4)
    , _tessellationGhostLayerScale(3.5)
    , _alphaScale(3.5)
    , _smoothingIterations(0)
    , _linePointInterval(1.2) {
}

void LineReconstructionDXA::setInputCrystalStructure(LatticeStructureType structure) {
    _inputCrystalStructure = structure;
}

void LineReconstructionDXA::setIdentificationMode(StructureAnalysis::Mode identificationMode) {
    _identificationMode = identificationMode;
}

void LineReconstructionDXA::setRmsd(float rmsd) {
    _rmsd = rmsd;
}

void LineReconstructionDXA::setStructureIdentificationOnly(bool structureIdentificationOnly) {
    _structureIdentificationOnly = structureIdentificationOnly;
}

void LineReconstructionDXA::setOnlyPerfectDislocations(bool flag) {
    _onlyPerfectDislocations = flag;
}

void LineReconstructionDXA::setMinClusterSize(int minClusterSize) {
    _minClusterSize = minClusterSize;
}

void LineReconstructionDXA::setCrystalPathSteps(int crystalPathSteps) {
    _crystalPathSteps = crystalPathSteps;
}

void LineReconstructionDXA::setTessellationGhostLayerScale(double tessellationGhostLayerScale) {
    _tessellationGhostLayerScale = tessellationGhostLayerScale;
}

void LineReconstructionDXA::setAlphaScale(double alphaScale) {
    _alphaScale = alphaScale;
}

void LineReconstructionDXA::setSmoothingIterations(int smoothingIterations) {
    _smoothingIterations = smoothingIterations;
}

void LineReconstructionDXA::setLinePointInterval(double linePointInterval) {
    _linePointInterval = linePointInterval;
}

LineReconstructionDXAOptions LineReconstructionDXA::buildOptions() const {
    return LineReconstructionDXAOptions{
        .inputCrystalStructure = _inputCrystalStructure,
        .identificationMode = _identificationMode,
        .rmsd = _rmsd,
        .structureIdentificationOnly = _structureIdentificationOnly,
        .onlyPerfectDislocations = _onlyPerfectDislocations,
        .minClusterSize = _minClusterSize,
        .crystalPathSteps = _crystalPathSteps,
        .tessellationGhostLayerScale = _tessellationGhostLayerScale,
        .alphaScale = _alphaScale,
        .smoothingIterations = _smoothingIterations,
        .linePointInterval = _linePointInterval,
    };
}

json LineReconstructionDXA::compute(const LammpsParser::Frame& frame, const std::string& outputFile) {
    const auto startTime = std::chrono::high_resolution_clock::now();
    auto stageStart = startTime;
    const LineReconstructionDXAOptions options = buildOptions();

    json result;
    json stageMetrics = json::array();
    auto markStage = [&](std::string_view name) {
        const auto now = std::chrono::high_resolution_clock::now();
        const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - stageStart).count();
        stageMetrics.push_back({{"stage", name}, {"elapsed_ms", elapsed}});
        stageStart = now;
    };

    if(frame.natoms <= 0) {
        return AnalysisResult::failure("Invalid number of atoms: " + std::to_string(frame.natoms));
    }
    if(frame.positions.empty()) {
        return AnalysisResult::failure("No position data available");
    }
    if(!FrameAdapter::validateSimulationCell(frame.simulationCell)) {
        return AnalysisResult::failure("Invalid simulation cell");
    }

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions) {
        return AnalysisResult::failure("Failed to create position property");
    }
    markStage("create_position_property");

    std::vector<Matrix3> preferredOrientations{Matrix3::Identity()};
    auto structureTypes = std::make_unique<ParticleProperty>(frame.natoms, DataType::Int, 1, 0, true);
    AnalysisContext context(
        positions.get(),
        frame.simulationCell,
        options.inputCrystalStructure,
        nullptr,
        structureTypes.get(),
        std::move(preferredOrientations)
    );

    StructureAnalysis structureAnalysis(
        context,
        !options.onlyPerfectDislocations,
        options.identificationMode,
        options.rmsd
    );
    markStage("structure_analysis_setup");

    structureAnalysis.identifyStructures();
    markStage("identify_structures");

    if(!outputFile.empty() && options.identificationMode == StructureAnalysis::Mode::PTM) {
        _jsonExporter.exportPTMData(structureAnalysis.context(), frame.ids, outputFile);
        markStage("export_ptm");
    }

    if(!outputFile.empty()) {
        _jsonExporter.exportForStructureIdentification(frame, structureAnalysis, outputFile);
        markStage("stream_atoms_msgpack");
    }

    if(options.structureIdentificationOnly && !outputFile.empty()) {
        result["is_failed"] = false;
        result["stage_metrics"] = std::move(stageMetrics);
        result["total_time"] = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - startTime
        ).count();
        return result;
    }

    result = DXA::LineReconstructionDXAPipeline::run(
        frame,
        outputFile,
        options,
        structureAnalysis,
        context,
        _jsonExporter,
        markStage
    );

    result["total_time"] = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - startTime
    ).count();
    result["stage_metrics"] = std::move(stageMetrics);
    return result;
}

}  // namespace Volt
