#include <volt/lrdxa_service.h>
#include <volt/lrdxa_pipeline.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/analysis_context.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string_view>

namespace Volt {

using namespace Volt::Particles;

namespace {

json readMsgpackFile(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if(!input) {
        throw std::runtime_error("Failed to open MessagePack file: " + path.string());
    }

    std::vector<std::uint8_t> bytes(
        (std::istreambuf_iterator<char>(input)),
        std::istreambuf_iterator<char>()
    );
    return json::from_msgpack(bytes);
}

void copyIfExists(const std::filesystem::path& source, const std::filesystem::path& destination) {
    if(std::filesystem::exists(source)) {
        std::filesystem::copy_file(source, destination, std::filesystem::copy_options::overwrite_existing);
    }
}

std::string shellEscape(const std::string& value) {
    std::string escaped = "'";
    for(char ch : value) {
        if(ch == '\'') {
            escaped += "'\\''";
        } else {
            escaped += ch;
        }
    }
    escaped += "'";
    return escaped;
}

std::filesystem::path detectCurrentExecutable() {
    std::error_code ec;
    const auto exe = std::filesystem::read_symlink("/proc/self/exe", ec);
    return ec ? std::filesystem::path{} : exe;
}

std::filesystem::path findOpenDxaExecutable() {
    if(const char* envBinary = std::getenv("VOLTLABS_OPENDXA_BINARY")) {
        std::filesystem::path candidate(envBinary);
        if(std::filesystem::exists(candidate)) {
            return candidate;
        }
    }

    const auto selfExecutable = detectCurrentExecutable();
    if(!selfExecutable.empty()) {
        const auto candidate = selfExecutable.parent_path().parent_path().parent_path() / "OpenDXA" / "build" / "Release" / "opendxa";
        if(std::filesystem::exists(candidate)) {
            return candidate;
        }
    }

    const auto sourcePath = std::filesystem::path(__FILE__);
    const auto repoRoot = sourcePath.parent_path().parent_path().parent_path();
    const auto repoCandidate = repoRoot / "OpenDXA" / "build" / "Release" / "opendxa";
    if(std::filesystem::exists(repoCandidate)) {
        return repoCandidate;
    }

    return std::filesystem::path("opendxa");
}

const char* crystalStructureName(LatticeStructureType structureType) {
    switch(structureType) {
        case LATTICE_FCC: return "FCC";
        case LATTICE_BCC: return "BCC";
        case LATTICE_HCP: return "HCP";
        case LATTICE_SC: return "SC";
        case LATTICE_CUBIC_DIAMOND: return "CUBIC_DIAMOND";
        case LATTICE_HEX_DIAMOND: return "HEX_DIAMOND";
        default: return "FCC";
    }
}

bool shouldUseOpenDxaFallback(LatticeStructureType structureType) {
    return structureType == LATTICE_SC || structureType == LATTICE_BCC;
}

json ensureBurgersMetadata(const json& entry) {
    json burgers = entry.value("burgers", json::object());
    const json vector = entry.value("burgers_vector", json::array({0.0, 0.0, 0.0}));
    burgers["vector"] = burgers.value("vector", vector);
    burgers["vector_local"] = burgers.value("vector_local", entry.value("burgers_vector_local", vector));
    burgers["vector_global"] = burgers.value("vector_global", entry.value("burgers_vector_global", vector));
    burgers["magnitude"] = burgers.value("magnitude", entry.value("magnitude", 0.0));
    return burgers;
}

json adaptOpenDxaSegment(const json& source, int structureType) {
    json segment = source;
    segment["cluster_id"] = source.value("cluster_id", 0);
    segment["structure_type"] = source.value("structure_type", structureType);
    segment["stage"] = source.value("stage", 0);
    segment["burgers"] = ensureBurgersMetadata(source);
    return segment;
}

json adaptOpenDxaLine(const json& source, int structureType) {
    json line = source;
    line["line_id"] = source.value("segment_id", source.value("line_id", 0));
    line["is_closed"] = source.value("is_closed", false);
    line["cluster_id"] = source.value("cluster_id", 0);
    line["structure_type"] = source.value("structure_type", structureType);
    line["dislocation_type_id"] = source.value("dislocation_type_id", 0);
    line["burgers"] = ensureBurgersMetadata(source);
    return line;
}

void writeOpenDxaFallbackOutputs(const std::filesystem::path& tempBase, const std::string& outputFile, int structureType) {
    const auto dislocationsPath = tempBase.string() + "_dislocations.msgpack";
    const json openDxaDislocations = readMsgpackFile(dislocationsPath);
    const json sourceSegments = openDxaDislocations.at("sub_listings").value("dislocation_segments", json::array());

    json lineArray = json::array();
    json segmentArray = json::array();
    for(const auto& sourceSegment : sourceSegments) {
        segmentArray.push_back(adaptOpenDxaSegment(sourceSegment, structureType));
        lineArray.push_back(adaptOpenDxaLine(sourceSegment, structureType));
    }

    json dislocationsJson = {
        {"export", openDxaDislocations.value("export", json::object())},
        {"main_listing", {
            {"dislocations", lineArray.size()},
            {"total_length", openDxaDislocations.at("main_listing").value("total_length", 0.0)},
            {"total_points", openDxaDislocations.at("main_listing").value("total_points", 0)}
        }},
        {"sub_listings", {{"dislocations", lineArray}}}
    };
    JsonUtils::writeJsonMsgpackToFile(dislocationsJson, outputFile + "_dislocations.msgpack", false);

    json segmentsJson = {
        {"export", openDxaDislocations.value("export", json::object())},
        {"main_listing", {
            {"dislocation_segments", segmentArray.size()},
            {"dislocations", segmentArray.size()}
        }},
        {"sub_listings", {{"dislocation_segments", segmentArray}}}
    };
    JsonUtils::writeJsonMsgpackToFile(segmentsJson, outputFile + "_dislocation_segments.msgpack", false);

    json emptyUnassignedEdges = {
        {"export", {{"LineReconstructionUnassignedEdgeExporter", {{"unassigned_edges", json::array()}}}}},
        {"main_listing", {{"unassigned_edges", 0}}},
        {"sub_listings", {{"unassigned_edges", json::array()}}}
    };
    JsonUtils::writeJsonMsgpackToFile(emptyUnassignedEdges, outputFile + "_unassigned_edges.msgpack", false);

    copyIfExists(tempBase.string() + "_defect_mesh.msgpack", outputFile + "_interface_mesh.msgpack");
    copyIfExists(tempBase.string() + "_simulation_cell.msgpack", outputFile + "_simulation_cell.msgpack");
    copyIfExists(tempBase.string() + "_atoms.msgpack", outputFile + "_atoms.msgpack");
    copyIfExists(tempBase.string() + "_structure_analysis_stats.msgpack", outputFile + "_structure_analysis_stats.msgpack");
    copyIfExists(tempBase.string() + "_ptm_data.msgpack", outputFile + "_ptm_data.msgpack");
}

json runOpenDxaFallback(const std::string& inputFilename, const std::string& outputFile, const LineReconstructionDXAOptions& options) {
    const auto tempRoot = std::filesystem::temp_directory_path();
    const auto tempTag = std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    const auto tempDir = tempRoot / ("lrdxa-sc-" + tempTag);
    const auto tempBase = tempDir / "opendxa_sc";
    std::filesystem::create_directories(tempDir);

    const auto opendxaBinary = findOpenDxaExecutable();

    std::ostringstream command;
    command
        << shellEscape(opendxaBinary.string()) << ' '
        << shellEscape(inputFilename) << ' '
        << shellEscape(tempBase.string()) << ' '
        << "--identificationMode " << shellEscape(options.identificationMode == StructureAnalysis::Mode::PTM ? "PTM" : "CNA") << ' '
        << "--crystalStructure " << shellEscape(crystalStructureName(options.inputCrystalStructure)) << ' '
        << "--rmsd " << options.rmsd << ' '
        << "--lineSmoothingLevel " << options.smoothingIterations << ' '
        << "--linePointInterval " << options.linePointInterval;

    if(options.onlyPerfectDislocations) {
        command << " --onlyPerfectDislocations true";
    }

    if(std::system(command.str().c_str()) != 0) {
        std::filesystem::remove_all(tempDir);
        return AnalysisResult::failure("OpenDXA fallback failed");
    }

    const std::filesystem::path finalBase = outputFile.empty() ? tempBase : std::filesystem::path(outputFile);
    if(!outputFile.empty()) {
        writeOpenDxaFallbackOutputs(tempBase, outputFile, options.inputCrystalStructure);
    }

    const json adaptedDislocations = readMsgpackFile(finalBase.string() + "_dislocations.msgpack");
    const auto& mainListing = adaptedDislocations.at("main_listing");
    json result = {
        {"is_failed", false},
        {"algorithm", "line-reconstruction-dxa"},
        {"backend", "opendxa_fallback"},
        {"lines", mainListing.value("dislocations", 0)},
        {"segments", mainListing.value("dislocations", 0)},
        {"total_length", mainListing.value("total_length", 0.0)}
    };
    std::filesystem::remove_all(tempDir);
    return result;
}

}  // namespace

LineReconstructionDXA::LineReconstructionDXA()
    : _inputCrystalStructure(LATTICE_FCC)
    , _inputFilename()
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

void LineReconstructionDXA::setInputFilename(std::string filename) {
    _inputFilename = std::move(filename);
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

    if(shouldUseOpenDxaFallback(options.inputCrystalStructure)) {
        result = runOpenDxaFallback(_inputFilename, outputFile, options);
        markStage("opendxa_fallback_analysis");
        result["total_time"] = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - startTime
        ).count();
        result["stage_metrics"] = std::move(stageMetrics);
        return result;
    }

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
