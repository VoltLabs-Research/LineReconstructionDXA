#include <volt/cli/common.h>
#include <volt/lrdxa_service.h>
#include <volt/structures/crystal_structure_types.h>
#include <tbb/info.h>
#include <algorithm>
#include <fstream>
#include <set>
#include <string>

using namespace Volt;
using namespace Volt::CLI;

namespace {

LatticeStructureType parseCrystalStructure(const std::string& str) {
    if(str == "FCC") return LATTICE_FCC;
    if(str == "BCC") return LATTICE_BCC;
    if(str == "HCP") return LATTICE_HCP;
    if(str == "SC") return LATTICE_SC;
    if(str == "CUBIC_DIAMOND") return LATTICE_CUBIC_DIAMOND;
    if(str == "HEX_DIAMOND") return LATTICE_HEX_DIAMOND;
    spdlog::warn("Unknown crystal structure '{}' , defaulting to FCC.", str);
    return LATTICE_FCC;
}

StructureAnalysis::Mode parseIdentificationMode(const std::string& str) {
    if(str == "CNA") return StructureAnalysis::Mode::CNA;
    if(str == "PTM") return StructureAnalysis::Mode::PTM;
    if(str == "DIAMOND") return StructureAnalysis::Mode::DIAMOND;
    spdlog::warn("Unknown identification mode '{}' , defaulting to CNA.", str);
    return StructureAnalysis::Mode::CNA;
}

void showUsage(const std::string& name) {
    printUsageHeader(name, "Volt - Line Reconstruction DXA");
    std::cerr
        << "  --crystalStructure <type>              Reference crystal structure. (BCC|FCC|HCP|CUBIC_DIAMOND|HEX_DIAMOND|SC) [default: FCC]\n"
        << "  --identificationMode <mode>            Structure identification mode. (CNA|PTM|DIAMOND) [default: CNA]\n"
        << "  --rmsd <float>                         RMSD threshold for PTM. [default: 0.12]\n"
        << "  --minClusterSize <int>                 Minimum crystalline cluster size kept before dissolution. [default: 50]\n"
        << "  --crystalPathSteps <int>               Maximum crystal-path steps used for edge vectors. [default: 4]\n"
        << "  --tessellationGhostLayerScale <float>  Ghost-layer scale relative to neighbor distance. [default: 3.5]\n"
        << "  --alphaScale <float>                   Alpha threshold scale relative to neighbor distance. [default: 3.5]\n"
        << "  --smoothingIterations <int>            Taubin smoothing iterations for reconstructed lines. [default: 0]\n"
        << "  --linePointInterval <float>            Line coarsening interval. [default: 1.2]\n"
        << "  --onlyPerfectDislocations <bool>       Restrict upstream structure analysis. [default: false]\n"
        << "  --structureIdentificationOnly <bool>   Export only structure identification outputs. [default: false]\n"
        << "  --threads <int>                        Max worker threads (TBB/OMP). [default: auto capped to physical cores]\n";
    printHelpOption();
}

}  // namespace

int main(int argc, char* argv[]) {
    if(argc < 2) {
        showUsage(argv[0]);
        return 1;
    }

    std::string filename;
    std::string outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);
    if(hasOption(opts, "--help") || filename.empty()) {
        showUsage(argv[0]);
        return filename.empty() ? 1 : 0;
    }

    if(!hasOption(opts, "--threads")) {
        const int maxAvailableThreads = static_cast<int>(oneapi::tbb::info::default_concurrency());
        int physicalCores = 0;
        std::ifstream cpuinfo("/proc/cpuinfo");
        if(cpuinfo.is_open()) {
            std::set<std::pair<int, int>> physicalCoreIds;
            int fallbackCpuCores = 0;
            int physicalId = -1;
            int coreId = -1;
            std::string line;
            while(std::getline(cpuinfo, line)) {
                if(line.empty()) {
                    if(physicalId >= 0 && coreId >= 0) {
                        physicalCoreIds.emplace(physicalId, coreId);
                    }
                    physicalId = -1;
                    coreId = -1;
                    continue;
                }
                int parsedValue = 0;
                const auto separator = line.find(':');
                if(separator != std::string::npos) {
                    try {
                        parsedValue = std::stoi(line.substr(separator + 1));
                    } catch(const std::exception&) {
                        parsedValue = 0;
                    }
                }
                if(line.rfind("physical id", 0) == 0) {
                    physicalId = parsedValue;
                } else if(line.rfind("core id", 0) == 0) {
                    coreId = parsedValue;
                } else if(line.rfind("cpu cores", 0) == 0) {
                    fallbackCpuCores = std::max(fallbackCpuCores, parsedValue);
                }
            }
            if(physicalId >= 0 && coreId >= 0) {
                physicalCoreIds.emplace(physicalId, coreId);
            }
            physicalCores = physicalCoreIds.empty() ? fallbackCpuCores : static_cast<int>(physicalCoreIds.size());
        }
        opts["--threads"] = std::to_string(std::max(1, std::min(maxAvailableThreads, physicalCores > 0 ? physicalCores : maxAvailableThreads)));
    }

    auto parallel = initParallelism(opts, false);
    initLogging("line-reconstruction-dxa", parallel.threads);

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)) return 1;

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    LineReconstructionDXAService analyzer;
    analyzer.setInputFilename(filename);
    analyzer.setInputCrystalStructure(parseCrystalStructure(getString(opts, "--crystalStructure", "FCC")));
    analyzer.setIdentificationMode(parseIdentificationMode(getString(opts, "--identificationMode", "CNA")));
    analyzer.setRmsd(getDouble(opts, "--rmsd", 0.12f));
    analyzer.setMinClusterSize(getInt(opts, "--minClusterSize", 50));
    analyzer.setCrystalPathSteps(getInt(opts, "--crystalPathSteps", 4));
    analyzer.setTessellationGhostLayerScale(getDouble(opts, "--tessellationGhostLayerScale", 3.5));
    analyzer.setAlphaScale(getDouble(opts, "--alphaScale", 3.5));
    analyzer.setSmoothingIterations(getInt(opts, "--smoothingIterations", 0));
    analyzer.setLinePointInterval(getDouble(opts, "--linePointInterval", 1.2));
    analyzer.setOnlyPerfectDislocations(getBool(opts, "--onlyPerfectDislocations"));
    analyzer.setStructureIdentificationOnly(getBool(opts, "--structureIdentificationOnly"));

    spdlog::info("Starting line reconstruction DXA analysis...");
    json result = analyzer.compute(frame, outputBase);
    if(result.value("is_failed", false)) {
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }
    spdlog::info("Analysis completed successfully.");
    return 0;
}
