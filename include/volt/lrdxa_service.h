#pragma once

#include <volt/core/volt.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/lrdxa_options.h>
#include <volt/lrdxa_json_exporter.h>
#include <nlohmann/json.hpp>
#include <string>

namespace Volt {

using json = nlohmann::json;

class LineReconstructionDXA {
public:
    LineReconstructionDXA();

    void setInputCrystalStructure(LatticeStructureType structure);
    void setInputFilename(std::string filename);
    void setIdentificationMode(StructureAnalysis::Mode identificationMode);
    void setRmsd(float rmsd);
    void setStructureIdentificationOnly(bool structureIdentificationOnly);
    void setOnlyPerfectDislocations(bool flag);
    void setMinClusterSize(int minClusterSize);
    void setCrystalPathSteps(int crystalPathSteps);
    void setTessellationGhostLayerScale(double tessellationGhostLayerScale);
    void setAlphaScale(double alphaScale);
    void setSmoothingIterations(int smoothingIterations);
    void setLinePointInterval(double linePointInterval);

    json compute(const LammpsParser::Frame& frame, const std::string& outputFile = "");

private:
    LineReconstructionDXAOptions buildOptions() const;

private:
    LatticeStructureType _inputCrystalStructure;
    std::string _inputFilename;
    StructureAnalysis::Mode _identificationMode;
    float _rmsd;
    bool _structureIdentificationOnly;
    bool _onlyPerfectDislocations;
    int _minClusterSize;
    int _crystalPathSteps;
    double _tessellationGhostLayerScale;
    double _alphaScale;
    int _smoothingIterations;
    double _linePointInterval;

    mutable LineReconstructionJsonExporter _jsonExporter;
};

}  // namespace Volt

namespace Volt {
using LineReconstructionDXAService = LineReconstructionDXA;
}
