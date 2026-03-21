#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/structures/crystal_structure_types.h>

namespace Volt {

struct LineReconstructionDXAOptions {
    LatticeStructureType inputCrystalStructure = LATTICE_FCC;
    StructureAnalysis::Mode identificationMode = StructureAnalysis::Mode::CNA;
    float rmsd = 0.12f;
    bool structureIdentificationOnly = false;
    bool onlyPerfectDislocations = false;
    int minClusterSize = 50;
    int crystalPathSteps = 4;
    double tessellationGhostLayerScale = 3.5;
    double alphaScale = 3.5;
    int smoothingIterations = 0;
    double linePointInterval = 1.2;
};

} 
