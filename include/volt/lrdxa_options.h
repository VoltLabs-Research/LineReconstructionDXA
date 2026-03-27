#pragma once

namespace Volt {

struct LineReconstructionDXAOptions {
    int crystalPathSteps = 4;
    double tessellationGhostLayerScale = 3.5;
    double alphaScale = 3.5;
    int smoothingIterations = 0;
    double linePointInterval = 1.2;
};

} 
