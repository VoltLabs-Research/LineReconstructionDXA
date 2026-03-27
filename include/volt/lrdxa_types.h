#pragma once

#include <volt/core/volt.h>
#include <volt/helpers/half_edge_mesh.h>
#include <volt/helpers/cluster_vector.h>
#include <volt/structures/crystal_structure_types.h>

namespace Volt::DXA {

struct LineReconstructionSegment {
    Point3 position1;
    Point3 position2;
    ClusterVector burgersVector;
    int clusterId = 0;
    int structureType = LATTICE_OTHER;
    int stage = 0;
};

struct LineReconstructionUnassignedEdge {
    Point3 position1;
    Point3 position2;
    int atom1 = -1;
    int atom2 = -1;
    int stage = 0;
};

struct LineReconstructionDislocationLine {
    std::vector<Point3> points;
    ClusterVector burgersVector;
    int clusterId = 0;
    int structureType = LATTICE_OTHER;
    bool isClosed = false;
    int dislocationTypeId = 0;
};

struct LineReconstructionInterfaceMeshVertex {};

struct LineReconstructionInterfaceMeshFace {
    int region = 0;
};

struct LineReconstructionInterfaceMeshEdge {};

using LineReconstructionInterfaceMesh = HalfEdgeMesh<LineReconstructionInterfaceMeshEdge, LineReconstructionInterfaceMeshFace, LineReconstructionInterfaceMeshVertex>;

}
