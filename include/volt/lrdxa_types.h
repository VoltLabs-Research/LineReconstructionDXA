#pragma once

#include <volt/core/volt.h>
#include <volt/geometry/half_edge_mesh.h>
#include <volt/structures/cluster_vector.h>
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

// A stitched dislocation line: a connected polyline of points with a single Burgers vector.
struct LineReconstructionDislocationLine {
    std::vector<Point3> points;
    ClusterVector burgersVector;
    int clusterId = 0;
    int structureType = LATTICE_OTHER;
    bool isClosed = false;  // true if the line forms a closed loop
    int dislocationTypeId = 0;  // FCC type: 0=unknown, 1=perfect, 2=Shockley, 3=Frank, 4=stair-rod, 5=Hirth
};

struct LineReconstructionInterfaceMeshVertex {};

struct LineReconstructionInterfaceMeshFace {
    int region = 0;
};

struct LineReconstructionInterfaceMeshEdge {};

using LineReconstructionInterfaceMesh = HalfEdgeMesh<LineReconstructionInterfaceMeshEdge, LineReconstructionInterfaceMeshFace, LineReconstructionInterfaceMeshVertex>;

}  // namespace Volt::DXA
