# LineReconstructionDXA

`LineReconstructionDXA` reconstructs dislocation lines from an upstream cluster package and DXA-style geometric parameters.

## CLI

Usage:

```bash
line-reconstruction-dxa <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--clusters-table <path>` | Yes | Path to `*_clusters.table` exported upstream. | |
| `--clusters-transitions <path>` | Yes | Path to `*_cluster_transitions.table` exported upstream. | |
| `--crystalStructure <type>` | No | Reference crystal structure: `BCC`, `FCC`, `HCP`, `CUBIC_DIAMOND`, `HEX_DIAMOND`, `SC`. | `FCC` |
| `--crystalPathSteps <int>` | No | Maximum crystal-path steps used for edge vectors. | `4` |
| `--tessellationGhostLayerScale <float>` | No | Ghost-layer scale relative to neighbor distance. | `3.5` |
| `--alphaScale <float>` | No | Alpha threshold scale relative to neighbor distance. | `3.5` |
| `--smoothingIterations <int>` | No | Taubin smoothing iterations for reconstructed lines. | `0` |
| `--linePointInterval <float>` | No | Line coarsening interval. | `1.2` |
| `--threads <int>` | No | Maximum worker threads. | auto capped to physical cores |
| `--help` | No | Print CLI help. | |

## Build With CoreToolkit

```bash
cd /path/to/voltlabs-ecosystem/tools/CoreToolkit
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/StructureIdentification
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/CommonNeighborAnalysis
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/PolyhedralTemplateMatching
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/OpenDXA
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/LineReconstructionDXA
conan create . -nr
```
