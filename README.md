# Interfacial Line Defect Analysis (ILDA)

Extracts interfacial dislocations and disconnections at a semi-coherent interface, inserting a line-segment representation annotated with Burgers vector and step height per segment (a DXA-style line-defect analysis).

## Install

```bash
vpm install @voltlabs/ilda
```

## CLI

```bash
ilda <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--grain-atoms` | yes | — | Per-atom annotated parquet (grain id, structure type, orientation) from the upstream grain-segmentation step. |
| `--grains` | yes | — | Per-grain table (structure type, orientation) from the upstream grain-segmentation step. |
| `--atomA` | no | `0` | Atom A particle identifier. |
| `--atomB` | no | `0` | Atom B particle identifier. |
| `--aA` | no | `1.0` | Lattice constant a (Grain A). |
| `--cA` | no | `0.0` | Lattice constant c (Grain A). |
| `--aB` | no | `1.0` | Lattice constant a (Grain B). |
| `--cB` | no | `0.0` | Lattice constant c (Grain B). |
| `--typeA` | no | `-1` | Structure type A (`-1` = derive). |
| `--typeB` | no | `-1` | Structure type B (`-1` = derive). |
| `--Rsphere` | no | `10.0` | Probe sphere radius. |
| `--htol` | no | `0.5` | Step height tolerance. |
| `--btol` | no | `0.01` | Burgers length tolerance. |
| `--angtol` | no | `5.0` | Burgers angular tolerance (deg). |
| `--distF` | no | `10.0` | Interface skin distance. |
| `--cis_tol` | no | `0.0` | Co-incidence site tolerance. |
| `--rmsd` | no | `0.1` | PTM RMSD cutoff. |
| `--estimateF` | no | `true` | Estimate the coherent reference frame. |
| `--single_circuit` | no | `false` | Trace a single circuit. |
| `--extract_lines` | no | `false` | Extract line segments. |
| `--selection_only` | no | `false` | Use only selected particles. |
| `--print_results` | no | `false` | Print a results summary. |
| `--circuitAtom1` | no | `0` | Single-circuit atom 1. |
| `--circuitAtom2` | no | `0` | Single-circuit atom 2. |
| `--n` | no | `0,0,1` | Interface plane normal (x, y, z). |
| `--xA` | no | `0,0,0` | Orientation x (Grain A); used when `--estimateF false`. |
| `--yA` | no | `0,0,0` | Orientation y (Grain A); used when `--estimateF false`. |
| `--xB` | no | `0,0,0` | Orientation x (Grain B); used when `--estimateF false`. |
| `--yB` | no | `0,0,0` | Orientation y (Grain B); used when `--estimateF false`. |
| `--EcohA` | no | — | Coherency strain (Grain A, 3x3 rows); used when `--estimateF false`. |
| `--EcohB` | no | — | Coherency strain (Grain B, 3x3 rows); used when `--estimateF false`. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_bonds.parquet` | Line Defects | BondExporter → glb |
| `{output_base}_disconnection_summary.parquet` | Disconnection Modes | — (listing-only) |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins
