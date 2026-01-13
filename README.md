# IMPACT (Integrated Multi-topology & Multi-criteria Polymeric Analysis of Crystallization & Cluster Toolkit)

IMPACT is a coarse-grained (CG) crystallization analysis toolkit for semi-crystalline polymers and polymeric nanocomposites. It supports structural order parameter (SOP) and DTT-based crystallinity analysis, optional cluster and conformation analysis, and volume-based crystallinity estimation. IMPACT can follow a growing dump file for real-time analysis and provides both fast (token-based) and defensive (line-based) input parsers.

## Features

- SOP-based crystallinity analysis with optional cluster and conformation analysis.
- DTT-based crystallinity analysis with optional cluster analysis.
- Volume-based crystallinity analysis with optional snapshot output.
- Follow mode for live analysis of a growing dump file.
- Fast token-based parser and slower line-based parser.

## Build

```bash
make
```

The build produces the `IMPACT` binary by default.

## Usage

```bash
IMPACT -in Input_file -out Output_folder Parameter
```

Parameters:

- `-sop R crystal_judge -c sop_R cluster_judge -cf`
  - Example: `-sop 1.44 0.8 -c 1.05 0.95 -cf`
- `-dtt crystal_judge crystal_size_judge -c dtt_R cluster_judge -cf`
  - Example: `-dtt 0.95 14 -c 2.0 0.625 -cf`
- `-ig typenumber`
- `-v Lx_probe (d)`
- `-ti time_interval` (default: 1)
- `-follow` (follow a growing dump file for realtime analysis)
- `-poll milliseconds` (follow polling interval; default: 1000)
- `-slowio` (line-based parser)
- `-fastio` (token-based parser)

Notes:

- The atom input file must follow this format sequence: `id mol type x y z ix iy iz`.
- The simulation box must be triclinic.
- `-cf` (conformation analysis) is only applicable to linear polymers.

## Examples

SOP analysis with cluster and conformation:

```bash
IMPACT -in input.dump -out output -sop 1.44 0.8 -c 1.05 0.95 -cf
```

DTT analysis with clustering:

```bash
IMPACT -in input.dump -out output -dtt 0.95 14 -c 2.0 0.625
```

Follow mode with fast parser:

```bash
IMPACT -in input.dump -out output -sop 1.44 0.8 -follow -fastio
```
