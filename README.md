# IMPACT: Integrated Multi-topology & Multi-criteria Polymeric Analysis of Crystallization & Cluster Toolkit

## Overview
IMPACT is a coarse-grained (CG) crystallization analysis toolkit for semi-crystalline polymers and polymeric nanocomposites. It supports structural order parameter (SOP) and trans-trans conformation (dtt) crystallinity analysis, optional cluster and polymer-chain conformation analysis, and volume-based crystallinity estimation.
This package is designed as a command-line workflow for crystallization analysis from molecular simulation dump files, with optional follow mode for live analysis of a growing trajectory.

---
## Analysis Modes

**SOP-based crystallinity**

  Uses structural order parameters to detect crystalline environments, and can be combined with optional cluster analysis and polymer conformation analysis.

**DTT-based crystallinity**

  Uses trans-trans conformation criteria for crystallinity identification, and can be combined with optional cluster analysis.

**Volume-based crystallinity**

  Estimates crystallinity from probe-volume sampling and can optionally write snapshot-style output.

**Follow mode**

  Supports live analysis of a growing dump file through `-follow`, with configurable polling interval and parser mode.

---
## Build

```bash
make
```

The build produces the `IMPACT` binary by default.

---
## Usage

```bash
IMPACT -in Input_file -out Output_folder Parameter
```

**Main options**

- `-sop R crystal_judge -c sop_R cluster_judge -cf`  
  Example: `-sop 1.44 0.8 -c 1.05 0.95 -cf`

- `-dtt crystal_judge crystal_size_judge -c dtt_R cluster_judge -cf`  
  Example: `-dtt 0.95 14 -c 2.0 0.625 -cf`

- `-ig typenumber`

- `-v Lx_probe (d)`

- `-ti time_interval` (default: 1)

- `-follow` (follow a growing dump file for realtime analysis)

- `-poll milliseconds` (follow polling interval; default: 1000)

- `-slowio` (line-based parser)

- `-fastio` (token-based parser)

**Notes**

- The atom input file must follow the format sequence `id mol type x y z ix iy iz`.
- The simulation box must be triclinic.
- `-cf` (conformation analysis) is only applicable to linear polymers.

**Examples**

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

---
## File Descriptions

- **`main.cpp`**  
  Program entry point and command-line workflow driver.

- **`input.cpp`**  
  Reads and prepares trajectory input data for the analysis routines.

- **`calculate.cpp`**  
  Core crystallinity-analysis routines.

- **`cluster.cpp`**  
  Optional cluster-analysis routines used together with crystallinity classification.

- **`conformation.cpp`**  
  Optional polymer-chain conformation analysis.

- **`volume_cry.cpp`**  
  Volume-based crystallinity estimation routines.

- **`crystal.h`**  
  Shared declarations and analysis data structures.

- **`tests/`**  
  Test or example inputs for validating the workflow.

- **`Makefile`**  
  Build instructions for compiling the `IMPACT` executable.

- **`LICENSE`**  
  The license text for the repository.

---
