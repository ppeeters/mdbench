# MDBNCH - Molecular Dynamics Benchmark

Translation of the classic Fortran 77 MDBNCH molecular dynamics benchmark to Julia.

## Overview

MDBNCH is a molecular dynamics benchmark for gold atoms using a many-body 'glue' (embedded atom) interaction potential. The original Fortran77 code (2141 lines, 30 program units) has been systematically translated to Julia.

## Quick Start

Run the Julia version:
```bash
julia mdbnch.jl
```

Expected output:
```
MDBNCH: A MOLECULAR DYNAMICS BENCHMARK, VERSION OF DECEMBER 17, 1988

MD BENCHMARK FOR 256 PARTICLES, 1000 STEPS.
...
COMPLETE BENCHMARK EXECUTION TIME : ~0.5 CP SECONDS.
```

## Files

- `mdbnch.f` - Original Fortran 77 source code (2141 lines)
- `mdbnch.jl` - Julia translation (875 lines)
- `DOCUMENTATION.md` - Detailed documentation
- `TRANSLATION_SUMMARY.md` - Translation methodology and status

## Translation Approach

The Julia translation captures the full structure of the Fortran benchmark:
- 19 COMMON blocks → unified `CommonBlocks` mutable struct
- Array bounds like `(-2:NM)` handled via index offsetting  
- BLOCK DATA → module-level constants
- DOUBLE PRECISION → Float64
- SECOND() timing → time()
- All potential functions fully translated

## Documentation

See [DOCUMENTATION.md](DOCUMENTATION.md) and [TRANSLATION_SUMMARY.md](TRANSLATION_SUMMARY.md) for complete details about the translation approach, algorithm, and implementation status.
