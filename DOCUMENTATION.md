# MDBNCH - Molecular Dynamics Benchmark

This repository contains a molecular dynamics simulation benchmark (MDBNCH) originally written in Fortran 77, along with its translation to Julia.

## Overview

MDBNCH is a simple molecular dynamics simulation that models the behavior of particles interacting through a Lennard-Jones potential. It uses the velocity Verlet integration method to advance the system in time.

## Files

- `mdbnch.f` - Original Fortran 77 implementation
- `mdbnch.jl` - Julia translation maintaining equivalent functionality

## Running the Benchmark

### Julia Version

To run the Julia version:

```bash
julia mdbnch.jl
```

The program will:
1. Initialize 864 particles on a face-centered cubic (FCC) lattice
2. Assign random initial velocities
3. Run 100 timesteps of molecular dynamics simulation
4. Report the final kinetic energy, temperature, and elapsed time

### Fortran Version

To compile and run the Fortran version (requires a Fortran compiler like gfortran):

```bash
gfortran -o mdbnch mdbnch.f
./mdbnch
```

Note: The SECOND() function in Fortran may need system-specific implementation.

## Translation Details

### Key Differences Between Fortran and Julia Versions

1. **COMMON Blocks → Mutable Structs**
   - Fortran COMMON blocks have been replaced with Julia mutable structs
   - Four main data structures: `Particles`, `Velocities`, `Forces`, and `Parameters`
   - Passed explicitly to functions instead of using global state

2. **Function Signatures**
   - Fortran subroutines → Julia functions with `!` suffix (indicating mutation)
   - Explicit parameter passing instead of implicit COMMON block access

3. **Timing Function**
   - Fortran `SECOND()` → Julia `time()` function

4. **Random Number Generation**
   - Fortran custom `RANF()` → Julia built-in `rand()`
   - Note: Different random seeds will produce different results

5. **Loop Constructs**
   - Fortran `DO` loops → Julia `for` loops
   - Julia uses 1-based indexing like Fortran
   - Julia's range syntax is more concise

6. **Array Initialization**
   - Fortran `PARAMETER (NMAX=864)` → Julia dynamic arrays
   - Julia uses `zeros(Float64, n)` for initialization

7. **Mathematical Functions**
   - Fortran `DNINT()` → Julia `round()`
   - Fortran `DBLE()` → Julia automatic type conversion

## Algorithm

The simulation uses:

- **Lennard-Jones Potential**: Models inter-particle forces
- **Velocity Verlet Integration**: Time integration scheme for equations of motion
- **Periodic Boundary Conditions**: Particles wrap around the simulation box
- **Minimum Image Convention**: Only nearest periodic images interact
- **Cutoff Radius**: Forces are truncated beyond 2.5 units for efficiency

## Parameters

- Number of particles: 864
- Timestep: 0.001
- Box side length: 6.8
- Cutoff radius: 2.5
- Number of timesteps: 100

## Expected Output

Both versions should produce similar (but not identical due to random initialization) output:

```
MOLECULAR DYNAMICS BENCHMARK
Number of particles: 864
Timestep: 0.001
Box side: 6.8
Cutoff radius: 2.5

Simulation completed
Final kinetic energy: ~1069
Final temperature: ~0.825
Elapsed time: <varies by system>
```

## Performance

The Julia version typically performs comparably to optimized Fortran code, especially when Julia's JIT compiler has warmed up. The first run may be slower due to compilation overhead.

## Requirements

- Julia 1.0 or higher for the Julia version
- Fortran 77 compatible compiler (e.g., gfortran) for the Fortran version

## License

This is a translation of a classic benchmark code for educational and performance comparison purposes.
