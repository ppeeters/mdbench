#!/bin/bash
###############################################################################
# Script to run MDBNCH implementations
###############################################################################

echo "MDBNCH - Molecular Dynamics Benchmark"
echo "======================================"
echo ""

# Check if Julia is installed
if command -v julia &> /dev/null; then
    echo "✓ Julia found: $(julia --version)"
    JULIA_AVAILABLE=1
else
    echo "✗ Julia not found"
    JULIA_AVAILABLE=0
fi

# Check if gfortran is installed
if command -v gfortran &> /dev/null; then
    echo "✓ gfortran found: $(gfortran --version | head -n1)"
    FORTRAN_AVAILABLE=1
else
    echo "✗ gfortran not found"
    FORTRAN_AVAILABLE=0
fi

echo ""

# Run Julia version
if [ $JULIA_AVAILABLE -eq 1 ]; then
    echo "Running Julia version..."
    echo "------------------------"
    julia mdbnch.jl
    echo ""
    
    echo "Running tests..."
    echo "------------------------"
    julia test_mdbnch.jl
    echo ""
    
    echo "Running comparison..."
    echo "------------------------"
    julia compare.jl
    echo ""
else
    echo "Skipping Julia version (Julia not installed)"
    echo ""
fi

# Compile and run Fortran version
if [ $FORTRAN_AVAILABLE -eq 1 ]; then
    echo "Compiling Fortran version..."
    echo "----------------------------"
    
    # Try to compile (may fail due to SECOND function)
    if gfortran -o mdbnch_fortran mdbnch.f 2>&1; then
        echo "✓ Compilation successful"
        echo ""
        echo "Running Fortran version..."
        echo "--------------------------"
        ./mdbnch_fortran
        echo ""
        
        # Clean up
        rm -f mdbnch_fortran
    else
        echo "✗ Compilation failed"
        echo "Note: The SECOND() function may need system-specific implementation"
        echo ""
    fi
else
    echo "Skipping Fortran version (gfortran not installed)"
    echo ""
fi

echo "======================================"
echo "Benchmark complete!"
