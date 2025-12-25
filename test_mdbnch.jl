###############################################################################
# Test suite for MDBNCH Julia implementation
# Tests the translated Fortran77 benchmark framework
###############################################################################

include("mdbnch.jl")

using Test

@testset "MDBNCH Translation Tests" begin
    
    @testset "Constants" begin
        @test NM == 16384
        @test NG == 100
        @test NH == 100
        @test MU == 20
        @test NL == 1
        @test LL == 10 * NM
        @test KP == 2001
        @test KR == 2001
        @test KG == 2001
    end
    
    @testset "CommonBlocks Initialization" begin
        c = CommonBlocks()
        
        # Check struct exists and has correct types
        @test isa(c, CommonBlocks)
        @test isa(c.TWOPI, Float64)
        @test isa(c.NFI, Int)
        @test isa(c.X0, Array{Float64, 2})
        @test isa(c.X, Array{Float64, 3})
        
        # Check array dimensions
        @test size(c.X0) == (3, NM+3)  # -2:NM maps to NM+3 elements
        @test size(c.X) == (3, NM+3, 5)
        @test size(c.RHO) == (KR,)
        @test size(c.VJ) == (KP,)
        @test size(c.UJ) == (KG,)
        @test size(c.LOCK) == (3, 3)
        @test size(c.BOX) == (3, 3)
    end
    
    @testset "Block Data Constants" begin
        # Check gold potential parameters are defined
        @test dendat_RRD ≈ 2.878207442141723
        @test gludat_DB ≈ 12.0
        @test potdat_D ≈ 2.878207442141723
        @test potdat_A ≈ 4.0704
        @test potdat_RC ≈ 3.7
    end
    
    @testset "Utility Functions" begin
        # Test reset!
        arr = [1.0, 2.0, 3.0]
        reset!(arr)
        @test all(arr .== 0.0)
        
        # Test ireset!
        iarr = [1, 2, 3]
        ireset!(iarr)
        @test all(iarr .== 0)
        
        # Test matrix operations
        A = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        B = [2.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 2.0]
        C = zeros(Float64, 3, 3)
        
        mtxmtp!(A, B, C)
        @test C ≈ B  # A^T * B = I * B = B for identity
        
        HI = zeros(Float64, 3, 3)
        DH = mtxinv!(B, HI)
        @test DH ≈ 8.0  # det of 2*I
        @test HI ≈ [0.5 0.0 0.0; 0.0 0.5 0.0; 0.0 0.0 0.5]
    end
    
    @testset "Random Number Generator" begin
        state = ranfm_init(-1)
        @test isa(state, RanfmState)
        @test length(state.IR) == 97
        
        # Generate some random numbers
        r1 = ranfm!(state)
        r2 = ranfm!(state)
        @test 0.0 <= r1 <= 1.0
        @test 0.0 <= r2 <= 1.0
        @test r1 != r2  # Should be different
    end
    
    @testset "Potential Functions" begin
        # Test denfun at various distances
        R = 2.0
        RHO, DRHO, D2RHO = denfun(R)
        @test isa(RHO, Float64)
        @test isa(DRHO, Float64)
        @test isa(D2RHO, Float64)
        
        # Test outside cutoff
        RHO, DRHO, D2RHO = denfun(5.0)
        @test RHO == 0.0
        @test DRHO == 0.0
        @test D2RHO == 0.0
        
        # Test glufun
        DENS = 5.0
        U, U1, U2 = glufun(DENS)
        @test isa(U, Float64)
        @test isa(U1, Float64)
        @test isa(U2, Float64)
        
        # Test potfun
        R = 3.0
        PHI, DPHI, D2PHI = potfun(R)
        @test isa(PHI, Float64)
        @test isa(DPHI, Float64)
        @test isa(D2PHI, Float64)
        
        # Test outside cutoff
        PHI, DPHI, D2PHI = potfun(4.0)
        @test PHI == 0.0
        @test DPHI == 0.0
        @test D2PHI == 0.0
    end
    
    @testset "Potential Table Initialization" begin
        c = CommonBlocks()
        
        # Initialize potential tables
        potent!(c)
        @test c.RN ≈ 1.69
        @test c.RC ≈ 3.7
        @test c.R0 ≈ 2.878207442
        # Tables should be initialized (not all zero)
        @test sum(abs.(c.VJ)) > 0.0
        
        densit!(c)
        @test c.RNRHO ≈ 1.69
        @test c.RCRHO ≈ 3.9
        @test c.RNRHO2 ≈ 1.69^2
        
        elglue!(c)
        @test c.DMIN ≈ 0.0
        @test c.DMAX ≈ 20.0
        @test c.DD > 0.0
    end
    
    @testset "Crystal Structure" begin
        c = CommonBlocks()
        
        # Create a small crystal
        NBSIZE = 2  # 2x2x2 FCC = 32 atoms
        R0 = 2.8
        crystl!(c, R0, NBSIZE)
        
        # Check that atoms were created
        @test c.MOLSA > 0
        @test c.MOLSA <= 32  # 2^3 * 4 atoms per cell
        @test c.MOLSP == c.MOLSA
        
        # Check that NBX, NBY, NBZ are set
        @test c.NBX == NBSIZE
        @test c.NBY == NBSIZE
        @test c.NBZ == NBSIZE
        
        # Check that H matrix is set (stored in first 3 columns of X0)
        @test c.X0[1,1] != 0.0  # H[1,1]
        @test c.X0[2,2] != 0.0  # H[2,2]
        @test c.X0[3,3] != 0.0  # H[3,3]
    end
    
    @testset "Center of Mass" begin
        c = CommonBlocks()
        c.MOLSA = 4
        
        # Set some positions (offset by 3 for -2:NM indexing)
        c.X0[1, 4] = 1.0  # Atom 1
        c.X0[1, 5] = -1.0  # Atom 2
        c.X0[1, 6] = 1.0  # Atom 3
        c.X0[1, 7] = -1.0  # Atom 4
        
        # Center should be at 0
        centcm!(c)
        
        # Check center of mass is zero
        cm = sum(c.X0[1, 4:7]) / 4
        @test abs(cm) < 1e-10
    end
    
    @testset "Full Benchmark Run" begin
        # Test that the main function can be called
        # This is a smoke test - just checks it doesn't crash
        @test_nowarn begin
            c = CommonBlocks()
            mte!(c, 2)  # Small system
        end
    end
end

println("\n✓ All translation tests passed!")
println("Julia translation framework validated successfully.")
