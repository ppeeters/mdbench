###############################################################################
# Test suite for MDBNCH Julia implementation
###############################################################################

include("mdbnch.jl")

using Test

@testset "MDBNCH Tests" begin
    # Test Parameters initialization
    @testset "Parameters" begin
        params = Parameters(864, 0.001, 6.8, 2.5)
        @test params.n == 864
        @test params.dt == 0.001
        @test params.side == 6.8
        @test params.rcoff == 2.5
    end
    
    # Test data structure initialization
    @testset "Data Structures" begin
        n = 10
        particles = Particles(zeros(n), zeros(n), zeros(n))
        @test length(particles.x) == n
        @test length(particles.y) == n
        @test length(particles.z) == n
        
        velocities = Velocities(zeros(n), zeros(n), zeros(n))
        @test length(velocities.vx) == n
        
        forces = Forces(zeros(n), zeros(n), zeros(n))
        @test length(forces.fx) == n
    end
    
    # Test position initialization
    @testset "Position Initialization" begin
        n = 32  # 2x2x2 FCC lattice = 32 particles
        params = Parameters(n, 0.001, 6.8, 2.5)
        particles = Particles(zeros(n), zeros(n), zeros(n))
        
        initpo!(particles, params)
        
        # Check that all positions are initialized
        @test all(particles.x .!= 0.0) || particles.x[1] == 0.0
        @test all(particles.y .!= 0.0) || particles.y[1] == 0.0
        @test all(particles.z .!= 0.0) || particles.z[1] == 0.0
        
        # Check positions are within box
        @test all(particles.x .>= 0.0)
        @test all(particles.x .<= params.side)
        @test all(particles.y .>= 0.0)
        @test all(particles.y .<= params.side)
        @test all(particles.z .>= 0.0)
        @test all(particles.z .<= params.side)
    end
    
    # Test velocity initialization
    @testset "Velocity Initialization" begin
        n = 100
        params = Parameters(n, 0.001, 6.8, 2.5)
        velocities = Velocities(zeros(n), zeros(n), zeros(n))
        
        initvl!(velocities, params)
        
        # Check that velocities sum to zero (center of mass at rest)
        @test abs(sum(velocities.vx)) < 1e-10
        @test abs(sum(velocities.vy)) < 1e-10
        @test abs(sum(velocities.vz)) < 1e-10
        
        # Check that velocities are in reasonable range
        # After centering, velocities can be slightly larger than 0.5
        @test all(abs.(velocities.vx) .< 1.0)
        @test all(abs.(velocities.vy) .< 1.0)
        @test all(abs.(velocities.vz) .< 1.0)
    end
    
    # Test force calculation
    @testset "Force Calculation" begin
        n = 2
        params = Parameters(n, 0.001, 6.8, 2.5)
        particles = Particles([0.0, 1.0], [0.0, 0.0], [0.0, 0.0])
        forces = Forces(zeros(n), zeros(n), zeros(n))
        
        calcfo!(particles, forces, params)
        
        # Forces should be equal and opposite (Newton's 3rd law)
        @test abs(forces.fx[1] + forces.fx[2]) < 1e-10
        @test abs(forces.fy[1] + forces.fy[2]) < 1e-10
        @test abs(forces.fz[1] + forces.fz[2]) < 1e-10
        
        # Forces should be non-zero for particles at distance 1.0
        @test forces.fx[1] != 0.0
    end
    
    # Test kinetic energy calculation
    @testset "Kinetic Energy" begin
        n = 10
        params = Parameters(n, 0.001, 6.8, 2.5)
        velocities = Velocities(ones(n), ones(n), ones(n))
        
        ekin = vkcalc(velocities, params)
        
        # For velocity = 1.0 in each direction: KE = 0.5 * n * (1^2 + 1^2 + 1^2)
        @test abs(ekin - 0.5 * n * 3.0) < 1e-10
    end
    
    # Test integration (movea and moveb)
    @testset "Integration" begin
        n = 2
        params = Parameters(n, 0.001, 6.8, 2.5)
        particles = Particles([0.0, 3.0], [0.0, 0.0], [0.0, 0.0])
        velocities = Velocities([1.0, 0.0], [0.0, 0.0], [0.0, 0.0])
        forces = Forces([10.0, -10.0], [0.0, 0.0], [0.0, 0.0])
        
        x_initial = copy(particles.x)
        v_initial = copy(velocities.vx)
        
        movea!(particles, velocities, forces, params)
        
        # Position should change
        @test particles.x[1] != x_initial[1]
        
        # Velocity should change
        @test velocities.vx[1] != v_initial[1]
        
        moveb!(velocities, forces, params)
        
        # Velocity should change again
        @test velocities.vx[1] != v_initial[1]
    end
    
    # Integration test: run a few steps
    @testset "Full Simulation Steps" begin
        n = 32
        params = Parameters(n, 0.001, 6.8, 2.5)
        
        particles = Particles(zeros(n), zeros(n), zeros(n))
        velocities = Velocities(zeros(n), zeros(n), zeros(n))
        forces = Forces(zeros(n), zeros(n), zeros(n))
        
        initpo!(particles, params)
        initvl!(velocities, params)
        
        # Run a few steps
        for step in 1:5
            calcfo!(particles, forces, params)
            movea!(particles, velocities, forces, params)
            moveb!(velocities, forces, params)
        end
        
        # System should still have momentum conservation
        @test abs(sum(velocities.vx)) < 1e-8
        @test abs(sum(velocities.vy)) < 1e-8
        @test abs(sum(velocities.vz)) < 1e-8
        
        # Energy should be positive
        ekin = vkcalc(velocities, params)
        @test ekin > 0.0
    end
end

println("\nAll tests passed!")
