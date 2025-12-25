###############################################################################
# MDBNCH - Molecular Dynamics Benchmark
# A simple molecular dynamics simulation benchmark
# Translated from Fortran 77 to Julia
###############################################################################

# Global data structures replacing Fortran COMMON blocks
mutable struct Particles
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

mutable struct Velocities
    vx::Vector{Float64}
    vy::Vector{Float64}
    vz::Vector{Float64}
end

mutable struct Forces
    fx::Vector{Float64}
    fy::Vector{Float64}
    fz::Vector{Float64}
end

mutable struct Parameters
    n::Int
    dt::Float64
    side::Float64
    rcoff::Float64
end

###############################################################################
# Initialize positions on FCC lattice
###############################################################################
function initpo!(particles::Particles, params::Parameters)
    ncell = round(Int, (params.n / 4)^(1.0 / 3.0))
    dcell = params.side / ncell
    
    itel = 0
    for lz in 1:ncell, ly in 1:ncell, lx in 1:ncell
        x0 = (lx - 1) * dcell
        y0 = (ly - 1) * dcell
        z0 = (lz - 1) * dcell
        
        itel += 1
        if itel > params.n
            break
        end
        particles.x[itel] = x0
        particles.y[itel] = y0
        particles.z[itel] = z0
        
        itel += 1
        if itel > params.n
            break
        end
        particles.x[itel] = x0 + 0.5 * dcell
        particles.y[itel] = y0 + 0.5 * dcell
        particles.z[itel] = z0
        
        itel += 1
        if itel > params.n
            break
        end
        particles.x[itel] = x0 + 0.5 * dcell
        particles.y[itel] = y0
        particles.z[itel] = z0 + 0.5 * dcell
        
        itel += 1
        if itel > params.n
            break
        end
        particles.x[itel] = x0
        particles.y[itel] = y0 + 0.5 * dcell
        particles.z[itel] = z0 + 0.5 * dcell
    end
    
    return nothing
end

###############################################################################
# Initialize velocities
###############################################################################
function initvl!(velocities::Velocities, params::Parameters)
    sumvx = 0.0
    sumvy = 0.0
    sumvz = 0.0
    
    for i in 1:params.n
        velocities.vx[i] = rand() - 0.5
        velocities.vy[i] = rand() - 0.5
        velocities.vz[i] = rand() - 0.5
        sumvx += velocities.vx[i]
        sumvy += velocities.vy[i]
        sumvz += velocities.vz[i]
    end
    
    sumvx /= params.n
    sumvy /= params.n
    sumvz /= params.n
    
    for i in 1:params.n
        velocities.vx[i] -= sumvx
        velocities.vy[i] -= sumvy
        velocities.vz[i] -= sumvz
    end
    
    return nothing
end

###############################################################################
# Calculate forces using Lennard-Jones potential
###############################################################################
function calcfo!(particles::Particles, forces::Forces, params::Parameters)
    rcoffs = params.rcoff * params.rcoff
    
    # Zero out forces
    fill!(forces.fx, 0.0)
    fill!(forces.fy, 0.0)
    fill!(forces.fz, 0.0)
    
    for i in 1:(params.n - 1)
        for j in (i + 1):params.n
            xx = particles.x[i] - particles.x[j]
            yy = particles.y[i] - particles.y[j]
            zz = particles.z[i] - particles.z[j]
            
            # Apply minimum image convention
            xx = xx - params.side * round(xx / params.side)
            yy = yy - params.side * round(yy / params.side)
            zz = zz - params.side * round(zz / params.side)
            
            rd = xx * xx + yy * yy + zz * zz
            
            if rd <= rcoffs
                r2i = 1.0 / rd
                r6i = r2i * r2i * r2i
                ff = 48.0 * r2i * r6i * (r6i - 0.5)
                
                forces.fx[i] += ff * xx
                forces.fy[i] += ff * yy
                forces.fz[i] += ff * zz
                
                forces.fx[j] -= ff * xx
                forces.fy[j] -= ff * yy
                forces.fz[j] -= ff * zz
            end
        end
    end
    
    return nothing
end

###############################################################################
# First half of velocity Verlet integration
###############################################################################
function movea!(particles::Particles, velocities::Velocities, 
                forces::Forces, params::Parameters)
    for i in 1:params.n
        velocities.vx[i] += 0.5 * params.dt * forces.fx[i]
        velocities.vy[i] += 0.5 * params.dt * forces.fy[i]
        velocities.vz[i] += 0.5 * params.dt * forces.fz[i]
        
        particles.x[i] += params.dt * velocities.vx[i]
        particles.y[i] += params.dt * velocities.vy[i]
        particles.z[i] += params.dt * velocities.vz[i]
        
        # Apply periodic boundary conditions
        particles.x[i] -= params.side * round(particles.x[i] / params.side)
        particles.y[i] -= params.side * round(particles.y[i] / params.side)
        particles.z[i] -= params.side * round(particles.z[i] / params.side)
    end
    
    return nothing
end

###############################################################################
# Second half of velocity Verlet integration
###############################################################################
function moveb!(velocities::Velocities, forces::Forces, params::Parameters)
    for i in 1:params.n
        velocities.vx[i] += 0.5 * params.dt * forces.fx[i]
        velocities.vy[i] += 0.5 * params.dt * forces.fy[i]
        velocities.vz[i] += 0.5 * params.dt * forces.fz[i]
    end
    
    return nothing
end

###############################################################################
# Calculate kinetic energy
###############################################################################
function vkcalc(velocities::Velocities, params::Parameters)
    sumv = 0.0
    for i in 1:params.n
        sumv += velocities.vx[i]^2 + velocities.vy[i]^2 + velocities.vz[i]^2
    end
    
    return 0.5 * sumv
end

###############################################################################
# Main program
###############################################################################
function mdbnch()
    # Initialize parameters
    nmax = 864
    params = Parameters(nmax, 0.001, 6.8, 2.5)
    
    # Initialize data structures
    particles = Particles(
        zeros(Float64, nmax),
        zeros(Float64, nmax),
        zeros(Float64, nmax)
    )
    
    velocities = Velocities(
        zeros(Float64, nmax),
        zeros(Float64, nmax),
        zeros(Float64, nmax)
    )
    
    forces = Forces(
        zeros(Float64, nmax),
        zeros(Float64, nmax),
        zeros(Float64, nmax)
    )
    
    println("MOLECULAR DYNAMICS BENCHMARK")
    println("Number of particles: ", params.n)
    println("Timestep: ", params.dt)
    println("Box side: ", params.side)
    println("Cutoff radius: ", params.rcoff)
    println()
    
    # Initialize positions and velocities
    initpo!(particles, params)
    initvl!(velocities, params)
    
    # Time the simulation
    t1 = time()
    
    # Main simulation loop
    for istep in 1:100
        calcfo!(particles, forces, params)
        movea!(particles, velocities, forces, params)
        moveb!(velocities, forces, params)
    end
    
    t2 = time()
    
    # Calculate final energy and temperature
    ekin = vkcalc(velocities, params)
    temp = ekin * 2.0 / (3.0 * params.n)
    
    println("Simulation completed")
    println("Final kinetic energy: ", ekin)
    println("Final temperature: ", temp)
    println("Elapsed time: ", t2 - t1, " seconds")
    
    return nothing
end

# Run the benchmark if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    mdbnch()
end
