module Hopping
export hopping_time, potentialTup, bump, P0, temperature, alpha, betad, xAxis
include("finiteDifferences.jl")
using PyPlot
close("all")

v_0 = 0.08 # Amplitude of the oscillations in the potential
ff = 3.0 # Forcing for the potential
LL = 1.0 # Length of one period
T0 = 15.0 # Temperature at the ends (i.e. bath temperature)

alpha = 0.0004
betad = 0.1

nPeriods = 6 # The number of periods to stretch out for
nPoints = 1000 # The number of points on the xAxis
dt = 1e-6 # Size of the time step (keep this much smaller than the grid
             # spacing)

# xAxis = linspace(-LL*nPeriods, LL*nPeriods, nPoints)
xAxis = linspace(-2LL, 2LL, nPoints)
dx = (xAxis[end] - xAxis[1])/nPoints # Grid spacing

# V(x, time) = ff*x^2 + 6ff # Potential
sigma = 0.1
bump = 0.0
# V(x) = 0.25x^3 + 3.2ff*exp(-((x-bump)^2)/(2sigma^2))
# V(x) = 28. - 56.35x + 37.1333x^2 - 9.65x^3 + 0.866667x^4
V(x) = 10*(1.9abs(x) - 1.1)^2 + 5x
potential = Float64[V(x) for x in xAxis]
potential0, potentialEnd = V(xAxis[1] - dx), V(xAxis[end] + dx)
potentialTup = (potential0, potential, potentialEnd)

# temperatureFun(x) = (T0 + 0.1*sin(2*(2pi/LL)*x))
temperatureFun(x) = T0
temperature = Float64[temperatureFun(x) for x in xAxis]

sigma  = 0.05
P0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x-0.6)^2)/(2sigma^2))
              for x in xAxis]
# P0 = ones(xAxis)
P0 /= discrete_quad(P0, xAxis[1], xAxis[end])
density = P0
# Calculate the initial energy of the system.
energy = energyFun(potential, P0, temperature, alpha, xAxis)
boltzmann_density = exp(-potential/T0)
boltzmann_density = boltzmann_density/discrete_quad(boltzmann_density,
                xAxis[1], xAxis[end])
# evolve_time = 1000dt
# @time begin
# density, temp = evolve_system(density, temperature, evolve_time, dt, potentialTup,
#                 alpha, betad, energy, xAxis)
# density2 = evolveP(P0, evolve_time, dt, potentialTup, ones(xAxis)*T0, xAxis)
# end
#
# plot(xAxis, 0.02potential, xAxis, density, xAxis, density2, xAxis,
#             temp, xAxis, ones(xAxis)*T0)
# legend(["Potential", "Coupled probability density",
#         "Uncoupled probability density", "Temperature", "Initial temperature"])
"""
    hopping_time(potentialTup::Tuple{Number, AbstractArray, Number},
                bump::Number, density::AbstractArray,
                temperature::AbstractArray, alpha::Number, betad::Number,
                xAxis::AbstractArray; dt=1e-4, tol=0.1)
Given an initial state of the system and a potential with a bump in it,
calculate how long it takes for the probability distribution to get over
the bump.
# Arguments:
* `potentialTup::Tuple{Number, AbstractArray, Number}`: A tuple containing the
potential.
* `bump::Number`: The location of the bump in the potetnial.
* `density::AbstractArray`: The initial density of the system.
* `temperature::AbstractArray`: The initial temperature of the system.
* `alpha::Number`: Dimensionless parameter that controls the coupling of the
system.
* `betad::Number`: Dimensionless parameter that determines how fast the
temperature gradients diffuse.
* `xAxis::AbstractArray`: The axis that we are working on.
* `dt=1e-4`: Size of one time step for simulation.
"""
function hopping_time(potentialTup::Tuple{Number, AbstractArray, Number},
                bump::Number, density::AbstractArray,
                temperature::AbstractArray, alpha::Number, betad::Number,
                xAxis::AbstractArray; dt=1e-4)
    bump_ind = indmin(abs(xAxis - bump))
    iters = 0
    # If the mean of the probability distribution crosses the bump, then the
    # system has crossed the bump.
    while discrete_quad(density.*xAxis, xAxis[1], xAxis[end]) > bump
        density = stepP(density, dt, potentialTup, temperature, xAxis)
        temperature = stepT(temperature, dt, density, potentialTup, alpha,
                            betad, energy, xAxis)
        iters += 1
        if iters > 10000
            println("alpha = $alpha , betad = $betad DNF")
            return dt*iters
        end
    end
    ret = dt*iters
end

hopping_time(alpha, betad) = hopping_time(potentialTup, bump, P0,
                temperature, alpha, betad, xAxis)

end # Module.
