include("finiteDifferences.jl")
using PyPlot
close("all")

v_0 = 0.08 # Amplitude of the oscillations in the potential
ff = 3.0 # Forcing for the potential
LL = 1.0 # Length of one period
T0 = 7.0 # Temperature at the ends (i.e. bath temperature)

alpha = 0.0004
beta = 0.1

nPeriods = 6 # The number of periods to stretch out for
nPoints = 600 # The number of points on the xAxis
dt = 1e-5 # Size of the time step (keep this much smaller than the grid
             # spacing)

# xAxis = linspace(-LL*nPeriods, LL*nPeriods, nPoints)
xAxis = linspace(5LL, 10LL, nPoints)
dx = (xAxis[end] - xAxis[1])/nPoints # Grid spacing

# V(x, time) = ff*x^2 + 6ff # Potential
sigma = 0.1
bump = 9
V(x) = 0.3x^3 + 3.0ff*exp(-((x-bump)^2)/(2sigma^2))
potential = Float64[V(x) for x in xAxis]
potential0, potentialEnd = V(xAxis[1] - dx), V(xAxis[end] + dx)
potentialTup = (potential0, potential, potentialEnd)

# temperatureFun(x) = (T0 + 0.1*sin(2*(2pi/LL)*x))
temperatureFun(x) = T0
temperature = Float64[temperatureFun(x) for x in xAxis]

sigma  = 0.05
P0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x-9.3)^2)/(2sigma^2))
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
#                 alpha, beta, energy, xAxis)
# density2 = evolveP(P0, evolve_time, dt, potentialTup, ones(xAxis)*T0, xAxis)
# end
#
# plot(xAxis, 0.02potential, xAxis, density, xAxis, density2, xAxis,
#             temp, xAxis, ones(xAxis)*T0)
# legend(["Potential", "Coupled probability density",
#         "Uncoupled probability density", "Temperature", "Initial temperature"])

function hopping_time(potentialTup::Tuple{Number, AbstractArray, Number},
                bump::Number, density::AbstractArray,
                temperature::AbstractArray, alpha::Number, beta::Number,
                xAxis::AbstractArray; dt=1e-4, tol=0.1)
    #=
    Given an initial state of the system and a potential with a bump in it,
    calculate how long it takes for the probability distribution to get over
    the bump.
    Parameters:
    potential:     An array containg the potential.
    bump:          The location of the bump in the potetnial.
    density:       The initial density of the system.
    temperature:   The initial temperature of the system.
    alpha:         Dimensionless parameter that controls the coupling of the
                   system.
    beta:          Dimensionless parameter that determines how fast the
                   temperature gradients diffuse.
    xAxis:         The axis that we are working on (dimensionless).
    dt:            Size of one time step for simulation.
    tol:           If the fraction of the probability density that is to the
                   left of the bump is less than tol, then we will consider
                   the particle to have crossed the barrier.
    =#
    # Calculate the location of the center of the initial distribution.
    bump_ind = indmin(abs(xAxis - bump))
    iters = 0
    # while discrete_quad(density[bump_ind:end], bump, xAxis[end]) > tol
    while discrete_quad(density.*xAxis, xAxis[1], xAxis[end]) > bump
        density = stepP(density, dt, potentialTup, temperature, xAxis)
        temperature = stepT(temperature, dt, density, potentialTup, alpha,
                            beta, energy, xAxis)
        iters += 1
        # println(discrete_quad(density.*xAxis, xAxis[1], xAxis[end]))
        if iters > 2000
            println("alpha = $alpha , beta = $beta DNF")
            return dt*iters
        end
    end
    ret = dt*iters
end
#
# hopping_time(potentialTup, bump, P0, temperature, alpha, beta, heat_capacity,
#                     xAxis; dt=1e-4, tol=0.7)
nPoints = 2
alphaMin = 0.0002
alphaMax = 0.0003
betaMin = 0.1
betaMax = 40.0
betaVec = linspace(betaMin, betaMax, nPoints)
alphaVec = linspace(alphaMin, alphaMax, nPoints)
@time a = Float64[hopping_time(potentialTup, bump, P0, temperature, alpha, beta,
                               xAxis)
                    for alpha in linspace(alphaMin, alphaMax, nPoints),
                    beta in linspace(betaMin, betaMax, nPoints)]
# @time a = map((alpha, beta) -> hopping_time(potentialTup, bump, P0, temperature, alpha, beta,
#                     heat_capacity, xAxis),
#                     linspace(0.0, 0.8, nPoints),
#                     linspace(1, 4, nPoints))

matshow(a)
xlabel(L"\beta")
ylabel(L"\alpha")
figure()
surf(betaVec, alphaVec, a)
xlabel(L"\beta")
ylabel(L"\alpha")
