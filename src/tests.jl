using FactCheck
include("finiteDifferences.jl")

# Initialize a particular starting system.
v_0 = 0.8 # Amplitude of the oscillations in the potential
ff = 3.0 # Forcing for the potential
LL = 1.0 # Length of one period
T0 = 1.0 # Temperature at the ends (i.e. bath temperature)

alpha = 0.00002
beta = 0.15

nPeriods = 4 # The number of periods to stretch out for.
nPoints = 200 # The number of points on the xAxis.

xAxis = linspace(0, LL*nPeriods, nPoints)
dx = (xAxis[end] - xAxis[1])/nPoints # Grid spacing

V(x, time) = ff*x - v_0*sin(2*(2pi/LL)*x) + 6ff # Potential
potential = Float64[V(x, 0) for x in xAxis]
potential0, potentialEnd = V(xAxis[1] - dx, 0), V(xAxis[end] + dx, 0)
potentialTup = (potential0, potential, potentialEnd)

temperatureFun(x) = (T0 + 0.1*sin(2*(2pi/LL)*x))
temperature = Float64[temperatureFun(x) for x in xAxis]

sigma  = 0.1
P0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x-3.0)^2)/(2sigma^2))
              for x in xAxis]
P0 /= discrete_quad(P0, xAxis[1], xAxis[end])

energy = energyFun(potential, P0, temperature, alpha, xAxis)
# Evolve the system using different time steps, for small time steps, halfing
# the time step should not cause the result of the simulation to change by much.
tol = 1e-3  # The amount that the two results are allowed to differ by.
evolveTime = 1.5  # The amount of time that we will evolve the system for.
timeSteps = evolveTime*[1/2048, 1/4096]  # Time steps used in the simulations.
timeStepsSystem = timeSteps  # Time steps for system evolution.
facts("Convergence tests.") do
    temperature_step1 = evolveT(temperature, evolveTime, timeSteps[1],
            potentialTup, P0, alpha, beta, energy, xAxis)
    temperature_step2 = evolveT(temperature, evolveTime, timeSteps[2],
            potentialTup, P0, alpha, beta, energy, xAxis)
    @fact (norm(temperature_step1 - temperature_step2)
        /mean([norm(temperature_step1), norm(temperature_step2)])
            --> less_than(tol) ) "Temperature evolution is not converging."

    density_step1 = evolveP(P0, evolveTime, timeSteps[1], potentialTup,
                        temperature, xAxis)
    density_step2 = evolveP(P0, evolveTime, timeSteps[2], potentialTup,
                        temperature, xAxis)
    @fact (norm(density_step1 - density_step2)
        /mean([norm(density_step1), mean(density_step2)])
            --> less_than(tol) ) "Probability evolution is not convering."
    # Evolve the system forward and check that the normed differences of the
    # tempearture and the probability density are changing by less than the
    # tolerance. The system evolution is much less convergent than the seperate
    # functions, so we need to decrease the step size.
    system_time_step1 = evolve_system(P0, temperature, evolveTime,
                            timeStepsSystem[1], potentialTup, alpha, beta,
                            energy, xAxis)

    system_time_step2 = evolve_system(P0, temperature, evolveTime,
                            timeStepsSystem[2], potentialTup, alpha, beta,
                            energy, xAxis)

    @fact (norm(system_time_step1[2] - system_time_step2[2])
        /mean([norm(system_time_step1[2]), norm(system_time_step2[2])])
            --> less_than(tol) ) """The temperature in the system evolution is
                                    not converging."""
    @fact (norm(system_time_step1[1] - system_time_step2[1])
        /mean([norm(system_time_step1[1]), norm(system_time_step2[1])])
            --> less_than(tol) ) """The probability density in the system
                                    evolution is not converging."""
end

# Define the physics for the hopping_time function.
LL = 1.0 # Length of one period
T0 = 7.0 # Temperature at the ends (i.e. bath temperature)

alpha = 0.0004
beta = 0.1

nPeriods = 6 # The number of periods to stretch out for
nPoints = 600 # The number of points on the xAxis
dt = 1e-4 # Size of the time step (keep this much smaller than the grid
             # spacing)
xAxis = linspace(-2LL, 2LL, nPoints)
dx = (xAxis[end] - xAxis[1])/nPoints # Grid spacing

sigma = 0.1
bump = 0.0
sigma = 0.1
V(x) = x^4 - 3*(exp((-(x - 0.5)^2)/sigma) + exp((-(x + 0.5)^2)/sigma)) + 3x

potential = Float64[V(x) for x in xAxis]
potential0, potentialEnd = V(xAxis[1] - dx), V(xAxis[end] + dx)
potentialTup = (potential0, potential, potentialEnd)

temperatureFun(x) = T0
temperature = Float64[temperatureFun(x) for x in xAxis]

sigma  = 0.05
P0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x-0.5)^2)/(2sigma^2))
              for x in xAxis]
# P0 = ones(xAxis)
P0 /= discrete_quad(P0, xAxis[1], xAxis[end])
density = P0

# Define a method of hopping_time for this particular potential.
hopping_time(alpha, beta) = hopping_time(potentialTup, bump, density, temperature, alpha, beta,
                xAxis; dt=1e-4)
facts("Hopping time") do
    @fact ( hopping_time(potentialTup, 1.0, density, temperature, alpha,
                    beta, xAxis; dt=1e-4)
                    --> 0.0)
                    "This system has already crossed the supposed bump."
end
