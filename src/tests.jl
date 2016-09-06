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
dV(x) = ff - v_0*2*(2pi/LL)*cos(2*(2pi/LL)*x)  # Derivative of the potential.
dpotential = [dV(xAxis[1] - dx) ; Float64[dV(x) for x in xAxis] ;
                dV(xAxis[end] + dx)]
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
# Time steps used in the simulations.
timeSteps = 0.1*evolveTime*[1/2048, 1/4096]
timeStepsSystem = timeSteps  # Time steps for system evolution.
facts("Convergence tests.") do
    temperature_step1 = evolveT(temperature, evolveTime, timeSteps[1],
            potential, dpotential, P0, alpha, beta, energy, xAxis)
    temperature_step2 = evolveT(temperature, evolveTime, timeSteps[2],
            potential, dpotential, P0, alpha, beta, energy, xAxis)
    @fact (norm(temperature_step1 - temperature_step2)
        /mean([norm(temperature_step1), norm(temperature_step2)])
            --> less_than(tol) ) "Temperature evolution is not converging."

    density_step1 = evolveP(P0, evolveTime, timeSteps[1], dpotential,
                        temperature, xAxis)
    density_step2 = evolveP(P0, evolveTime, timeSteps[2], dpotential,
                        temperature, xAxis)
    @fact (norm(density_step1 - density_step2)
        /mean([norm(density_step1), mean(density_step2)])
            --> less_than(tol) ) "Probability evolution is not convering."
    # Evolve the system forward and check that the normed differences of the
    # tempearture and the probability density are changing by less than the
    # tolerance. The system evolution is much less convergent than the seperate
    # functions, so we need to decrease the step size.
    system_time_step1 = evolve_system(P0, temperature, evolveTime,
                            timeStepsSystem[1], potential, dpotential, alpha, beta,
                            energy, xAxis)

    system_time_step2 = evolve_system(P0, temperature, evolveTime,
                            timeStepsSystem[2], potential, dpotential, alpha, beta,
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
