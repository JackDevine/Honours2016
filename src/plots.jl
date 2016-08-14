include("finiteDifferences.jl")
using PyPlot
close("all")

v_0 = 0.08 # Amplitude of the oscillations in the potential
ff = 3.0 # Forcing for the potential
LL = 1.0 # Length of one period
T0 = 15.0 # Temperature at the ends (i.e. bath temperature)

alpha = -0.000392
betad = 0.005

nPeriods = 6 # The number of periods to stretch out for
nPoints = 600 # The number of points on the xAxis
dt = 2e-5 # Size of the time step (keep this much smaller than the grid
             # spacing)

# xAxis = linspace(-LL*nPeriods, LL*nPeriods, nPoints)
xAxis = linspace(-2LL, 2LL, nPoints)
dx = (xAxis[end] - xAxis[1])/nPoints # Grid spacing

# V(x, time) = ff*x^2 + 6ff # Potential
sigma = 0.1
bump = 9
# V(x) = 0.2x^3 + 3.2ff*exp(-((x-bump)^2)/(2sigma^2))
# V(x) = 2 - 18.7333x + 28.7222x^2 - 14.3889x^3 + 2.44444x^4 - 0.0444444x^5
V(x) = 10*(1.9abs(x) - 1.1)^2
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
evolve_time = 5000dt
@time begin
density, temp = evolve_system(density, temperature, evolve_time, dt, potentialTup, alpha, betad, energy, xAxis)
density2 = evolveP(P0, evolve_time, dt, potentialTup, ones(xAxis)*T0, xAxis)
end

plot(xAxis, 0.05potential, xAxis, 5density, xAxis, 5density2, xAxis,
            temp, xAxis, ones(xAxis)*T0)
legend(["Potential", "Coupled probability density",
        "Uncoupled probability density", "Temperature", "Initial temperature"], loc=2)
