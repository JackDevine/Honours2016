# Simulate the system to make some videos, the results will be saved into a text
# file to be visualized by Mathematica.

include("finiteDifferences.jl")
using PyPlot
close("all")

# Directory for storing the results.
prettyVideoDirectory = (
"C:\\Users\\devja964\\Documents\\Project\\Code" )

v_0 = 0.08 # Amplitude of the oscillations in the potential
ff = 3.0 # Forcing for the potential
LL = 1.0 # Length of one period
T0 = 13.0 # Temperature at the ends (i.e. bath temperature)

alpha = -0.0026
beta = 0.0001

nPeriods = 6 # The number of periods to stretch out for
nPoints = 600 # The number of points on the xAxis
nSteps = 1000
dt = 1e-4 # Size of the time step (keep this much smaller than the grid
             # spacing)

# xAxis = linspace(-LL*nPeriods, LL*nPeriods, nPoints)
xAxis = linspace(0LL, 10LL, nPoints)
dx = (xAxis[end] - xAxis[1])/nPoints # Grid spacing

# V(x) = ff*x^2 + 6ff # Potential
sigma = 0.05
bump = 9.3
# V(x) = 0.3x^3 + 5.0ff*exp(-((x-bump)^2)/(2sigma^2))
V(x) = 50x + 10*sin((2pi/LL)*x)
# V(x) = 1.0
potential = Float64[V(x) for x in xAxis]
potential0, potentialEnd = V(xAxis[1] - dx), V(xAxis[end] + dx)
potentialTup = (potential0, potential, potentialEnd)

# temp_fun(x) = (T0 + 0.1*sin(2*(2pi/LL)*x))
tempFun(x) = T0
temperature = Float64[tempFun(x) for x in xAxis]

sigma  = 0.1
P0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x - 7)^2)/(2sigma^2))
              for x in xAxis]
# P0 = ones(xAxis)
P0 /= discrete_quad(P0, xAxis[1], xAxis[end])
density = P0

T = Array(Float64, nPoints, nSteps)
T[:, 1] = temperature
P = Array(Float64, nPoints, nSteps)
P[:, 1] = density
energy = energyFun(potential, density, temperature, alpha, xAxis)

for i = 2:nSteps
    density = stepP(density, dt, potentialTup, temperature, xAxis)
    # temperature = stepT(temperature, dt, density, potentialTup, alpha,
    #                     beta, energy, xAxis)
    T[:, i] = temperature
    P[:, i] = density
end

plot(xAxis, temperature, xAxis, density, xAxis, 0.02potential)
writedlm("$(prettyVideoDirectory)\\xAxis.txt", xAxis)
writedlm("$(prettyVideoDirectory)\\potential.txt", potential)
writedlm("$(prettyVideoDirectory)\\density.txt", P)
writedlm("$(prettyVideoDirectory)\\temperature.txt", T)
