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
T0 = 7.0 # Temperature at the ends (i.e. bath temperature)

alpha = 0.00016
beta = 0.000000001

nPeriods = 6 # The number of periods to stretch out for
nPoints = 600 # The number of points on the xAxis
nSteps = 2000
dt = 1e-5 # Size of the time step (keep this much smaller than the grid
             # spacing)

# xAxis = linspace(-LL*nPeriods, LL*nPeriods, nPoints)
xAxis = linspace(5LL, 10LL, nPoints)
dx = (xAxis[end] - xAxis[1])/nPoints # Grid spacing

# V(x) = ff*x^2 + 6ff # Potential
sigma = 0.05
bump = 9
V(x) = 0.3x^3 + 7.0ff*exp(-((x-bump)^2)/(2sigma^2))
# V(x) = 1.0
disV = Float64[V(x) for x in xAxis]
disV0, disVEnd = V(xAxis[1] - dx), V(xAxis[end] + dx)
disVTup = (disV0, disV, disVEnd)

# temp_fun(x) = (T0 + 0.1*sin(2*(2pi/LL)*x))
temp_fun(x) = T0
disTemp = Float64[temp_fun(x) for x in xAxis]
disTemp0, disTempEnd = temp_fun(xAxis[1] - dx), temp_fun(xAxis[end] + dx)
disTempTup = (disTemp0, disTemp, disTempEnd)

sigma  = 0.1
P0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x-9.3)^2)/(2sigma^2))
              for x in xAxis]
# P0 = ones(xAxis)
P0 /= discrete_quad(P0, xAxis[1], xAxis[end])
density = P0

T = Array(Float64, nPoints, nSteps)
T[:, 1] = disTemp
P = Array(Float64, nPoints, nSteps)
P[:, 1] = density

for i = 2:nSteps
    density = stepP(density, dt, disVTup, disTemp, xAxis)
    disTemp = stepT(disTemp, dt, density, disVTup, alpha,
                        beta, energy, xAxis)
    T[:, i] = disTemp
    P[:, i] = density
end

# plot(xAxis, disTemp, xAxis, density, xAxis, 0.02disV)
writedlm("$(prettyVideoDirectory)\\xAxis.txt", xAxis)
writedlm("$(prettyVideoDirectory)\\potential.txt", disV)
writedlm("$(prettyVideoDirectory)\\density.txt", P)
writedlm("$(prettyVideoDirectory)\\temperature.txt", T)
