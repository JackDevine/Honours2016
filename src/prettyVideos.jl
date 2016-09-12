# Simulate the system to make some videos, the results will be saved into a text
# file to be visualized by Mathematica.

include("finiteDifferences.jl")
using FiniteDifferences
using PyPlot
using ForwardDiff
close("all")

# Directory for storing the results.
prettyVideoDirectory = (
"/Users/jackdevine/Desktop/MolecularMotors/Code" )
params = """
slope = 0.0  # Slope of the potential.
LL = 1.0 # Length of one period
fName = "uncoupledFlatPotential.txt"  # The name of the output file.
T0 = 10.0  # Temperature.
alpha = -0.0026
beta = 2.0
v_0 = 0.0 # Amplitude of the oscillations in the potential

nPeriods = 6 # The number of periods to stretch out for
nPoints = 1000 # The number of points on the xAxis
nSteps = 5000
dt = 5e-6 # Size of the time step (keep this much smaller than the grid
             # spacing)
meanInit = 0.0  # The initial center of the probability distribution.
xi, xe = -3, 3
V(x) = slope*x + 1 + v_0*sin((2pi/LL)*x)
skip = 100  # Frames to skip.
"""

outfile = open("$(prettyVideoDirectory)/$fName", "w")
write(outfile, params)
close(outfile)
f = open("$(prettyVideoDirectory)/$fName")
params = readlines(f)
close(f)
for line in params
    ex = parse(line)
    eval(ex)
end
# xAxis = linspace(-LL*nPeriods, LL*nPeriods, nPoints)
xAxis = linspace(xi*LL, xe*LL, nPoints)
dx = (xAxis[end] - xAxis[1])/nPoints # Grid spacing

# V(x) = ff*x^2 + 6ff # Potential
# V(x) = 0.3x^3 + 5.0ff*exp(-((x-bump)^2)/(2sigma^2))
# V(x) = slope*x + 1 + v_0*sin((2pi/LL)*x)
# V(x) = 1.0
potential = Float64[V(x) for x in xAxis]
dV(x) = ForwardDiff.derivative(V, x)
dpotential = [dV(xAxis[1] - dx) ; Float64[dV(x) for x in xAxis] ;
                                                dV(xAxis[end] + dx)]

# temp_fun(x) = (T0 + 0.1*sin(2*(2pi/LL)*x))
tempFun(x) = T0
temperature = Float64[tempFun(x) for x in xAxis]

sigma  = 0.1
P0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x-meanInit)^2)/(2sigma^2))
              for x in xAxis]
# P0 = ones(xAxis)
P0 /= discrete_quad(P0, xAxis[1], xAxis[end])
density = P0

system = System()
system.density = density
system.potential = potential
system.dpotential = dpotential
system.temperature = temperature
system.xAxis = xAxis
system.energy = energyFun(system, alpha)

T = Array(Float64, nPoints, nSteps)
T[:, 1] = temperature
P = Array(Float64, nPoints, nSteps)
P[:, 1] = density
heat = Array(Float64, nPoints - 1, nSteps)
heat[:, 1] = ( (density[2:end].*discrete_derivative(potential, xAxis)
            + discrete_derivative(density.*potential, xAxis))
            .*discrete_derivative(potential, xAxis) )

for i = 2:nSteps
    system.density = stepP(system, dt)
    system.temperature = stepT(system, alpha, beta, dt)
    heat[:, i] = ( (density[2:end].*discrete_derivative(system.potential,
                                                                system.xAxis)
                 + discrete_derivative(system.density.*system.potential,
                                                                system.xAxis))
                 .*discrete_derivative(system.potential, system.xAxis) )
    T[:, i] = system.temperature
    P[:, i] = system.density
end

plot(xAxis[2:end], 0.000001heat[:, end], system.xAxis, system.density,
        system.xAxis, 0.02system.potential, system.xAxis, system.temperature)
writedlm("$(prettyVideoDirectory)/xAxis.txt", system.xAxis)
writedlm("$(prettyVideoDirectory)/potential.txt", system.potential)
writedlm("$(prettyVideoDirectory)/density.txt", P[:, 1:skip:end])
writedlm("$(prettyVideoDirectory)/temperature.txt", T[:, 1:skip:end])
writedlm("$(prettyVideoDirectory)/heat.txt", heat[:, 1:skip:end])
