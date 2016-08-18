@everywhere include("finiteDifferences.jl")
using PyPlot
close("all")

@everywhere begin
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
V(x) = 0.3x^3 + 5.0ff*exp(-((x-bump)^2)/(2sigma^2))
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

nPoints = 50
alphaMin = -0.002
alphaMax = -0.0001
betaMin = 0.000001
betaMax = 0.001
betaVec = linspace(betaMin, betaMax, nPoints)
alphaVec = linspace(alphaMin, alphaMax, nPoints)
alphaBeta = [[alpha, beta] for alpha in alphaVec, beta in betaVec]
end # @everywhere

# println("Array comprehension:")
# @time a = Float64[hopping_time(potentialTup, bump, P0, temperature, alpha, beta,
#                                xAxis)
#                     for alpha in alphaVec, beta in betaVec]
# println("map:")
# @time b = map(z -> hopping_time(z[1], z[2]), alphaBeta[:])
# b = reshape(a, nPoints, nPoints)

println("pmap:")
@time c = pmap(z -> hopping_time(z[1], z[2]), alphaBeta[:])
c = reshape(c, nPoints, nPoints)

# println("Parallel loop")
# d = SharedArray(Float64, nPoints, nPoints)
# @time @sync @parallel for i in 1:nPoints
#     for j in 1:nPoints
#         d[i, j] = hopping_time(alphaVec[i], betaVec[j])
#     end
# end


matshow(c)
xlabel(L"\beta")
ylabel(L"\alpha")
figure()
surf(betaVec, alphaVec, c)
xlabel(L"\beta")
ylabel(L"\alpha")
