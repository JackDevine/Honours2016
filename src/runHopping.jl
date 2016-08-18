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

xAxis = linspace(-2LL, 2LL, nPoints)
dx = (xAxis[end] - xAxis[1])/nPoints # Grid spacing

sigma = 0.1
bump = 0.0
V(x) = x^4 - 3*(exp(-(x - 0.5)^2/sigma) + exp(-(x + 0.5)^2/sigma)) + 3x
potential = Float64[V(x) for x in xAxis]
potential0, potentialEnd = V(xAxis[1] - dx), V(xAxis[end] + dx)
potentialTup = (potential0, potential, potentialEnd)

# temperatureFun(x) = (T0 + 0.1*sin(2*(2pi/LL)*x))
temperatureFun(x) = T0
temperature = Float64[temperatureFun(x) for x in xAxis]

sigma  = 0.05
P0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x-0.5)^2)/(2sigma^2))
              for x in xAxis]
# P0 = ones(xAxis)
P0 /= discrete_quad(P0, xAxis[1], xAxis[end])
density = P0
# Calculate the initial energy of the system.
energy = energyFun(potential, P0, temperature, alpha, xAxis)
boltzmann_density = exp(-potential/T0)
boltzmann_density = boltzmann_density/discrete_quad(boltzmann_density,
                   xAxis[1], xAxis[end])

nPoints = 5
alphaMin = -0.00001
alphaMax = -0.00005
betaMin = 0.001
betaMax = 0.5
betaVec = linspace(betaMin, betaMax, nPoints)
alphaVec = linspace(alphaMin, alphaMax, nPoints)
alphaBeta = [[alpha, beta] for alpha in alphaVec, beta in betaVec]
# Make a method of run_hopping for this particular potential.
hopping_time(alpha, beta) = hopping_time(potentialTup, bump, density,
 temperature, alpha, beta, xAxis; dt=1e-4)
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
