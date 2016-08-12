# @everywhere include("simulateSteadyState.jl")
# @everywhere using Hopping
using PyPlot
close("all")

@everywhere begin
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
