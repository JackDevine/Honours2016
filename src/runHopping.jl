# @everywhere include("simulateSteadyState.jl")
# @everywhere using Hopping
close("all")

@everywhere begin
nPoints = 500
alphaMin = -0.0003
alphaMax = -0.00029
betaMin = 0.000001
betaMax = 1.0
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
c = convert(Array{Float64, 2}, c)

matshow(c)
title("Hopping time for a bistable potential")
xlabel(L"\beta")
ylabel(L"\alpha")
figure()
surf(betaVec, alphaVec, c)
xlabel(L"\beta")
ylabel(L"\alpha")
figure()
surf(log(c))
figure()
matshow(log(c), extent=[betaVec[1], betaVec[end], 1e5*alphaVec[1],
 1e5*alphaVec[end]])
colorbar()
title("Hopping time for a bistable potential (log scale)")
xlabel(L"\beta")
ylabel(L"\alpha \times 10^5")
# @time c = pmap(alpha -> hopping_time(alpha, 0.05), alphaVec)
# c = convert(Array{Float64, 1}, c)
# plot(alphaVec, c)
