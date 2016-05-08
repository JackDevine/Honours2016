# Given the steady state functions, we want to know how the system behaves for
# different values of the parameters.

include("SteadyState.jl")

n = 25
current = Array(Float64, n)
kappa_vec = linspace(0.01, 1.2, n)
count = 1
@time begin
    for kappa in kappa_vec
        current[count], P0 = steady_params(T0, kappa, D, f)
        count += 1
    end
end
