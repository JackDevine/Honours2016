# Given the steady state functions, we want to know how the system behaves for
# different values of the parameters.

@everywhere include("SteadyState.jl")

T0 = 300
kappa = 0.3
D = 0.025
f = 1.0
L = 1.0
k_B = 0.01
v_0 = 1.0

n = 2
# current = Array(Float64, n)
kappa_vec = linspace(0.01, 1.2, n)
D_vec = linspace(0.01, 1.2, n)
# count = 1
# @time begin
#     for kappa in kappa_vec
#         current[count], P0 = steady_params(T0, kappa, D, f)
#         count += 1
#     end
# end

current_ser = Array(Float64, n)
@time current_ser = map(kappa -> steady_params(T0, kappa, D, f)[1], kappa_vec)

current_par = SharedArray(Float64, n)
@time current_par = pmap(kappa -> steady_params(T0, kappa, D, f)[1], kappa_vec)

current_par2 = SharedArray(Float64, n)

@time begin
    @parallel for i = 1:n
        current_par2[i] = steady_params(T0, kappa_vec[i], D, f)[1]
    end
end

@time [steady_params(T0, kappa_vec[i], D_vec[j], f)
                        for i in 1:n, j in 1:n]

@time @parallel [steady_params(T0, kk, DD, f)
                        for kk in kappa_vec, DD in D_vec]
