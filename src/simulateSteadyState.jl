include("finiteDifferences.jl")
using PyPlot
close("all")

v_0 = 0.08 # Amplitude of the oscillations in the potential
ff = 3.0 # Forcing for the potential
LL = 1.0 # Length of one period
T0 = 17.0 # Temperature at the ends (i.e. bath temperature)

alpha = 0.00001
beta = 3.0
heat_capacity = 15.0

n_periods = 6 # The number of periods to stretch out for
n_points = 600 # The number of points on the x_axis
dt = 0.0001 # Size of the time step (keep this much smaller than the grid
             # spacing)

# x_axis = linspace(-LL*n_periods, LL*n_periods, n_points)
x_axis = linspace(5LL, 10LL, n_points)
dx = (x_axis[end] - x_axis[1])/n_points # Grid spacing

# V(x, time) = ff*x^2 + 6ff # Potential
sigma = 0.1
V(x) = 0.1x^4 + 15ff*exp(-((x-9)^2)/(2sigma^2)) + ff
dis_V = Float64[V(x) for x in x_axis]
dis_V0, dis_V_end = V(x_axis[1] - dx), V(x_axis[end] + dx)
dis_V_tup = (dis_V0, dis_V, dis_V_end)

# temp_fun(x) = (T0 + 0.1*sin(2*(2pi/LL)*x))
temp_fun(x) = T0
dis_temp = Float64[temp_fun(x) for x in x_axis]
dis_temp0, dis_temp_end = temp_fun(x_axis[1] - dx), temp_fun(x_axis[end] + dx)
dis_temp_tup = (dis_temp0, dis_temp, dis_temp_end)

sigma  = 0.01
P_0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x-9.3)^2)/(2sigma^2))
              for x in x_axis]
# P_0 = ones(x_axis)
P_0 /= discrete_quad(P_0, x_axis[1], x_axis[end])
density = P_0
# Calculate the initial energy of the system.
energy = energy_fun(dis_V, P_0, dis_temp, x_axis)
boltzmann_density = exp(-dis_V/T0)
boltzmann_density = boltzmann_density/discrete_quad(boltzmann_density,
                x_axis[1], x_axis[end])
evolve_time = 1.5
@time begin
density, temp = evolve_system(density, dis_temp, evolve_time, dt, dis_V_tup,
                alpha, beta, heat_capacity, energy, x_axis)
density2 = evolveP(P_0, evolve_time, dt, dis_V_tup, ones(x_axis)*T0, x_axis)
end

plot(x_axis, 0.02dis_V, x_axis, density, x_axis, density2, x_axis,
            temp, x_axis, ones(x_axis)*T0)
legend(["Potential", "Coupled probability density",
        "Uncoupled probability density", "Temperature", "Initial temperature"])
