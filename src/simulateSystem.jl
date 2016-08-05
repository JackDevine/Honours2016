# Simulate the system for some given parameters
include("finiteDifferences.jl")
using PyPlot
close("all")

v_0 = 0.08 # Amplitude of the oscillations in the potential
ff = 3.0 # Forcing for the potential
LL = 1.0 # Length of one period
T0 = 3.0 # Temperature at the ends (i.e. bath temperature)

alpha = 0.0004
beta = 0.15

n_periods = 6 # The number of periods to stretch out for
n_points = 600 # The number of points on the x_axis
n_steps = 500 # The number of time steps to simulate for
dt = 0.0005 # Size of the time step (keep this much smaller than the grid
             # spacing)

x_axis = linspace(-LL*n_periods, LL*n_periods, n_points)
dx = (x_axis[end] - x_axis[1])/n_points # Grid spacing

V(x, time) = ff*x^2 - v_0*sin(2*(2pi/LL)*x) + 6ff # Potential
dis_V = Float64[V(x, 0) for x in x_axis]
dis_V0, dis_V_end = V(x_axis[1] - dx, 0), V(x_axis[end] + dx, 0)
dis_V_tup = (dis_V0, dis_V, dis_V_end)

# temp_fun(x) = (T0 + 0.1*sin(2*(2pi/LL)*x))
temp_fun(x) = T0
dis_temp = Float64[temp_fun(x) for x in x_axis]
dis_temp0, dis_temp_end = temp_fun(x_axis[1] - dx), temp_fun(x_axis[end] + dx)
dis_temp_tup = (dis_temp0, dis_temp, dis_temp_end)

sigma  = 0.1
P_0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x-2.0)^2)/(2sigma^2))
              for x in x_axis]
# P_0 = ones(x_axis)
P_0 /= discrete_quad(P_0, x_axis[1], x_axis[end])
density = P_0
# Calculate the initial energy of the system.
energy = energy_fun(dis_V, P_0, dis_temp, x_axis)
energy_vec = Array(Float64, n_steps)
energy_vec[1] = energy
potential_energy = Array(Float64, n_steps)
potential_energy[1] = discrete_quad(dis_V.*density, x_axis[1], x_axis[end])
thermal_energy = Array(Float64, n_steps)
thermal_energy[1] = discrete_quad(dis_temp, x_axis[1], x_axis[end])
change = Array(Float64, n_steps-1)
density_new = Array(Float64, length(x_axis))
@time begin
    for i = 2:n_steps
        t = i*dt

        # Use this code if you have a time dependant potential
 #         dis_V =
        # Float64[V(x, time_scale*t) for x in x_axis]
 #         dis_V_arr[:, i] =
        # dis_V
 #         dis_V_minus1, dis_V0, dis_V_end = V(x_axis[1] - 2dx,
        # 0), V(x_axis[1] - dx, 0), V(x_axis[end] + dx, 0)
 #         dis_V_tup =
        # (dis_V_minus1, dis_V0, dis_V, dis_V_end)
        density = stepP(density, dt, dis_V_tup, dis_temp, x_axis)
        dis_temp = stepT(dis_temp, dt, density, dis_V_tup, alpha,
                            beta, energy, x_axis)
        energy_vec[i] = energy_fun(dis_V, density, dis_temp, x_axis)
        # energy_vec[i] = discrete_quad(density.*dis_V, x_axis[1], x_axis[end]) + discrete_quad(dis_temp, x_axis[1], x_axis[end])
        potential_energy[i] = discrete_quad(density.*dis_V, x_axis[1], x_axis[end])
        thermal_energy[i] = discrete_quad(dis_temp, x_axis[1], x_axis[end])
    end
end
# plot(1:n_steps, energy_vec, 1:n_steps, potential_energy, 1:n_steps, thermal_energy)
# legend(["Energy", "pot", "therm"])
plot(x_axis, 0.1dis_V, x_axis, density, x_axis, dis_temp)
