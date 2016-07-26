# Simulate the system for some given parameters
include("finiteDifferences.jl")
using PyPlot

v_0 = 0.08 # Amplitude of the oscillations in the potential
ff = 3.0 # Forcing for the potential
LL = 1.0 # Length of one period
T0 = 3.0 # Temperature at the ends (i.e. bath temperature)

alpha = 0.0004
beta = 0.15

n_periods = 4 # The number of periods to stretch out for
n_points = 5000 # The number of points on the x_axis
n_steps = 300 # The number of time steps to simulate for
dt = 0.0005 # Size of the time step (keep this much smaller than the grid
             # spacing)

x_axis = linspace(-LL*n_periods, LL*n_periods, n_points)
dx = (x_axis[end] - x_axis[1])/n_points # Grid spacing

V(x, time) = ff*x^2 - v_0*sin(2*(2pi/LL)*x) + 6ff # Potential
dis_V = Float64[V(x, 0) for x in x_axis]
dis_V_minus1, dis_V0, dis_V_end = V(x_axis[1] - 2dx, 0), V(x_axis[1] - dx, 0),
                                    V(x_axis[end] + dx, 0)
dis_V_tup = (dis_V0, dis_V, dis_V_end)

# temp_fun(x) = (T0 + 0.1*sin(2*(2pi/LL)*x))
temp_fun(x) = T0
dis_temp = Float64[temp_fun(x) for x in x_axis]
dis_temp0, dis_temp_end = temp_fun(x_axis[1] - dx), temp_fun(x_axis[end] + dx)
dis_temp_tup = (dis_temp0, dis_temp, dis_temp_end)

sigma  = 0.1
P_0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x-1.0)^2)/(2sigma^2))
              for x in x_axis]
# P_0 = ones(x_axis)
P_0 /= discrete_quad(P_0, x_axis[1], x_axis[end])
temp = Array(Float64, n_points)
P = Array(Float64, n_points, n_steps)
P[:, 1] = P_0
# Calculate the initial energy of the system.
energy = energy_fun(dis_V, P_0, dis_temp, x_axis)
T = Array(Float64, n_points, n_steps)
T[:, 1] = dis_temp

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

        P[:, i] = stepP(P[:, i-1], dt, dis_V_tup, dis_temp, x_axis)

        dis_temp = stepT(dis_temp, dt, P[:, i-1], dis_V_tup, alpha,
                            beta, energy, x_axis)
        T[:, i] = dis_temp
    end
end

n = 1
time_step = n_steps
plot(x_axis, 0.1dis_V, x_axis, P[:, time_step], x_axis, T[:, time_step])
