using FactCheck
include("finiteDifferences.jl")

# Initialize a particular starting system.
v_0 = 0.8 # Amplitude of the oscillations in the potential
ff = 3.0 # Forcing for the potential
LL = 1.0 # Length of one period
T0 = 1.0 # Temperature at the ends (i.e. bath temperature)

alpha = 0.0004
beta = 0.15

n_periods = 4 # The number of periods to stretch out for
n_points = 2000 # The number of points on the x_axis
n_steps = 2000 # The number of time steps to simulate for
dt = 0.00005 # Size of the time step (keep this much smaller than the grid
             # spacing)

x_axis = linspace(0, LL*n_periods, n_points)
dx = (x_axis[end] - x_axis[1])/n_points # Grid spacing

V(x, time) = ff*x - v_0*sin(2*(2pi/LL)*x) + 6ff # Potential
dis_V = Float64[V(x, 0) for x in x_axis]
dis_V_minus1, dis_V0, dis_V_end = V(x_axis[1] - 2dx, 0), V(x_axis[1] - dx, 0),
                                    V(x_axis[end] + dx, 0)
dis_V_tup = (dis_V0, dis_V, dis_V_end)

temp_fun(x) = (T0 + 0.1*sin(2*(2pi/LL)*x))
dis_temp = Float64[temp_fun(x) for x in x_axis]
dis_temp0, dis_temp_end = temp_fun(x_axis[1] - dx), temp_fun(x_axis[end] + dx)
dis_temp_tup = (dis_temp0, dis_temp, dis_temp_end)

sigma  = 0.1
P_0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x-3.0)^2)/(2sigma^2))
              for x in x_axis]
P_0 /= discrete_quad(P_0, x_axis[1], x_axis[end])

energy = discrete_quad(dis_V.*P_0, x_axis[1], x_axis[end])
               + discrete_quad(dis_temp, x_axis[1], x_axis[end])
# Evolve the system using different time steps, for small time steps, halfing
# the time step should not cause the result of the simulation to change by much.
tol = 1e-3  # The amount that the two results are allowed to differ by.
evolve_time = 0.5  # The amount of time that we will evolve the system by.
time_steps = evolve_time*[1/2048, 1/4096]  # Time steps used in the simulations.
facts("Convergence tests.") do
    @fact (norm(evolveT(dis_temp, evolve_time, time_steps[1], dis_V_tup, P_0,
    alpha, beta, energy, x_axis)
    - evolveT(dis_temp, evolve_time, time_steps[2], dis_V_tup, P_0, alpha, beta,
    energy, x_axis))
    --> less_than(tol) ) "Temperature evolution is not converging."

    @fact (norm(evolveP(P_0, evolve_time, time_steps[1], dis_V_tup, dis_temp,
                                x_axis)
    - evolveP(P_0, evolve_time, time_steps[2], dis_V_tup, dis_temp, x_axis) )
    --> less_than(tol) ) "Probability evolution is not convering."
    # Evolve the system forward and check that the normed differences of the
    # tempearture and the probability density are changing by less than the
    # tolerance.
    system_time_step1 = evolve_system(P_0, dis_temp, evolve_time, time_steps[1],
                                dis_V_tup, alpha, beta, energy, x_axis)

    system_time_step2 = evolve_system(P_0, dis_temp, evolve_time, time_steps[2],
                                dis_V_tup, alpha, beta, energy, x_axis)
    @fact (norm(system_time_step1[1] - system_time_step2[1])
            --> less_than(tol) ) """The temperature in the system evolution is
                                    not converging."""
    @fact (norm(system_time_step1[2] - system_time_step2[2])
            --> less_than(tol) ) """The probability density in the system
                                    evolution is not converging."""
end
