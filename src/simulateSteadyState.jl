include("finiteDifferences.jl")
using PyPlot
close("all")

v_0 = 0.08 # Amplitude of the oscillations in the potential
ff = 3.0 # Forcing for the potential
LL = 1.0 # Length of one period
T0 = 7.0 # Temperature at the ends (i.e. bath temperature)

alpha = 0.0004
beta = 0.1

n_periods = 6 # The number of periods to stretch out for
n_points = 600 # The number of points on the x_axis
dt = 1e-5 # Size of the time step (keep this much smaller than the grid
             # spacing)

# x_axis = linspace(-LL*n_periods, LL*n_periods, n_points)
x_axis = linspace(5LL, 10LL, n_points)
dx = (x_axis[end] - x_axis[1])/n_points # Grid spacing

# V(x, time) = ff*x^2 + 6ff # Potential
sigma = 0.1
bump = 9
V(x) = 0.3x^3 + 3.0ff*exp(-((x-bump)^2)/(2sigma^2))
dis_V = Float64[V(x) for x in x_axis]
dis_V0, dis_V_end = V(x_axis[1] - dx), V(x_axis[end] + dx)
dis_V_tup = (dis_V0, dis_V, dis_V_end)

# temp_fun(x) = (T0 + 0.1*sin(2*(2pi/LL)*x))
temp_fun(x) = T0
dis_temp = Float64[temp_fun(x) for x in x_axis]
dis_temp0, dis_temp_end = temp_fun(x_axis[1] - dx), temp_fun(x_axis[end] + dx)
dis_temp_tup = (dis_temp0, dis_temp, dis_temp_end)

sigma  = 0.05
P_0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x-9.3)^2)/(2sigma^2))
              for x in x_axis]
# P_0 = ones(x_axis)
P_0 /= discrete_quad(P_0, x_axis[1], x_axis[end])
density = P_0
# Calculate the initial energy of the system.
energy = energy_fun(dis_V, P_0, dis_temp, alpha, x_axis)
boltzmann_density = exp(-dis_V/T0)
boltzmann_density = boltzmann_density/discrete_quad(boltzmann_density,
                x_axis[1], x_axis[end])
evolve_time = 1000dt
@time begin
density, temp = evolve_system(density, dis_temp, evolve_time, dt, dis_V_tup,
                alpha, beta, energy, x_axis)
density2 = evolveP(P_0, evolve_time, dt, dis_V_tup, ones(x_axis)*T0, x_axis)
end

plot(x_axis, 0.02dis_V, x_axis, density, x_axis, density2, x_axis,
            temp, x_axis, ones(x_axis)*T0)
legend(["Potential", "Coupled probability density",
        "Uncoupled probability density", "Temperature", "Initial temperature"])

function hopping_time(dis_V_tup::Tuple{Number, AbstractArray, Number}, bump::Number,
                density::AbstractArray, temperature::AbstractArray,
                alpha::Number, beta::Number, heat_capacity::Number,
                x_axis::AbstractArray; dt=1e-4, tol=0.1)
    # Given an initial state of the system and a potential with a bump in it,
    # calculate how long it takes for the probability distribution to get over
    # the bump.
    # Parameters:
    # dis_V:         An array containg the potential.
    # bump:          The location of the bump in the potetnial.
    # density:       The initial density of the system.
    # temperature:   The initial temperature of the system.
    # alpha:         Dimensionless parameter that controls the coupling of the
    #                system.
    # beta:          Dimensionless parameter that determines how fast the
    #                temperature gradients diffuse.
    # heat_capacity: Heat capacity of the environment (dimensionless).
    # x_axis:        The axis that we are working on (dimensionless).
    # dt:            Size of one time step for simulation.
    # tol:           If the fraction of the probability density that is to the
    #                left of the bump is less than tol, then we will consider
    #                the particle to have crossed the barrier.

    # Calculate the location of the center of the initial distribution.
    bump_ind = indmin(abs(x_axis - bump))
    iters = 0
    # while discrete_quad(density[bump_ind:end], bump, x_axis[end]) > tol
    while discrete_quad(density.*x_axis, x_axis[1], x_axis[end]) > bump
        density = stepP(density, dt, dis_V_tup, temperature, x_axis)
        temperature = stepT(temperature, dt, density, dis_V_tup, alpha,
                            beta, heat_capacity, energy, x_axis)
        iters += 1
        # println(discrete_quad(density.*x_axis, x_axis[1], x_axis[end]))
        if iters > 2000
            println("alpha = $alpha , beta = $beta DNF")
            return dt*iters
        end
    end
    ret = dt*iters
end
#
# hopping_time(dis_V_tup, bump, P_0, dis_temp, alpha, beta, heat_capacity,
#                     x_axis; dt=1e-4, tol=0.7)
nPoints = 50
alphaMin = 0.0002
alphaMax = 0.04
betaMin = 0.1
betaMax = 3.0
betaVec = linspace(betaMin, betaMax, nPoints)
alphaVec = linspace(alphaMin, alphaMax, nPoints)
# @time a = Float64[hopping_time(dis_V_tup, bump, P_0, dis_temp, alpha, beta,
#                     heat_capacity, x_axis)
#                     for alpha in linspace(alphaMin, alphaMax, nPoints),
#                     beta in linspace(betaMin, betaMax, nPoints)]
# @time a = map((alpha, beta) -> hopping_time(dis_V_tup, bump, P_0, dis_temp, alpha, beta,
#                     heat_capacity, x_axis),
#                     linspace(0.0, 0.8, nPoints),
#                     linspace(1, 4, nPoints))

# matshow(a)
# xlabel(L"\beta")
# ylabel(L"\alpha")
# figure()
# surf(betaVec, alphaVec, a)
# xlabel(L"\beta")
# ylabel(L"\alpha")
