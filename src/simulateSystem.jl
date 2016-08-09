# Simulate the system for some given parameters
include("finiteDifferences.jl")
using PyPlot
close("all")

v_0 = 0.08 # Amplitude of the oscillations in the potential
ff = 3.0 # Forcing for the potential
LL = 1.0 # Length of one period
T0 = 3.0 # Temperature at the ends (i.e. bath temperature)

alpha = 0.004
beta = 0.015

nPeriods = 6 # The number of periods to stretch out for
nPoints = 600 # The number of points on the x_axis
nSteps = 1000 # The number of time steps to simulate for
dt = 0.0005 # Size of the time step (keep this much smaller than the grid
             # spacing)

x_axis = linspace(-LL*nPeriods, LL*nPeriods, nPoints)
dx = (x_axis[end] - x_axis[1])/nPoints # Grid spacing

V(x, time) = ff*x^2 - v_0*sin(2*(2pi/LL)*x) + 6ff # Potential
potential = Float64[V(x, 0) for x in x_axis]
potential0, potentialEnd = V(x_axis[1] - dx, 0), V(x_axis[end] + dx, 0)
potentialTup = (potential0, potential, potentialEnd)

# temperatureFun(x) = (T0 + 0.1*sin(2*(2pi/LL)*x))
temperatureFun(x) = T0
temperature = Float64[temperatureFun(x) for x in x_axis]
temperature0, temperatureEnd = temperatureFun(x_axis[1] - dx),
                temperatureFun(x_axis[end] + dx)
temperatureTup = (temperature0, temperature, temperatureEnd)

sigma  = 0.1
P0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x-2.0)^2)/(2sigma^2))
              for x in x_axis]
# P0 = ones(x_axis)
P0 /= discrete_quad(P0, x_axis[1], x_axis[end])
density = P0
# Calculate the initial energy of the system.
energy = energyFun(potential, P0, temperature, alpha, x_axis)
energyVec = Array(Float64, nSteps)
energyVec[1] = energy
potentialEnergy = Array(Float64, nSteps)
potentialEnergy[1] = discrete_quad(potential.*density, x_axis[1], x_axis[end])
thermalEnergy = Array(Float64, nSteps)
thermalEnergy[1] = discrete_quad(temperature, x_axis[1], x_axis[end])
@time begin
    for i = 2:nSteps
        t = i*dt

        # Use this code if you have a time dependant potential
 #         potential =
        # Float64[V(x, time_scale*t) for x in x_axis]
 #         potential_arr[:, i] =
        # potential
 #         potential_minus1, potential0, potential_end = V(x_axis[1] - 2dx,
        # 0), V(x_axis[1] - dx, 0), V(x_axis[end] + dx, 0)
 #         potentialTup =
        # (potential_minus1, potential0, potential, potential_end)
        density = stepP(density, dt, potentialTup, temperature, x_axis)
        temperature = stepT(temperature, dt, density, potentialTup, alpha,
                            beta, energy, x_axis)
        energyVec[i] = energyFun(potential, density, temperature, alpha, x_axis)
        # energyVec[i] = discrete_quad(density.*potential, x_axis[1], x_axis[end]) + discrete_quad(temperature, x_axis[1], x_axis[end])
        potentialEnergy[i] = discrete_quad(density.*potential,
                            x_axis[1], x_axis[end])
        thermalEnergy[i] = discrete_quad(temperature, x_axis[1], x_axis[end])
    end
end
# plot(1:nSteps, energyVec, 1:nSteps, potentialEnergy, 1:nSteps, thermalEnergy)
# legend(["Energy", "pot", "therm"])
plot(x_axis, 0.1potential, x_axis, density, x_axis, temperature)
