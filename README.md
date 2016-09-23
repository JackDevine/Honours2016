# Honours2016
Honours project on Brownian dynamics, the source code is in the `src` directory and includes unit tests for the important functions in the project, the most important file is called finiteDifferences.jl, here you will find my implementation pf finite differences.

In the `notebooks` directory, you will find some notebooks that walk through examples on how to use the code and some important results from the project.

If you want to run the code on your own machine, then you will need to install [Julia](http://julialang.org/), once you have Julia, there are some packages that you will need to run the code, simply go to the Julia command line and type.

```julia
julia> Pkg.add("PyPlot")  # Used for plotting.
julia> Pkg.add("ForwardDiff")  # Used for automatically finding the derivatives of functions.
julia> Pkg.add("SymPy")  # Used to obtain analytical solutions to the Smoluchowski equation.
```

Once you have these prerequisites, simply clone this repo and start playing. You may need to change some paths in some of the notebooks, but that won't be too hard.

# Usage
The main point of the code in this repository is to perform finite differencing, this is done with two main functions called `stepP` and `stepT`. These functions perform finite differences representing the probability density or the temperature respectively. I have made a user type called `System` which contains the data for a system. You can use this type as an argument to either `stepP` or `stepT`, an example is shown below
```julia
using FiniteDifferences
using ForwardDiff
using PyPlot

# Demonstration of periodic boundary conditions.
alpha = 5e-3
beta = 1e-2
T0 = 1.0
nSteps = 15000  # The number of steps in time that we will take.
nPoints = 1000
dt = 1e-4
xAxis = linspace(0, 5, nPoints)
dx = (xAxis[end] - xAxis[1])/nPoints
V(x) = 0.5sin(2pi*x) + 0.2cos(2*2pi*x) + 0.8  # The potential for the system.
dV(x) = ForwardDiff.derivative(V, x)  # The derivative of the potential.

# Make the initial probability distribution.
sigma  = 0.4
P0 = Float64[(1/(sigma*sqrt(2pi)))*exp(-((x - 3)^2)/(2sigma^2))
              for x in xAxis]
P0 /= discrete_quad(P0, xAxis[1], xAxis[end])
density = P0

temperature = T0*ones(nPoints)

# Initialize the system.
potential = Float64[V(x) for x in xAxis]
dpotential = [dV(xAxis[1] - dx) ; Float64[dV(x) for x in xAxis] ; dV(xAxis[end] + dx)]

system = System()
system.density = density
system.potential = potential
system.dpotential = dpotential
system.temperature = temperature
system.xAxis = xAxis
system.energy = energyFun(system, alpha)

# Make some arrays for storing the results as we go.
P = Array(Float64, nPoints, nSteps)
T = Array(Float64, nPoints, nSteps)
P[:, 1] = system.density
T[:, 1] = system.temperature

for i in 2:nSteps
    system.temperature = stepT(system, alpha, beta, dt; bndType = :neumann)
#     system.energy, system.temperature = stepT(system, alpha, beta, dt; bndType = :dirichlet)
    system.density = stepP(system, dt, bndType = :periodic)
    
    P[:, i] = system.density
    T[:, i] = system.temperature
end

plot(xAxis, 0.5system.potential, linewidth = 2)
plot(xAxis, 0.8T0*ones(xAxis), color = "red", linewidth = 2)
plot(xAxis, P0, color = "black", linewidth = 2)
fill_between(xAxis, P0, color = "black", alpha = 0.2)
title("Initial configuration")
xlabel(L"x / L", fontsize = 30)
ylabel("Potential, Probability, temperature", fontsize = 20)
legend(["Potential", "Temperature", "Probaility"], loc=2, frameon = false)
ylim([0, 1.1])

figure()
plot(xAxis, 0.5system.potential, linewidth = 2)
plot(xAxis, 0.8system.temperature, color = "red", linewidth = 2)
plot(xAxis, system.density, color = "black", linewidth = 2)
fill_between(xAxis, system.density, color = "black", alpha = 0.2)
title("An example of periodic boundary conditions")
xlabel(L"x / L", fontsize = 30)
ylabel("Potential, Probability, temperature", fontsize = 30)
ylim([0, 1.1])
legend(["Potential", "Temperature", "Probaility"], loc=2)
```
This will result in the Figure A.1 in the thesis 
