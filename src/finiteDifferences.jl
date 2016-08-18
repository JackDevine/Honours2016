# Useful functions for doing finite differences on the system, all of these
# functions use the dimensionalized equations.
"""
    discrete_quad(vec::AbstractArray, start::Number, fin::Number)
Do discrete quadrature on a vector vec where the points in vec are evenly
spaced from start to fin, if vec is even, use the Simpson rule, otherwise use
the trapezoid rule.
# Examples
```julia
julia> discrete_quad(sin(linspace(0, pi, 1000)), 0, pi)
julia> 1.9979983534190815
```
"""
function discrete_quad(vec::AbstractArray, start::Number, fin::Number)
    nn = length(vec)  # Number of points that we are discretizing over.
    h = (fin - start)/nn
    if length(vec) % 2 == 1
        # Use the trapezoid rule.
        return (h/2)*(vec[1] + vec[end] + 2*sum(vec[2:end-1]))
    end
    # Otherwise use the simpson rule.
    s = vec[1] + vec[end]
    for i in 1:2:nn
        s += 4*vec[i]
    end
    for i in 2:2:nn-1
        s += 2*vec[i]
    end
    s*h/3
end

"""
    discrete_derivative(vec::AbstractArray, xAxis::AbstractArray)
Calculate the derivative of vec along the axis xAxis, note that the
vector returned is one element shorter than the vectors given.
"""
function discrete_derivative(vec::AbstractArray, xAxis::AbstractArray)
    diff(vec)./diff(xAxis)
end
"""
    energyFun(potential::AbstractArray, density::AbstractArray,
            temperature::AbstractArray, alpha::Number,
            xAxis::AbstractArray)
Calcualte the energy of the sytem.
"""
function energyFun(potential::AbstractArray, density::AbstractArray,
            temperature::AbstractArray, alpha::Number,
            xAxis::AbstractArray)
    left_bnd = xAxis[1]
    right_bnd = xAxis[end]
    (discrete_quad(potential.*density, left_bnd, right_bnd)
        + (1/alpha)*discrete_quad(temperature, left_bnd, right_bnd))
end
"""
    stepP(P::AbstractArray, dt::Number,
            potentialTup::Tuple{Number, AbstractArray, Number},
            temperature::AbstractArray, xAxis::AbstractArray)
Evolve the probability density P forward by an amount dt.
# Arguments:
* `P::AbstractArray`: Initial probability density, this is a vector of the
same length as the xAxis (see below).
* `dt::Number`: Amount of time to simulate forward by.
* `potentialTup::Tuple{Number, AbstractArray, Number}`:  A tuple containing
information on the potential, the first element is the potential to the left of
the left boundary, the second element is a vector containing the values of the
potential on the xAxis,the third element is the potential to the right of the
right boundary.
* `temperature::AbstractArray`: A vector of the tempeature.
* xAxis: A vector containing the x coordinates of the points that we are doing
finite differencing on, the points must be equally spaced.
Returns the probability density evolved forward by an amount dt, assumes
periodic boundary conditions.
* `bndType::Symbol`: Specify the solution or the derivative at the boudary,
    * `:dirichlet`: The density remains consant at the boudaries.
    * `:periodic`: The density is the same at the left and right sides.
    * `:neumann`: The derivative is zero at the boundaries.
"""
function stepP(P::AbstractArray, dt::Number,
            potentialTup::Tuple{Number, AbstractArray, Number},
            temperature::AbstractArray, xAxis::AbstractArray,
            bndType::Symbol=:dirichlet)
    # Augment the discrete temperature and the discrete potential since we will
    # be evaluating them at points beyond the boundary.
    T0 = temperature[1]
    # The derivative of the temperature is zero at the boundaries.
    temperature = [temperature[1] ; temperature ; temperature[end]]

    dis_V0, dis_V, dis_V_end = potentialTup
    dis_V = [dis_V0 ; dis_V ; dis_V_end]

    n_points = length(xAxis)
    dx = (xAxis[end] - xAxis[1])/n_points  # Grid spacing
    rr = dt/(2dx^2)  # Useful constant for making the matrix
    # First we will make a matrix A that represents the forward Euler step.
    # Create the diagonals of the matrix, in order to linearize the problem, we
    # had to assume that T(x, t + dt) = T(x, t). We could improve on this by
    # using an explicit method to estimate T(x, t + dt).
    diag_minus1 = rr*(0.25*(dis_V[3:end-1] - dis_V[1:end-3])
                 - temperature[1:end-3])
    diag0 = rr*(-(dis_V[3:end] - 2dis_V[2:end-1] + dis_V[1:end-2])
           + 2*temperature[2:end-1]) + 1
    diag1 = rr*(-0.25*(dis_V[4:end] - dis_V[2:end-2]) - temperature[4:end])
    if bndType == :dirichlet
        A = Tridiagonal(diag_minus1, diag0, diag1)
        B = Tridiagonal(-diag_minus1, 2 - diag0, -diag1)
    elseif bndType == :periodic
        A = spdiagm((diag_minus1, diag0, diag1), (-1, 0, 1))
        A[1, end] = diag_minus1[1]
        A[end, 1] = diag1[end]
        B = 2speye(size(A)...) - A
    elseif bndType == :neumann
        # Derivative is zero at the boundaries.
        A = spdiagm((diag_minus1, diag0, diag1), (-1, 0, 1))
        B = 2speye(size(A)...) - A
        A[1, 2] = -A[1, 1]
        A[end, end - 1] = -A[end, end]
        B[1, 2] = -B[1, 1]
        B[end, end - 1] = -B[end, end]
    end
    P = A\(B*P)  # Since A*P^{n+1} = B*P^n
    # Return the normalized probability density
    P = P/discrete_quad(P, xAxis[1], xAxis[end])
end
"""
Evolve the temperature forward by an amount dt.
# Arguments:
* `temperature::AbstractArray`:   Initial temperature, this is a vector of the
same length as the xAxis (see below).
* `dt::Number`:            Amount of time to simulate forward by.
* `dis_density::AbstractArray`:   A vector of the probability density,
must be the same length as the xAxis (see below).
* `potentialTup::Tuple{Number, AbstractArray, Number}`: A tuple containing
information on the potential, the first element is the potential to the left of
the left boundary, the second element is a vector containing the values of the
potential on the xAxis, the third element is the potential to the right of the
right boundary, the third element is the potential to the right of the right
boundary.
* `alpha::Number`: Dimensionless parameter that describes how much the heat
from the motor affects the temperature.
* `beta::Number`: Dimensionless parameter that describes how quickly the
temperature diffuses to a constant value.
* `init_energy::Number`:   The initial energy of the system, must be in
dimensionless units, i.e. the actual energy divided by E0, where E0 is the
potential energy gained by moving one period up the potential (characteristic
energy of the system).
* `xAxis::AbstractArray`: A vector containing the x coordinates of the points
that we are doing finite differencing on, the points must be equally spaced.
Returns the temperature evolved forward by an amount dt as well as the
updated energy of the system, assumes periodic boundary conditions.
* `bndType::Symbol=:neumann`: The type of boundary condition for the temperature
can includes `:neumann`, `:dirichlet`, `:periodic`.
"""
function stepT(temperature::AbstractArray, dt::Number,
            dis_density::AbstractArray,
            potentialTup::Tuple{Number, AbstractArray, Number}, alpha::Number,
            beta::Number, energy::Number,
            xAxis::AbstractArray, bndType::Symbol=:neumann)
    # Augment the discrete density and the discrete potential since we will be
    # evaluating them at points beyond the boundary.
    dis_V0, dis_V, dis_V_end = potentialTup
    dis_V = [dis_V0 ; dis_V ; dis_V_end]
    # The density is periodic.
    dis_density = [dis_density[end] ; dis_density ; dis_density[1]]

    n_points = length(xAxis)
    dx = (xAxis[end] - xAxis[1])/n_points
    rr = dt/(2dx^2)
    # The inhomogeniety at the end of the equation.
    in_homo = rr*(alpha/4)*dis_density[2:end-1].*(dis_V[3:end] - dis_V[1:end-2])
    # The diagonals of the matrix.
    diag_minus1 = rr*((alpha/4)*dis_density[1:end-3]
                .*(dis_V[3:end-1] - dis_V[1:end-3]) - beta)
    diag0 = (2*beta*rr + 1)*ones(xAxis)
    diag1 = rr*(-(alpha/4)*dis_density[4:end]
            .*(dis_V[4:end] - dis_V[3:end-1]) - beta)
    # Add the inhomogeneity to the diagonals.
    diag_minus1 += in_homo[2:end]
    diag0 += in_homo
    diag1 += in_homo[1:end-1]

    if bndType == :neumann
        # Derivative is zero at the boundaries.
        A = spdiagm((diag_minus1, diag0, diag1), (-1, 0, 1))
        B = 2speye(size(A)...) - A
        A[1, 2] = -A[1, 1]
        A[end, end - 1] = -A[end, end]
        B[1, 2] = -B[1, 1]
        B[end, end - 1] = -B[end, end]
    elseif bndType == :periodic
        A = spdiagm((diag_minus1, diag0, diag1), (-1, 0, 1))
        A[1, end] = diag_minus1[1]
        A[end, 1] = diag1[end]
        B = 2speye(size(A)...) - A
    elseif bndType == :dirichlet
        A = Tridiagonal(diag_minus1, diag0, diag1)
        B = Tridiagonal(-diag_minus1, 2 - diag0, -diag1)
    end
    # Update the temperature using the Crank Nicolson scheme.
    temperature = A\(B*temperature)
    # temperature = A\temperature
    # The scaling of the temperature.
    potential_energy = discrete_quad(dis_V.*dis_density,
                            xAxis[1], xAxis[end])
    scaling = (energy - potential_energy)/
                ((1/alpha)*discrete_quad(temperature, xAxis[1], xAxis[end]))
    # Return the scaled temperature.
    temperature*scaling
end
"""
    evolveP(density::AbstractArray, evolveTime::Number, dt::Number,
            potentialTup::Tuple{Number, AbstractArray, Number},
            temperature::AbstractArray, xAxis::AbstractArray)
Evolve the probability density forward by an amount of time evolveTime using a
time step dt.
# Arguments:
* `density::AbstractArray`: The initial value for the probability density, must
be a vector of the same size of the xAxis.
* `evolveTime::Number`:  The amount of time to evolve the probability density
for in the dimensionless time unit.
* `dt::Number`:          The time step that we are using for the evolution.
* `potentialTup::Tuple{Number, AbstractArray, Number}`: A tuple containing
information on the potential, the first element is the potential to the left of
the left boundary, the second element is a vector containing the values of the
potential on the xAxis, the third element is the potential to the right of the
right boundary, the third element is the potential to the right of the right
boundary.
* `temperature::Number`: A vector containing the discretized temperature, must
be the same length as the xAxis.
* `xAxis::AbstractArray`: A vector describing the axis that we are discretizing
over in the dimensionless coordinates.
"""
function evolveP(density::AbstractArray, evolveTime::Number, dt::Number,
            potentialTup::Tuple{Number, AbstractArray, Number},
            temperature::AbstractArray, xAxis::AbstractArray)
    n_steps = round(Int, evolveTime/dt)
    for i = 1:n_steps
        density = stepP(density, dt, potentialTup, temperature, xAxis)
    end
    density
end
"""
    evolveT(temperature::AbstractArray, evolveTime::Number, dt::Number,
            potentialTup::Tuple{Number, AbstractArray, Number},
            density::AbstractArray, alpha::Number, beta::Number,
            energy::Number, xAxis::AbstractArray)
Evolve the dicrete temperature forward by an amount evolveTime using a time
step dt. This function keeps the probability density constant.
# Arguments:
* `temperature::AbstractArray`: A vector containing the discretized
temperature, must be the same length as the xAxis.
* `evolveTime::Number`: The amount of time to evolve the probability density for
 in the dimensionless time unit.
* `dt::Number`: The time step that we are using for the evolution.
* `potentialTup::::Tuple{Number, AbstractArray, Number}`: A tuple containing
information on the potential, the first element is the potential to the left of
the left boundary, the second element is a vector containing the values of the
potential on the xAxis, the third element is the potential to the right of the
right boundary, the third element is the potential to the right of the right
boundary.
* `density::AbstractArray`: The initial value for the probability density, must be a
             vector of the same size of the xAxis.
* `alpha::Number`: Dimensionless parameter that describes the coupling between
             the probability density and the temperature.
* `beta::Number`: Dimensionless parameter that describes how quickly the
temperature diffuses to a constant value.
* `energy::Number`: The dimensionless energy of the system.
* `xAxis::AbstractArray`: A vector describing the axis that we are discretizing
over in the dimensionless coordinates.
"""
function evolveT(temperature::AbstractArray, evolveTime::Number, dt::Number,
            potentialTup::Tuple{Number, AbstractArray, Number},
            density::AbstractArray, alpha::Number, beta::Number,
            energy::Number, xAxis::AbstractArray)
    n_steps = round(Int, evolveTime/dt)
    for i = 1:n_steps
        # Update temperature.
        temperature = stepT(temperature, dt, density, potentialTup, alpha,
                    beta, energy, xAxis)
    end
    temperature
end
"""
    evolve_system(density::AbstractArray, temperature::AbstractArray,
            evolveTime::Number, dt::Number,
            potentialTup::Tuple{Number, AbstractArray, Number}, alpha::Number,
            beta::Number, energy::Number,
            xAxis::AbstractArray)
Evolve the entire coupled system forward by an amount evolveTime using a time
step dt.
# Arguments:
density:     The initial value for the probability density, must be a
             vector of the same size of the xAxis.
temperature: A vector containing the discretized temperature, must be the
             same length as the xAxis.
evolveTime:  The amount of time to evolve the probability density for in
             the dimensionless time unit.
dt:          The time step that we are using for the evolution.
potentialTup:A tuple containing information on the potential, the first
             element is the potential to the left of
             the left boundary, the second element is a vector
             containing the values of the potential on the xAxis,
             the third element is the potential to the right of the
             right boundary, the third element is the potential to the
             right of the right boundary.
alpha:       Dimensionless parameter that describes the coupling between
             the probability density and the temperature.
beta:        Dimensionless parameter that describes how quickly the
             temperature diffuses to a constant value.
energy:      The dimensionless energy of the system.
xAxis:      A vector describing the axis that we are discretizing over in
             the dimensionless coordinates.
"""
function evolve_system(density::AbstractArray, temperature::AbstractArray,
            evolveTime::Number, dt::Number,
            potentialTup::Tuple{Number, AbstractArray, Number}, alpha::Number,
            beta::Number, energy::Number,
            xAxis::AbstractArray)
        n_steps = round(Int, evolveTime/dt)
    for i = 1:n_steps
        temperature = stepT(temperature, dt, density, potentialTup, alpha, beta,
                                energy, xAxis)
        density = stepP(density, dt, potentialTup, temperature, xAxis)
    end
    density, temperature
end
"""
    steady_state(density::AbstractArray, temperature::AbstractArray,
            dt::Number, potentialTup::Tuple{Number, AbstractArray, Number},
            alpha::Number, beta::Number, energy::Number,
            xAxis::AbstractArray, tol::Number)
Evolve the system forward until it is changing by less than the tolerance.
Returns the steady state probability density function and the temperature.
# Arguments:
* `density::AbstractArray`: The initial value for the probability density, must
be a vector of the same size of the xAxis.
* `temperature::AbstractArray`: A vector containing the discretized temperature,
 must be the same length as the xAxis.
* `dt::Number`: The time step that we are using for the evolution.
* `potentialTup::Tuple{Number, AbstractArray, Number}`: A tuple containing
information on the potential, the first element is the potential to the left of
the left boundary, the second element is a vector containing the values of the
potential on the xAxis, the third element is the potential to the right of the
right boundary, the third element is the potential to the right of the right
boundary.
* `alpha::Number`: Dimensionless parameter that describes the coupling between
the probability density and the temperature.
* `beta::Number`: Dimensionless parameter that describes how quickly the
temperature diffuses to a constant value.
* `energy::Number`: The dimensionless energy of the system.
* `xAxis::Number`: A vector describing the axis that we are discretizing over in
the dimensionless coordinates.
* `tol::Number`: The distance between vectors that we will tolerate before
saying that the system has reached the steady state.
"""
function steady_state(density::AbstractArray, temperature::AbstractArray,
            dt::Number, potentialTup::Tuple{Number, AbstractArray, Number},
            alpha::Number, beta::Number, energy::Number,
            xAxis::AbstractArray, tol::Number)
    density_new = density
    density_old = density + 2tol
    temperatureNew = temperature
    temperatureOld = temperature + 2tol

    system_energy = energyFun(dis_V, density, temperature, xAxis)
    # Step the system forward until both the probability density and the
    # temperature are changing by less than the tolerance.
    iters = 1
    while (norm(density_new - density_old) > tol ||
            norm(temperatureNew - temperatureOld) > tol)
        density_old = density_new
        temperatureOld = temperatureNew

        temperatureNew = stepT(temperature, dt, density, potentialTup,
                            alpha, beta, system_energy, xAxis)
        density_new = stepP(density, dt, potentialTup, temperature, xAxis)
        if iters > 2000
            return (norm(density_new - density_old),
                    norm(temperatureNew - temperatureOld))
        end
        iters += 1
    end
    density_new, temperatureNew
end

"""
    hopping_time(potentialTup::Tuple{Number, AbstractArray, Number},
                bump::Number, density::AbstractArray,
                temperature::AbstractArray, alpha::Number, beta::Number,
                xAxis::AbstractArray; dt=1e-4, tol=0.1)
Given an initial state of the system and a potential with a bump in it,
calculate how long it takes for the probability distribution to get over
the bump.
# Arguments:
* `potentialTup::Tuple{Number, AbstractArray, Number}`: A tuple containing the
potential.
* `bump::Number`: The location of the bump in the potetnial.
* `density::AbstractArray`: The initial density of the system.
* `temperature::AbstractArray`: The initial temperature of the system.
* `alpha::Number`: Dimensionless parameter that controls the coupling of the
system.
* `beta::Number`: Dimensionless parameter that determines how fast the
temperature gradients diffuse.
* `xAxis::AbstractArray`: The axis that we are working on.
* `dt=1e-4`: Size of one time step for simulation.
"""
function hopping_time(potentialTup::Tuple{Number, AbstractArray, Number},
                bump::Number, density::AbstractArray,
                temperature::AbstractArray, alpha::Number, beta::Number,
                xAxis::AbstractArray; dt=1e-4)
    bump_ind = indmin(abs(xAxis - bump))
    iters = 0
    energy = energyFun(potentialTup[2], density, temperature, alpha, xAxis)
    # If the mean of the probability distribution crosses the bump, then the
    # system has crossed the bump.
    while discrete_quad(density.*xAxis, xAxis[1], xAxis[end]) > bump
        density = stepP(density, dt, potentialTup, temperature, xAxis)
        temperature = stepT(temperature, dt, density, potentialTup, alpha,
                            beta, energy, xAxis)
        iters += 1
        if iters > 5000
            println("alpha = $alpha , beta = $beta DNF")
            return dt*iters
        end
    end
    ret = dt*iters
end
