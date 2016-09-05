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
Calculate the energy of the system.
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
* `dtpotential::AbstractArray`: A vector of the derivative of the potential,
contains points just outside the boundaries so it needs to be the length of
the xAxis plus 2.
* `temperature::AbstractArray`: A vector of the tempeature.
* xAxis: A vector containing the x coordinates of the points that we are doing
finite differencing on, the points must be equally spaced.
* `bndType::Symbol`: Specify the solution or the derivative at the boudary,
    * `:dirichlet`: The density remains constant at the boudaries.
    * `:periodic`: The density is the same at the left and right sides.
    * `:neumann`: The derivative is zero at the boundaries.
"""
function stepP(P::AbstractArray, dt::Number,
            dpotential::AbstractArray,
            temperature::AbstractArray, xAxis::AbstractArray;
            bndType::Symbol=:absorbing, normalization::Bool = true, returnMatrix::Bool = false)
    # Augment the discrete temperature and the discrete potential since we will
    # be evaluating them at points beyond the boundary.
    T0 = temperature[1]
    # The derivative of the temperature is zero at the boundaries.
    temperature = [temperature[1] ; temperature ; temperature[end]]

    n_points = length(xAxis)
    dx = (xAxis[end] - xAxis[1])/n_points  # Grid spacing.
    rr = dt/(dx^2)  # Useful constant for making the matrix
    # First we will make a matrix A that represents the forward Euler step.
    # Create the diagonals of the matrix, in order to linearize the problem, we
    # had to assume that T(x, t + dt) = T(x, t). We could improve on this by
    # using an explicit method to estimate T(x, t + dt).
    diag_minus1 = rr*(temperature[1:end-3] - 0.5dx*dpotential[1:end-3])
    diag0 = -rr*2*temperature[2:end-1]
    diag1 = rr*(temperature[4:end] + 0.5dx*dpotential[4:end])
    if bndType == :absorbing
        A = Tridiagonal(diag_minus1, diag0, diag1)
        II = Tridiagonal(zeros(diag1), ones(diag0), zeros(diag1))
        P = (II - 0.5A)\((II + 0.5A)*P)  # Since (1 - A/2)*P^{n+1} = (1 + A/2)*P^n.
    elseif bndType == :periodic
        A = spdiagm((diag_minus1, diag0, diag1), (-1, 0, 1))
        A[1, end] = diag_minus1[1]
        A[end, 1] = diag1[end]
        P = (speye(A) - 0.5A)\((speye(A) + 0.5A)*P)  # Since (1 - A/2)*P^{n+1} = (1 + A/2)*P^n.
    elseif bndType == :neumann
        # Derivative is zero at the boundaries.
        A = spdiagm((diag_minus1, diag0, diag1), (-1, 0, 1))
        B = speye(A) + 0.5A
        A = speye(A) - 0.5A
        A[1, 2] = -A[1, 1]
        A[end, end - 1] = -A[end, end]
        B[1, 2] = -B[1, 1]
        B[end, end - 1] = -B[end, end]
        P = A\(B*P)  # Since (1 - A/2)*P^{n+1} = (1 + A/2)*P^n.
    end
    if returnMatrix
        return A
    end
    if normalization
        # Return the normalized probability density.
        P = P/discrete_quad(P, xAxis[1], xAxis[end])
    end
    P
end

current(potential, density, temperature) = (density[2:end].*discrete_derivative(potential, xAxis)
                                          + discrete_derivative(density.*temperature, xAxis) )

""""
stepT(temperature::AbstractArray, dt::Number,
            density::AbstractArray, potential::AbstractArray,
            dpotential::AbstractArray, alpha::Number,
            beta::Number, energy::Number,
            xAxis::AbstractArray; bndType::Symbol=:neumann)
Evolve the temperature forward by an amount dt.
# Arguments:
* `temperature::AbstractArray`: Initial temperature, this is a vector of the
same length as the xAxis (see below).
* `dt::Number`: Amount of time to simulate forward by.
* `density::AbstractArray`: A vector of the probability density,
must be the same length as the xAxis (see below).
* `potential::AbstractArray`: A vector of the potential.
* `dtpotential::AbstractArray`: A vector of the derivative of the potential,
contains points just outside the boundaries so it needs to be the length of
the xAxis plus 2.
* `alpha::Number`: Dimensionless parameter that describes how much the heat
from the Brownian particle affects the temperature.
* `beta::Number`: Dimensionless parameter that describes how quickly the
temperature diffuses to a constant value.
* `energy::Number`: The initial energy of the system, must be in
dimensionless units, i.e. the actual energy divided by E0, where E0 is the
characteristic energy of the system.
* `xAxis::AbstractArray`: A vector containing the x coordinates of the points
that we are doing finite differencing on, the points must be equally spaced.
Returns the temperature evolved forward by an amount dt as well as the
updated energy of the system, assumes periodic boundary conditions.
* `bndType::Symbol=:neumann`: The type of boundary condition for the temperature
can include `:neumann`, `:dirichlet`, `:periodic`.
"""
function stepT(temperature::AbstractArray, dt::Number,
            density::AbstractArray, potential::AbstractArray,
            dpotential::AbstractArray, alpha::Number,
            beta::Number, energy::Number,
            xAxis::AbstractArray; bndType::Symbol=:neumann)
    # Augment the discrete density since we will be
    # evaluating them at points beyond the boundary.
    # The density is periodic.
    density = [density[end] ; density ; density[1]]

    n_points = length(xAxis)
    dx = (xAxis[end] - xAxis[1])/n_points
    rr = dt/(dx^2)
    # The inhomogeniety at the end of the equation.
    in_homo = -rr*alpha*density[2:end-1].*(dpotential[2:end-1].^2)*dx^2
    # The diagonals of the matrix.
    diag_minus1 = rr*(-0.5*alpha*density[1:end-3].*dpotential[2:end-2]*dx + beta)
    diag0 = -2*beta*rr*ones(xAxis)
    diag1 = rr*(0.5*alpha*density[4:end].*dpotential[3:end-1]*dx + beta)
    # Add the inhomogeneity to the diagonals.
    diag_minus1 += in_homo[2:end]
    diag0 += in_homo
    diag1 += in_homo[1:end-1]

    if bndType == :absorbing
        A = Tridiagonal(diag_minus1, diag0, diag1)
        II = Tridiagonal(zeros(diag1), ones(diag0), zeros(diag1))
        temperature = (II - 0.5A)\((II + 0.5A)*temperature)
    elseif bndType == :periodic
        A = spdiagm((diag_minus1, diag0, diag1), (-1, 0, 1))
        A[1, end] = diag_minus1[1]
        A[end, 1] = diag1[end]
        temperature = (speye(A) - 0.5A)\((speye(A) + 0.5A)*temperature)
    elseif bndType == :neumann
        # Derivative is zero at the boundaries.
        A = spdiagm((diag_minus1, diag0, diag1), (-1, 0, 1))
        B = speye(A) + 0.5A
        A = speye(A) - 0.5A
        A[1, 2] = -A[1, 1]
        A[end, end - 1] = -A[end, end]
        B[1, 2] = -B[1, 1]
        B[end, end - 1] = -B[end, end]
        temperature = A\(B*temperature)
    end
    # The scaling of the temperature.
    potential_energy = discrete_quad(potential.*density[2:end-1],
                            xAxis[1], xAxis[end])
    scaling = (energy - potential_energy)/
                ((1/alpha)*discrete_quad(temperature, xAxis[1], xAxis[end]))
    # Return the scaled temperature.
    temperature*scaling
end
"""
    evolveP(density::AbstractArray, evolveTime::Number, dt::Number,
            dpotential::AbstractArray, temperature::AbstractArray,
            xAxis::AbstractArray ; probBndType::Symbol = :periodic, tempBndType = :neumann)
Evolve the probability density forward by an amount of time evolveTime using a
time step dt.
# Arguments:
* `density::AbstractArray`: The initial value for the probability density, must
be a vector of the same size of the xAxis.
* `evolveTime::Number`:  The amount of time to evolve the probability density
for in the dimensionless time unit.
* `dt::Number`:          The time step that we are using for the evolution.
* `dtpotential::AbstractArray`: A vector of the derivative of the potential,
contains points just outside the boundaries so it needs to be the length of
the xAxis plus 2.
* `temperature::AbstractArray`: A vector containing the discretized temperature, must
be the same length as the xAxis.
* `xAxis::AbstractArray`: A vector describing the axis that we are discretizing
over in the dimensionless coordinates.
* `probBndType::Symbol = :periodic` The type of boundary condition on the probability
density.
"""
function evolveP(density::AbstractArray, evolveTime::Number, dt::Number,
            dpotential::AbstractArray, temperature::AbstractArray,
            xAxis::AbstractArray ; probBndType::Symbol = :periodic)
    n_steps = round(Int, evolveTime/dt)
    for i = 1:n_steps
        density = stepP(density, dt, dpotential, temperature, xAxis ; bndType = probBndType)
    end
    density
end
"""
    evolveT(temperature::AbstractArray, evolveTime::Number, dt::Number,
            potential::AbstractArray, dpotential::AbstractArray,
            density::AbstractArray, alpha::Number, beta::Number,
            energy::Number, xAxis::AbstractArray ; tempBndType::Symbol = :neumann)
Evolve the dicrete temperature forward by an amount evolveTime using a time
step dt. This function keeps the probability density constant.
# Arguments:
* `temperature::AbstractArray`: A vector containing the discretized
temperature, must be the same length as the xAxis.
* `evolveTime::Number`: The amount of time to evolve the probability density for
 in the dimensionless time unit.
* `dt::Number`: The time step that we are using for the evolution.
* `potential::AbstractArray`: A vector of the potential.
* `dtpotential::AbstractArray`: A vector of the derivative of the potential,
contains points just outside the boundaries so it needs to be the length of
the xAxis plus 2.
* `density::AbstractArray`: The initial value for the probability density, must be a
             vector of the same size of the xAxis.
* `alpha::Number`: Dimensionless parameter that describes the coupling between
             the probability density and the temperature.
* `beta::Number`: Dimensionless parameter that describes how quickly the
temperature diffuses to a constant value.
* `energy::Number`: The dimensionless energy of the system.
* `xAxis::AbstractArray`: A vector describing the axis that we are discretizing
over in the dimensionless coordinates.
* `tempBndType::Symbol = :neumann` The type of boundary condition on the temperature.
"""
function evolveT(temperature::AbstractArray, evolveTime::Number, dt::Number,
            potential::AbstractArray, dpotential::AbstractArray,
            density::AbstractArray, alpha::Number, beta::Number,
            energy::Number, xAxis::AbstractArray ; tempBndType::Symbol = :neumann)
    n_steps = round(Int, evolveTime/dt)
    for i = 1:n_steps
        # Update temperature.
        temperature = stepT(temperature, dt, density, potential, dpotential, alpha,
                beta, energy, xAxis ; bndType = tempBndType)
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
* `density::AbstractArray`: The initial value for the probability density, must
 be a vector of the same size of the xAxis.
* `temperature::AbstractArray`: A vector containing the discretized
temperature, must be the same length as the xAxis.
* `evolveTime::Number`: The amount of time to evolve the probability density
 for in the dimensionless time unit.
* `dt::Number`: The time step that we are using for the evolution.
* `potential::AbstractArray`: A vector of the potential.
* `dtpotential::AbstractArray`: A vector of the derivative of the potential,
contains points just outside the boundaries so it needs to be the length of
the xAxis plus 2.
* `alpha::Number`: Dimensionless parameter that describes the coupling between
 the probability density and the temperature.
* `beta::Number`: Dimensionless parameter that describes how quickly the
 temperature diffuses to a constant value.
* `energy::Number`: The dimensionless energy of the system.
* `xAxis::AbstractArray`: A vector describing the axis that we are discretizing
 over in the dimensionless coordinates.
* `tempBndType::Symbol = :neumann` The type of boundary condition on the temperature.
* `probBndType::Symbol = :poeriodic`: The type of boundary condition for the probability
distribution.
"""
function evolve_system(density::AbstractArray, temperature::AbstractArray,
            evolveTime::Number, dt::Number,
            potential::AbstractArray,
            dpotential::AbstractArray, alpha::Number,
            beta::Number, energy::Number,
            xAxis::AbstractArray ; tempBndType::Symbol = :neumann, probBndType::Symbol = :periodic)
        n_steps = round(Int, evolveTime/dt)
    for i = 1:n_steps
        temperature = stepT(temperature, dt, density, potential, dpotential, alpha, beta,
                        energy, xAxis ; bndType = tempBndType)
        density = stepP(density, dt, dpotential, temperature, xAxis ; bndType = probBndType)
    end
    density, temperature
end

"""
    hermite_coeff(points::AbstractArray, fvals::AbstractArray, dfvals::AbsrtactArray)
Use Hermite interpolation to fit the data, points is an array of the points in the x axis being used, fvals are the
function values at those points, dfvals are the derivatives at those points.
# Examples
Say that you have a function that you want to interpolate, you know the function and its first derivative at the
following points:

x     = [0.0, 1.0, 2.0]

f(x)  = [0.0, 8.0, 2.0]

f'(x) = [0.0, 0.0, 0.0]

You can find the coefficients of the polynomial that interpolates this data as follows.
```julia
julia> hermite_coeff([0.0, 1.0, 2.0], [0.0, 8.0, 2.0], [0.0, 0.0, 0.0])
6-element Array{Float64,1}:
   0.0
   0.0
  35.5
 -40.5
  14.5
  -1.5
```
"""
function hermite_coeff(points::AbstractArray, fvals::AbstractArray, dfvals::AbstractArray)
    polynomialOrder = length(fvals) + length(dfvals)
    A = Array(Float64, length(fvals) + length(dfvals), polynomialOrder)
    itr = 1
    for ii in 1:length(points)
        A[ii, :] = [points[itr]^n for n in 0:(polynomialOrder-1)]
        itr += 1
    end
    itr = 1
    for ii in (1:length(dfvals)) + length(fvals)
        A[ii, :] = [n > 0 ? n*points[itr]^(n-1) : 0.0 for n in 0:(polynomialOrder-1)]
        itr += 1
    end
    A\[fvals ; dfvals]
end
