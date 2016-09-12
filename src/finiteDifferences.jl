module FiniteDifferences
export stepP, stepT, energyFun, hermite_coeff, discrete_quad,
        discrete_derivative, System, measure_kramers, kramers_rate,
        kramers_rate_analytical, @constant, evolveP, evolveT, discrete_quad,
        discrete_derivative


"""
    @constant(varname, varvalue)
A macro for making the data in varvalue available to all of the workers
through varname.
"""
macro constant(varname, varvalue)
  tmp = eval(varvalue)
  quote
    for i in procs()
            @spawnat i global const $varname = $tmp
    end
  end
end
# Useful functions for doing finite differences on the system, all of these
# functions use the dimensionalized equations.
"""
    System
Contains all of the data for the system in arrays.
# Fields
* `potential::AbstractArray`
* `dpotential::AbstractArray`
* `density::AbstractArray`
* `temperature::AbstractArray`
* `xAxis::AbstractArray`
* `energy::Number`
"""
type System
    potential::AbstractArray
    dpotential::AbstractArray
    density::AbstractArray
    temperature::AbstractArray
    xAxis::AbstractArray
    energy::Number
end

"""
Initialize an empty system.
"""
System() = System([], [], [], [], [], 0)

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

function energyFun(system::System, alpha::Number)
        energyFun(system.potential, system.density, system.temperature,
                alpha, system.xAxis)
end
"""
    stepP(P::AbstractArray, dt::Number,
                dpotential::AbstractArray,
                temperature::AbstractArray, xAxis::AbstractArray;
                bndType::Symbol=:absorbing, normalization::Bool = true,
                returnMatrix::Bool = false)
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
* `normalization::Bool`: Whether or not to normalize the solution after
doing a finite differences step forward.
* `returnMatrix::Bool`: If true then the function will return the finite
differences matrix.
"""
function stepP(P::AbstractArray, dt::Number,
            dpotential::AbstractArray,
            temperature::AbstractArray, xAxis::AbstractArray;
            bndType::Symbol=:absorbing, normalization::Bool = true,
            returnMatrix::Bool = false)
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
    diag_minus1 = rr*(-0.5dpotential[1:end-3]*dx +
            0.25*(-temperature[3:end-1] + temperature[1:end-3]
            + 4temperature[2:end-2]))
    diag0 = -2rr*temperature[2:end-1]
    diag1 = rr*(0.5*dpotential[4:end]*dx + 0.25*(temperature[4:end]
                - temperature[2:end-2] + 4temperature[3:end-1]))
    if bndType == :absorbing
        # The density is zero at both boundaries.
        A = Tridiagonal(diag_minus1, diag0, diag1)
        II = Tridiagonal(zeros(diag1), ones(diag0), zeros(diag1))
        P = (II - 0.5A)\((II + 0.5A)*P)
    elseif bndType == :periodic
        # The density is the same at both boundaries.
        A = spdiagm((diag_minus1, diag0, diag1), (-1, 0, 1))
        A[1, end] = diag_minus1[1]
        A[end, 1] = diag1[end]
        P = (speye(A) - 0.5A)\((speye(A) + 0.5A)*P)
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

"""
    stepP(system::System, dt::Number;
            bndType::Symbol = :absorbing, normalization::Bool = true,
            returnMatrix::Bool = false)
"""
function stepP(system::System, dt::Number;
            bndType::Symbol = :absorbing, normalization::Bool = true,
            returnMatrix::Bool = false)
    # Return the value from the more terse method of stepP.
    if returnMatrix
        return stepP(system.density, dt, system.dpotential, system.temperature,
                        system.xAxis; bndType = bndType,
                        normalization = normalization, returnMatrix = true)
    end
    stepP(system.density, dt, system.dpotential, system.temperature,
            system.xAxis; bndType = bndType,
            normalization = normalization, returnMatrix = false)
end

"""
    current(potential, density, temperature)
Calculate the current of the system, where the current is the value inside
the derivative on the RHS of the Smoluchowski equation.
"""
function current(potential, density, temperature)
    (density[2:end].*discrete_derivative(potential, xAxis)
    + discrete_derivative(density.*temperature, xAxis) )
end
"""
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
* `returnMatrix::Bool`: If true then the function will return the finite
differences matrix.
"""
function stepT(temperature::AbstractArray, dt::Number,
            density::AbstractArray, potential::AbstractArray,
            dpotential::AbstractArray, alpha::Number,
            beta::Number, energy::Number,
            xAxis::AbstractArray; bndType::Symbol=:neumann,
            returnMatrix::Bool = false)
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
    diag_minus1 = rr*(-0.5*alpha*density[1:end-3].*dpotential[2:end-2]*dx
                        + beta)
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
        # The value at the left boundary is the same as the value at the right
        # boundary.
        A = spdiagm((diag_minus1, diag0, diag1), (-1, 0, 1))
        A[1, end] = diag_minus1[1]
        A[end, 1] = diag1[end]
        temperature = (speye(A) - 0.5A)\((speye(A) + 0.5A)*temperature)
    elseif bndType == :dirichlet  # The value at the boundary remains fixed.
        # Calculate the heat due to the motor at the boundaries, in the case
        # where we apply Dirichlet boundary conditions, the energy is no longer
        # conserved, so we have to subtract the energy lost as heat for each
        # step.
        heatLeft = (density[2]*dpotential[2] + temperature[1]*(density[3]
                                - density[1])/(2dx))*dpotential[2]
        heatRight = (density[end-1]*dpotential[end-1]
                          + temperature[end]*(density[end] -
                          density[end-2])/(2dx))*dpotential[end-1]
        heat = heatRight - heatLeft
        # Calculate the heat due to temperature gradients at the boundaries.
        heatLeft = (temperature[2] - temperature[1])/(2dx)
        heatRight = (temperature[end] - temperature[end-1])/(2dx)
        heat += heatRight - heatLeft
        energy += heat
        # The value at the boundary remains fixed.
        II = Tridiagonal(zeros(diag1), ones(diag0), zeros(diag1))
        diag0[1] = 0.0
        diag0[end] = 0.0
        diag1[1] = 0.0
        diag_minus1[end] = 0.0
        A = Tridiagonal(diag_minus1, diag0, diag1)
        diag0[1] = 0.0
        diag0[end] = 0.0
        B = Tridiagonal(diag_minus1, diag0, diag1)
        temperatureLeft, temperatureRight = temperature[1], temperature[end]
        temperature = (II - 0.5A)\((II + 0.5B)*temperature)
        # The scaling of the temperature.
        potential_energy = discrete_quad(potential.*density[2:end-1],
                            xAxis[1], xAxis[end])
        scaling = (energy - potential_energy)/
                ((1/alpha)*discrete_quad(temperature, xAxis[1], xAxis[end]))
        temperature = temperature*scaling
        temperature[1], temperature[end] = temperatureLeft, temperatureRight
        return energy, temperature
    elseif bndType == :neumann
        # Derivative is zero at the boundaries.
        II = Tridiagonal(zeros(diag1), ones(diag0), zeros(diag1))
        diag1[1] = -diag0[1]
        diag_minus1[end] = -diag0[end]
        A = Tridiagonal(diag_minus1, diag0, diag1)
        temperature = (II - 0.5A)\((II + 0.5A)*temperature)
    end
    if returnMatrix
        return A
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
    stepT(system::System, alpha::Number, beta::Number, dt::Number;
            bndType::Symbol = :absorbing, normalization::Bool = true,
            returnMatrix::Bool = false)
"""
function stepT(system::System, alpha::Number, beta::Number, dt::Number;
            bndType::Symbol = :absorbing, returnMatrix::Bool = false)
    # Return the value from the more terse method of stepT.
    if returnMatrix
        return stepT(system.temperature, dt, system.density, system.potential,
                system.dpotential, alpha, beta, system.energy, system.xAxis;
                bndType = bndType, returnMatrix = true)
    end
    stepT(system.temperature, dt, system.density, system.potential,
            system.dpotential, alpha, beta, system.energy, system.xAxis;
            bndType = bndType, returnMatrix = false)
end

"""
    evolveP(density::AbstractArray, evolveTime::Number, dt::Number,
            dpotential::AbstractArray, temperature::AbstractArray,
            xAxis::AbstractArray ; probBndType::Symbol = :periodic,
            tempBndType = :neumann)
Evolve the probability density forward by an amount of time evolveTime using a
time step dt.
# Arguments:
* `density::AbstractArray`: The initial value for the probability density, must
be a vector of the same size of the xAxis.
* `evolveTime::Number`: The amount of time to evolve the probability density
for in the dimensionless time unit.
* `dt::Number`: The time step that we are using for the evolution.
* `dtpotential::AbstractArray`: A vector of the derivative of the potential,
contains points just outside the boundaries so it needs to be the length of
the xAxis plus 2.
* `temperature::AbstractArray`: A vector containing the discretized
temperature, must be the same length as the xAxis.
* `xAxis::AbstractArray`: A vector describing the axis that we are discretizing
over in the dimensionless coordinates.
* `probBndType::Symbol = :periodic` The type of boundary condition on the
probability density.
"""
function evolveP(density::AbstractArray, evolveTime::Number, dt::Number,
            dpotential::AbstractArray, temperature::AbstractArray,
            xAxis::AbstractArray ; probBndType::Symbol = :periodic)
    n_steps = round(Int, evolveTime/dt)
    for i = 1:n_steps
        density = stepP(density, dt, dpotential, temperature, xAxis;
                        bndType = probBndType)
    end
    density
end
"""
    evolveT(temperature::AbstractArray, evolveTime::Number, dt::Number,
            potential::AbstractArray, dpotential::AbstractArray,
            density::AbstractArray, alpha::Number, beta::Number,
            energy::Number, xAxis::AbstractArray ;
            tempBndType::Symbol = :neumann)
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
* `density::AbstractArray`: The initial value for the probability density, must
be a vector of the same size of the xAxis.
* `alpha::Number`: Dimensionless parameter that describes the coupling between
             the probability density and the temperature.
* `beta::Number`: Dimensionless parameter that describes how quickly the
temperature diffuses to a constant value.
* `energy::Number`: The dimensionless energy of the system.
* `xAxis::AbstractArray`: A vector describing the axis that we are discretizing
over in the dimensionless coordinates.
* `tempBndType::Symbol = :neumann` The type of boundary condition on the
temperature.
"""
function evolveT(temperature::AbstractArray, evolveTime::Number, dt::Number,
            potential::AbstractArray, dpotential::AbstractArray,
            density::AbstractArray, alpha::Number, beta::Number,
            energy::Number, xAxis::AbstractArray;
            tempBndType::Symbol = :neumann)
    n_steps = round(Int, evolveTime/dt)
    for i = 1:n_steps
        # Update temperature.
        temperature = stepT(temperature, dt, density, potential, dpotential,
                alpha, beta, energy, xAxis ; bndType = tempBndType)
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
* `tempBndType::Symbol = :neumann` The type of boundary condition on the
temperature.
* `probBndType::Symbol = :poeriodic`: The type of boundary condition for the
probability
distribution.
"""
function evolve_system(density::AbstractArray, temperature::AbstractArray,
            evolveTime::Number, dt::Number,
            potential::AbstractArray,
            dpotential::AbstractArray, alpha::Number,
            beta::Number, energy::Number,
            xAxis::AbstractArray ; tempBndType::Symbol = :neumann,
            probBndType::Symbol = :periodic)
        n_steps = round(Int, evolveTime/dt)
    for i = 1:n_steps
        temperature = stepT(temperature, dt, density, potential, dpotential,
                                alpha, beta,
                        energy, xAxis ; bndType = tempBndType)
        density = stepP(density, dt, dpotential, temperature, xAxis ; bndType =
        probBndType)
    end
    density, temperature
end

"""
    hermite_coeff(points::AbstractArray, fvals::AbstractArray,
         dfvals::AbsrtactArray)
Use Hermite interpolation to fit the data, points is an array of the points in
the x axis being used, fvals are the
function values at those points, dfvals are the derivatives at those points.
# Examples
Say that you have a function that you want to interpolate, you know the
function and its first derivative at the
following points:

x     = [0.0, 1.0, 2.0]

f(x)  = [0.0, 8.0, 2.0]

f'(x) = [0.0, 0.0, 0.0]

You can find the coefficients of the polynomial that interpolates this data as
 follows.
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
function hermite_coeff(points::AbstractArray, fvals::AbstractArray,
     dfvals::AbstractArray)
    polynomialOrder = length(fvals) + length(dfvals)
    A = Array(Float64, length(fvals) + length(dfvals), polynomialOrder)
    itr = 1
    for ii in 1:length(points)
        A[ii, :] = [points[itr]^n for n in 0:(polynomialOrder-1)]
        itr += 1
    end
    itr = 1
    for ii in (1:length(dfvals)) + length(fvals)
        A[ii, :] = [n > 0 ? n*points[itr]^(n-1) : 0.0 for n in
         0:(polynomialOrder-1)]
        itr += 1
    end
    A\[fvals ; dfvals]
end

"""
    kramers_rate(pRight::AbstractArray, timeVec::AbstractArray)
Given the probability of the particle being in the top well, calculate the
Kramers rate for the system by measuring the mean of the log of the slope of
the probability of being in the right well.
# Examples
```julia
julia> timeVec = linspace(0, 2, 1000)
       pRight = Float64[exp(-2*t) for t in timeVec]
       kramers_rate(pRight, timeVec)  # 2.0
```
"""
function kramers_rate(pRight::AbstractArray, timeVec::AbstractArray)
    slope = discrete_derivative(log(pRight), timeVec)
    mean(abs(slope[end]))
end


"""
    measure_kramers(wellPositions::AbstractArray, system::System,
                        alpha::Number, beta::Number, dt::Number,
                        nSteps::Integer)
Given the system in an intial state as well as the values of `alpha` and
`beta`, calculate the krammers rate by evolving the system forward and running
`kramers_rate` on `pRight`, where `pRight` is the proability of being in the
right well.
"""
function measure_kramers(wellPositions::AbstractArray, system::System,
                        alpha::Number, beta::Number, dt::Number,
                        nSteps::Integer; probBndType::Symbol = :absorbing,
                        tempBndType::Symbol = :neumann)
    nPoints = length(system.density)
    energy = energyFun(system, alpha)
    systemLocal = System(system.potential, system.dpotential, system.density,
                            system.temperature, system.xAxis, energy)

    hIndex = round(Int, nPoints/2)
    pRight = Array(Float64, nSteps)
    pRight[1] = discrete_quad(system.density[hIndex:end],
                            system.xAxis[hIndex], system.xAxis[end])

    if tempBndType == :dirichlet
        for i in 2:nSteps
            systemLocal.density = stepP(systemLocal, dt; bndType = probBndType)
            systemLocal.energy, systemLocal.temperature = stepT(systemLocal,
                                    alpha, beta, dt; bndType = :dirichlet)
            pRight[i] = discrete_quad(systemLocal.density[hIndex:end],
                            systemLocal.xAxis[hIndex], systemLocal.xAxis[end])
        end
    end
    for i in 2:nSteps
        systemLocal.density = stepP(systemLocal, dt; bndType = probBndType)
        systemLocal.temperature = stepT(systemLocal, alpha,
                                                beta, dt; bndType = tempBndType)
        pRight[i] = discrete_quad(systemLocal.density[hIndex:end],
                        systemLocal.xAxis[hIndex], systemLocal.xAxis[end])
    end
    kramers_rate(pRight, (1:nSteps)*dt)
end

"""
    measure_kramers(wellPositions::AbstractArray, system::System,
                        dt::Number, nSteps::Integer)
Given the system in an intial state, calculate the krammers rate by evolving
the uncoupled system forward and running `kramers_rate` on `pRight`, where
`pRight` is the proability of being in the right well.
"""
function measure_kramers(wellPositions::AbstractArray, system::System,
                            dt::Number, nSteps::Integer;
                            probBndType::Symbol = :absorbing)
    nPoints = length(system.density)
    systemLocal = System(system.potential, system.dpotential, system.density,
                            system.temperature, system.xAxis, 0.0)

    hIndex = round(Int, nPoints/2)
    pRight = Array(Float64, nSteps)
        pRight[1] = discrete_quad(system.density[hIndex:end],
                            system.xAxis[hIndex], system.xAxis[end])

    for i in 2:nSteps
        systemLocal.density = stepP(systemLocal, dt; bndType = probBndType)
        pRight[i] = discrete_quad(systemLocal.density[hIndex:end],
                        systemLocal.xAxis[hIndex], systemLocal.xAxis[end])
    end
    kramers_rate(pRight, (1:nSteps)*dt)
end


"""
    measure_kramers(wellPositions::AbstractArray, coeff::AbstractArray,
                        initDensity::AbstractArray, temperature::AbstractArray,
                        xAxis::AbstractArray,
                        dt::Number, nSteps::Integer;
                        probBndType = :absorbing)
Measure the Kramers rate using the coefficients in `coeff` to create the
potential.
"""
function measure_kramers(wellPositions::AbstractArray, coeff::AbstractArray,
                        initDensity::AbstractArray, temperature::AbstractArray,
                        xAxis::AbstractArray,
                        dt::Number, nSteps::Integer;
                        probBndType = :absorbing)
    dx = (xAxis[end] - xAxis[1])/length(xAxis)
    function V(x, coeff)
        acc = 0
        for ii = 1:length(coeff)
            acc += coeff[ii]*x^(ii-1)
        end
        acc
    end
    V(x) = V(x, coeff)
    dV(x) = ForwardDiff.derivative(V, x)

    potential = Float64[V(x) for x in xAxis]
    dpotential = [dV(xAxis[1] - dx) ; Float64[dV(x) for x in xAxis] ;
                                        dV(xAxis[end] + dx)]

    nPoints = length(initDensity)
    energy = 0.0
    system = System(potential, dpotential, initDensity,
                            temperature, xAxis, energy)

    hIndex = round(Int, nPoints/2)
    pRight = Array(Float64, nSteps)
        pRight[1] = discrete_quad(system.density[hIndex:end],
                            system.xAxis[hIndex], system.xAxis[end])

    for i in 2:nSteps
        system.density = stepP(system, dt; bndType = probBndType)
        pRight[i] = discrete_quad(system.density[hIndex:end],
                        system.xAxis[hIndex], system.xAxis[end])
    end
    kramers_rate(pRight, (1:nSteps)*dt)
end

"""
    kramers_rate_analytical(coeff::AbstractArray, wellPositions::AbstractArray,
                                temperature::Number)
Given the coefficients of the polynomial representing the potential, calculate
the kramers rate from the upper well to the lower well.
# Examples
```julia
    wellPositions = [-2.0, 0.0, 2.0]
    coeff = hermite_coeff(wellPositions, [0.0, 7.0, 2.0], [0.0, 0.0, 0.0])
    kramers_rate_analytical(coeff, wellPositions, 1.0)  # 0.0061340150666157394
```
"""
function kramers_rate_analytical(coeff::AbstractArray,
                        wellPositions::AbstractArray, temperature::Number)
    potentialA = sum([coeff[i+1]*wellPositions[3]^i
                                for i in 0:(length(coeff) - 1)])
    potentialB = sum([coeff[i+1]*wellPositions[2]^i
                                for i in 0:(length(coeff) - 1)])
    potentialC = sum([coeff[i+1]*wellPositions[1]^i
                                for i in 0:(length(coeff) - 1)])
    ddpotentialA = sum([i*(i - 1)*coeff[i+1]*wellPositions[3]^(i - 2)
                                for i in 2:(length(coeff) - 1)])
    ddpotentialB = sum([i*(i - 1)*coeff[i+1]*wellPositions[2]^(i - 2)
                                for i in 2:(length(coeff) - 1)])
    ddpotentialC = sum([i*(i - 1)*coeff[i+1]*wellPositions[1]^(i - 2)
                                for i in 2:(length(coeff) - 1)])
    Eforwards = potentialB - potentialA
    Ebackwards = potentialB - potentialC
    (
    ((sqrt(abs(ddpotentialA*ddpotentialB)))/(2pi))*exp(-Eforwards/temperature)
     -
    ((sqrt(abs(ddpotentialC*ddpotentialB)))/(2pi))*exp(-Ebackwards/temperature)
    )
end

end  # module
