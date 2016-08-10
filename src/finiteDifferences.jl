# Useful functions for doing finite differences on the system, all of these
# functions use the dimensionalized equations.

function discrete_quad(vec::AbstractArray, start::Number, fin::Number)
    # Use the trapezoid rule to do discrete quadrature on a vector vec where the
    # points in vec are evenly spaced from start to fin.
    h = (fin - start)/length(vec)
    (h/2)*(vec[1] + vec[end] + 2*sum(vec[2:end-1]))
end

function discrete_derivative(vec::AbstractArray, x_axis::AbstractArray)
    # Calculate the derivative of vec along the axis x_axis, note that the
    # vector returned is one element shorter than the vectors given.
    diff(vec)./diff(x_axis)
end

function energyFun(potential::AbstractArray, density::AbstractArray,
            temperature::AbstractArray, alpha::Number,
            x_axis::AbstractArray)
    # Given the potential, the probability density, the temperature of the
    # system and the x_axis that we are working on, calcualte the energy.
    left_bnd = x_axis[1]
    right_bnd = x_axis[end]
    (discrete_quad(potential.*density, left_bnd, right_bnd)
        + (1/alpha)*discrete_quad(temperature, left_bnd, right_bnd))
end

function stepP(P::AbstractArray, dt::Number,
            potentialTup::Tuple{Number, AbstractArray, Number},
            temperature::AbstractArray, x_axis::AbstractArray)
    #=
    Evolve the probability density P forward by an amount dt.
    Parameters:
    P:             Initial probability density, this is a vector of the same
                   length as the x_axis (see below).
    dt:            Amount of time to simulate forward by.
    potentialTup:  A tuple containing information on the potential, the first
                   element is the potential to the left of
                   the left boundary, the second element is a vector
                   containing the values of the potential on the x_axis,
                   the third element is the potential to the right of the
                   right boundary.
    temperature:   A vector of the tempeature.
    x_axis:        A vector containing the x coordinates of the points that we
                   are doing finite differencing on, the
                   points must be equally spaced.
    Returns the probability density evolved forward by an amount dt, assumes
    periodic boundary conditions.
    =#
    # Augment the discrete temperature and the discrete potential since we will
    # be evaluating them at points beyond the boundary.
    T0 = temperature[1]
    # The derivative of the temperature is zero at the boundaries.
    temperature = [temperature[1] ; temperature ; temperature[end]]

    dis_V0, dis_V, dis_V_end = potentialTup
    dis_V = [dis_V0 ; dis_V ; dis_V_end]

    n_points = length(x_axis)
    dx = (x_axis[end] - x_axis[1])/n_points  # Grid spacing
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

    A = spdiagm((diag_minus1, diag0, diag1), (-1, 0, 1))
    # Periodic boundary conditions.
    # A[1, end] = diag_minus1[1]
    # A[end, 1] = diag1[end]

    B = 2speye(size(A)...) - A  # B represents the backward Euler step.
    P = A\(B*P)  # Since A*P^{n+1} = B*P^n
    # Return the normalized probability density
    P = P/discrete_quad(P, x_axis[1], x_axis[end])
end

function stepT(temperature::AbstractArray, dt::Number, dis_density::AbstractArray,
            potentialTup::Tuple{Number, AbstractArray, Number}, alpha::Number,
            beta::Number, energy::Number,
            x_axis::AbstractArray)
    #=
    Evolve the temperature forward by an amount dt.
    Parameters:
    temperature:   Initial temperature, this is a vector of the same length as
                   the x_axis (see below).
    dt:            Amount of time to simulate forward by.
    dis_density:   A vector of the probability density, must be the same
                   length as the x_axis (see below).
    potentialTup:  A tuple containing information on the potential, the first
                   element is the potential to the left of
                   the left boundary, the second element is a vector
                   containing the values of the potential on the x_axis,
                   the third element is the potential to the right of the
                   right boundary, the third element is the potential to the
                   right of the right boundary.
    alpha:         Dimensionless parameter that describes how much the heat
                   from the motor affects the temperature.
    beta:          Dimensionless parameter that describes how quickly the
                   temperature diffuses to a constant value.
    init_energy:   The initial energy of the system, must be in dimensionless
                   units, i.e. the actual energy divided by E0, where E0 is
                   the potetnial energy gained by moving one period up the
                   potential (characteristic energy of the system).
    x_axis:        A vector containing the x coordinates of the points that we
                   are doing finite differencing on, the points must be
                   equally spaced.
    Returns the temperature evolved forward by an amount dt as well as the
    updated energy of the system, assumes periodic boundary conditions.
    =#
    # Augment the discrete density and the discrete potential since we will be
    # evaluating them at points beyond the boundary.
    dis_V0, dis_V, dis_V_end = potentialTup
    dis_V = [dis_V0 ; dis_V ; dis_V_end]
    # The density is periodic.
    dis_density = [dis_density[end] ; dis_density ; dis_density[1]]

    n_points = length(x_axis)
    dx = (x_axis[end] - x_axis[1])/n_points
    rr = dt/(2dx^2)
    # The inhomogeniety at the end of the equation.
    in_homo = rr*(alpha/4)*dis_density[2:end-1].*(dis_V[3:end] - dis_V[1:end-2])
    # The diagonals of the matrix.
    diag_minus1 = rr*((alpha/4)*dis_density[1:end-3]
                .*(dis_V[3:end-1] - dis_V[1:end-3]) - beta)
    diag0 = (2*beta*rr + 1)*ones(x_axis)
    diag1 = rr*(-(alpha/4)*dis_density[4:end]
            .*(dis_V[4:end] - dis_V[3:end-1]) - beta)
    # Add the inhomogeneity to the diagonals.
    diag_minus1 += in_homo[2:end]
    diag0 += in_homo
    diag1 += in_homo[1:end-1]
    A = spdiagm((diag_minus1, diag0, diag1), (-1, 0, 1))
    # Periodic boundary conditions.
    # A[1, end] = diag_minus1[1]
    # A[end, 1] = diag1[end]
    # Neumann boundary conditions (derivative is zero at the boundaries).
    temperature[1], temperature[end] = 0.0, 0.0
    A[1, 2] = -A[1, 1]
    A[end, end - 1] = -A[end, end]
    B = 2speye(size(A)...) - A
    # Update the temperature using the Crank Nicolson scheme.
    # temperature = A\(B*temperature)
    temperature = A\temperature
    # The scaling of the temperature.
    potential_energy = discrete_quad(dis_V.*dis_density,
                            x_axis[1], x_axis[end])
    scaling = (energy - potential_energy)/
                ((1/alpha)*discrete_quad(temperature, x_axis[1], x_axis[end]))
    # Return the scaled temperature.
    temperature*scaling
end

function evolveP(density::AbstractArray, evolveTime::Number, dt::Number,
            potentialTup::Tuple{Number, AbstractArray, Number},
            temperature::AbstractArray, x_axis::AbstractArray)
    #=
    evolveP will evolve the probability density forward by an amount of time
    evolveTime using a time step dt.
    Parameters:
    density:     The initial value for the probability density, must be a
                 vector of the same size of the x_axis.
    evolveTime:  The amount of time to evolve the probability density for in
                 the dimensionless time unit.
    dt:          The time step that we are using for the evolution.
    potentialTup:A tuple containing information on the potential, the first
                 element is the potential to the left of
                 the left boundary, the second element is a vector
                 containing the values of the potential on the x_axis,
                 the third element is the potential to the right of the
                 right boundary, the third element is the potential to the
                 right of the right boundary.
    temperature: A vector containing the discretized temperature, must be the
                 same length as the x_axis.
    x_axis:      A vector describing the axis that we are discretizing over in

    =#             the dimensionless coordinates.
    n_steps = round(Int, evolveTime/dt)
    for i = 1:n_steps
        density = stepP(density, dt, potentialTup, temperature, x_axis)
    end
    density
end

function evolveT(temperature::AbstractArray, evolveTime::Number, dt::Number,
            potentialTup::Tuple{Number, AbstractArray, Number},
            density::AbstractArray, alpha::Number, beta::Number,
            energy::Number, x_axis::AbstractArray)
    #=
    evolveT will evolve the dicrete temperature forward by an amount
    evolveTime using a time step dt. This function keeps the probability
    density constant.
    Parameters:
    temperature:    A vector containing the discretized temperature, must be the
                 same length as the x_axis.
    evolveTime: The amount of time to evolve the probability density for in
                 the dimensionless time unit.
    dt:          The time step that we are using for the evolution.
    potentialTup:A tuple containing information on the potential, the first
                 element is the potential to the left of
                 the left boundary, the second element is a vector
                 containing the values of the potential on the x_axis,
                 the third element is the potential to the right of the
                 right boundary, the third element is the potential to the
                 right of the right boundary.
    density:     The initial value for the probability density, must be a
                 vector of the same size of the x_axis.
    alpha:       Dimensionless parameter that describes the coupling between
                 the probability density and the temperature.
    beta:        Dimensionless parameter that describes how quickly the
                 temperature diffuses to a constant value.
    energy:      The dimensionless energy of the system.
    x_axis:      A vector describing the axis that we are discretizing over in
                 the dimensionless coordinates.
    =#
    n_steps = round(Int, evolveTime/dt)
    for i = 1:n_steps
        # Update temperature.
        temperature = stepT(temperature, dt, density, potentialTup, alpha,
                    beta, energy, x_axis)
    end
    temperature
end

function evolve_system(density::AbstractArray, temperature::AbstractArray,
            evolveTime::Number, dt::Number,
            potentialTup::Tuple{Number, AbstractArray, Number}, alpha::Number,
            beta::Number, energy::Number,
            x_axis::AbstractArray)
    #=
    evolve_system will evolve the entire coupled system forward by an amount
    evolveTime using a time step dt.
    Parameters:
    density:     The initial value for the probability density, must be a
                 vector of the same size of the x_axis.
    temperature: A vector containing the discretized temperature, must be the
                 same length as the x_axis.
    evolveTime:  The amount of time to evolve the probability density for in
                 the dimensionless time unit.
    dt:          The time step that we are using for the evolution.
    potentialTup:A tuple containing information on the potential, the first
                 element is the potential to the left of
                 the left boundary, the second element is a vector
                 containing the values of the potential on the x_axis,
                 the third element is the potential to the right of the
                 right boundary, the third element is the potential to the
                 right of the right boundary.
    alpha:       Dimensionless parameter that describes the coupling between
                 the probability density and the temperature.
    beta:        Dimensionless parameter that describes how quickly the
                 temperature diffuses to a constant value.
    energy:      The dimensionless energy of the system.
    x_axis:      A vector describing the axis that we are discretizing over in
                 the dimensionless coordinates.
    =#
        n_steps = round(Int, evolveTime/dt)
    for i = 1:n_steps
        temperature = stepT(temperature, dt, density, potentialTup, alpha, beta,
                                energy, x_axis)
        density = stepP(density, dt, potentialTup, temperature, x_axis)
    end
    density, temperature
end

function steady_state(density::AbstractArray, temperature::AbstractArray,
            dt::Number, potentialTup::Tuple{Number, AbstractArray, Number},
            alpha::Number, beta::Number, energy::Number,
            x_axis::AbstractArray, tol::Number)
    #=
    Evolve the system forward until it is changing by less than the tolerance.
    Returns the steady state probability density function and the temperature.
    Parameters:
    density:     The initial value for the probability density, must be a
                 vector of the same size of the x_axis.
    temperature: A vector containing the discretized temperature, must be the
                 same length as the x_axis.
    dt:          The time step that we are using for the evolution.
    potentialTup:A tuple containing information on the potential, the first
                 element is the potential to the left of
                 the left boundary, the second element is a vector
                 containing the values of the potential on the x_axis,
                 the third element is the potential to the right of the
                 right boundary, the third element is the potential to the
                 right of the right boundary.
    alpha:       Dimensionless parameter that describes the coupling between
                 the probability density and the temperature.
    beta:        Dimensionless parameter that describes how quickly the
                 temperature diffuses to a constant value.
    energy:      The dimensionless energy of the system.
    x_axis:      A vector describing the axis that we are discretizing over in
                 the dimensionless coordinates.
    tol:         The distance between vectors that we will tolerate before
                 saying that the system has reached the steady state.
    =#
    density_new = density
    density_old = density + 2tol
    temperatureNew = temperature
    temperatureOld = temperature + 2tol

    system_energy = energyFun(dis_V, density, temperature, x_axis)
    # Step the system forward until both the probability density and the
    # temperature are changing by less than the tolerance.
    iters = 1
    while (norm(density_new - density_old) > tol ||
            norm(temperatureNew - temperatureOld) > tol)
        density_old = density_new
        temperatureOld = temperatureNew

        temperatureNew = stepT(temperature, dt, density, potentialTup,
                            alpha, beta, system_energy, x_axis)
        density_new = stepP(density, dt, potentialTup, temperature, x_axis)
        if iters > 2000
            return (norm(density_new - density_old),
                    norm(temperatureNew - temperatureOld))
        end
        iters += 1
    end
    density_new, temperatureNew
end
