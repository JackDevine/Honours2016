#Functions for finite differencing

function sparse_stepP(P, A, time_step)
    #Function that takes a probability ditribution P, a finite differencing #matrix A and a time_step and steps the distribution forward in time by an #amount equal to time_step. This function is one dimensional and assumes #that the temperature and the potential do not change with time.
    return (speye(size(A)...) - 0.5*Δt*A)\((speye(size(A)...) + 0.5*Δt*A)*P)
end

#This version does not use global variables
const k_B = 1.38e-3
#TODO should I evaluate the temperature and the potential functions at all
#points on the grid? This way I could pass both as an array and I could
#gaurantee that all points are only evaluated once.
function sparse_stepP(P::Array, time_step::Number, V::Function,
     T_fun::Function, n_periods::Integer, time::Number)
    #Accepts an initial density, a time step a ptoetial (as a function), the
    #temperature (as a function), the number of periods and the time that we are
    #at. The time is for the temperature and the potential, incase they are time
    #dependant.
    J = length(P)
    delta = 1.0/J
    diag1 = Float64[V(delta*(j+1), time) - V(delta*j, time) + k_B*T_fun(delta*j, time) for j in 1:J-1]

    diag0 = Float64[V(delta*(j-1), time) - V(delta*j, time) - 2k_B*T_fun(delta*j, time) for j in 1:J]

    diag_minus1 = Float64[k_B*T_fun(delta*j, time) for j in 2:J]

    A = spdiagm((diag0, diag1, diag_minus1), (0, 1, -1))
    return (speye(size(A)...) - 0.5*time_step*A)\((speye(size(A)...) + 0.5*time_step*A)*P)
end

function discrete_quad(vec::Array, start::Number, fin::Number)
    #Use the trapezoid rule to do quadrature on a vector vec
    #where the points in vec are evenly spaced from start to fin
    h = (fin - start)/length(vec)
    (h/2)*(vec[1] + vec[end] + 2*sum(vec[2:end-1]))
end

#Functions for sthochastic simulations

#An approximation of the dirac delta for calculating derivatives.
dirac_delta(x, sigma) = (1/(sigma*sqrt(pi)))*exp(-(x^2)/(2sigma^2))

function simulate_coupled_particle(position::Number, x_axis, Δt::Number,
    μ::Function, dis_dV::Array, dis_ddV::Array, dis_T::Array, kappa::Number,
    D::Number)
    #Simulate a particle that is coupled to the environment. Use this method for
    #a discretized version of the temperature where the temperature is
    #calculated at all points at each step, so that the discretized temperature
    #can be used for finite differencing. This means that the x_axis will be
    #restricted. The particle current at time t_n and position x is one if
    #(X(t_n) < x && X(t_n+1) > x) || ...

    Δ = 1/length(dis_T)
    Y = sqrt(2*k_B*Δt)*randn()
    position_index = indmin(abs(x_axis - position))

    new_position = position + μ(position)*Δt +
            sqrt(2*k_B*Δt)*randn()*sqrt(dis_T[position_index])
    new_position_index = indmin(abs(x_axis - new_position))

    epsilon = 0.0005
    #Should I make the current fuzzy? If it is, then that will affect the
    #derivative of the current
    current = Float64[(position < x < new_position) ? 1 : (new_position < x <
            position) ? -1 : 0 for x in x_axis]
    d_current =  Float64[sign(new_position - position)*(
            dirac_delta(x - position, epsilon) -
            dirac_delta(x - new_position, epsilon)) for x in x_axis]

    dis_T = dis_T - kappa*Δt*(d_current.*dis_dV + current.*dis_ddV)

    #Now apply diffusion to the temperature, since we are not given the matrix
    #for the finite differencing, we will have to generate it ourselves.
    diag_minus1 = Float64[D for j in 2:n_points]
    diag0 = Float64[-2D for j in 1:n_points]
    diag1 = Float64[D for j in 1:n_points-1]

    A = spdiagm((diag_minus1, diag0, diag1), (-1, 0, 1))
    A[1, end] = D
    A[end, 1] = D
    dis_T = (speye(size(A)...) - 0.5*Δt*A)\((speye(size(A)...) + 0.5*Δt*A)*dis_T)

    return new_position, dis_T
end
