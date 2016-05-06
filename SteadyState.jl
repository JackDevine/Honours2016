# Deal with the steady state by using the analytical solution given by
# Gardiner

# Fisrt we create the objective function that needs to be optimized.

using JuMP
using NLopt
using FastAnonymous

const T0 = 300
kappa = 0.3
D = 0.025
f = 1.0
const L = 1.0
const k_B = 0.01
const v_0 = 1.0

v0(x) =  v_0*sin((2*pi/L)*x)
dv0(x) = v_0*(2*pi/L)*cos((2*pi/L)*x)
VV(x) = f*x + v0(x)
dVV(x) = f + (2*pi/L)*v_0*cos((2*pi/L)*x)
function temperature(xx, J_s)
    # Return the steady state temperature for a given postion xx and a trial
    # steady state current J_s
    T0 + (kappa*J_s/D)*(v0(xx) - v0(0))
end

function exp_int(xx, J_s)
    # A helper function for the objective function to calculate the exponential
    # integrals
    quadgk(yy -> dVV(yy)/(k_B*temperature(yy, J_s)), 0, xx)[1]
end

function psi(xx, J_s)
    # Helper function for the density
    exp(-quadgk(yy -> dVV(yy)/(2k_B*temperature(yy, J_s)), 0, xx)[1])
end

function density(xx, J_s, P0)
    # Given the current J_s and the density at the boundaires P0, calculate the
    # density at xx
    P0*( ((temperature(L, J_s)/psi(L, J_s))*quadgk(yy -> 1/psi(yy, J_s), 0,
        xx)[1] +
    (temperature(0, J_s)/psi(0, J_s))*quadgk(yy -> 1/psi(yy, J_s), xx, L)[1] ) /
    (temperature(xx, J_s)/psi(xx, J_s))*quadgk(yy -> 1/psi(yy, J_s), 0, L)[1] )
end

function objective_fun(J_s, P0)
    # An abjective function for J_s, this should be 0 for the correct J_s
    temp1 = 2k_B*((temperature(L, J_s)/psi(L, J_s) - (temperature(0, J_s)
          /psi(0,J_s)) ))
    temp2 = (quadgk(yy -> 1/psi(yy, J_s), 0, L)[1])^-1
    (J_s - temp1*P0*temp2)^2
end

registerNLFunction(:objective_fun, 2, objective_fun, autodiff=true)

area(J_s, P0) = quadgk(yy -> density(yy, J_s, P0), 0, L )[1] # The integral of
                                                    #the density over one period
area_tol = 1e-10 # The amount that we will allow the area to differ from 1 by.
registerNLFunction(:area, 2, area, autodiff=true)

function steady_params(T0, kappa, D, f)
    # Given the parameters for the environment estimate the steady state current
    # and the value of the density at the boundaries.
    m = Model(solver=NLoptSolver(algorithm=:LD_SLSQP))
#     m = Model()

    @defVar(m, current)
    @defVar(m, 0 <= P0 <= 1)

    @setNLObjective(m, Min, objective_fun(current, P0) )
    @addNLConstraint(m, area(current, P0) == 1 ) # The density is normalized

    setValue(current, 1.0)
    setValue(P0, 0.1)

    status = solve(m)
    getValue(current), getValue(P0)
end
