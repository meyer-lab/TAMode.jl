using TAMode
using Plots
using OrdinaryDiffEq
using DiffEqDevTools
using IterativeSolvers

params = TAMode.param(ones(15) * 0.5)
solInit = TAMode.getAutocrine(params)
params.gasCur = 100.0
prob = ODEProblem(TAMode.TAMreact, solInit, (0.0, 1000.0), params)

setups = [
    Dict(:alg => Rodas4()),
    Dict(:alg => Rodas42()),
    Dict(:alg => Rodas4P()),
    Dict(:alg => Ros4LStab()),
    Dict(:alg => Rodas5()),
    Dict(:alg => Rosenbrock23()), # slow
    Dict(:alg => TRBDF2(autodiff = false, linsolve = LinSolveGMRES())), # slow
    Dict(:alg => ABDF2(autodiff = false, linsolve = LinSolveGMRES())), # slow
    # Dict(:alg => Exprb43()), # slow
    # Dict(:alg => Exprb32()), # slow
]

test_sol = solve(prob, AutoTsit5(Rodas5()), reltol = 1.0e-16, abstol = 1.0e-16)

tols = 1.0 ./ 10.0 .^ LinRange(9, 12, 4)

wp = WorkPrecisionSet(prob, tols, tols, setups; save_everystep = false, appxsol = test_sol, maxiters = Int(1e6))

plot(wp)
