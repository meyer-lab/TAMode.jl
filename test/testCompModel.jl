
CompParams = ones(12) * 0.5


function noReact!(du, u)
    return TAMode.TAMreact(du, u, TAMode.compParamm(CompParams), 0.0, reaction = false)
end


@testset "Derivative is zero with no reaction and flat concentration." begin
    u = ones(27 * TAMode.compSize)
    du = ones(27 * TAMode.compSize)

    noReact!(du, u)

    @test_broken all(aboutZero.(du))
end


@testset "Look at Jacobian of function." begin
    u = zeros(27 * TAMode.compSize)
    du = zeros(27 * TAMode.compSize)

    JJ = ForwardDiff.jacobian(noReact!, du, u)
    JJ2 = ForwardDiff.jacobian(noReact!, du, u)

    @test all(JJ .=== JJ2)
    @test all(triu(JJ, 28) .== 0.0)
    @test all(tril(JJ, -28) .== 0.0)

    for ii = 2:26
        @test all(diag(JJ, ii) .== 0.0)
        @test all(diag(JJ, -ii) .== 0.0)
    end
end


# @testset "In the absence of a reaction we just see diffusion." begin
#     solInit = zeros(2700)
#     solInit[1] = 1000.0
# 
#     pp = TAMode.compParamm(CompParams)
# 
#     function func(du, u, r, t)
#         TAMode.TAMreact(du, u, r, t, reaction = false)
#     end
# 
#     prob = ODEProblem(func, solInit, maximum(tps), pp)
#     solut = solve(prob, AutoTsit5(Rodas5()); saveat = tps, reltol = TAMode.solTol, isoutofdomain = TAMode.domainDef).u
# 
#     #TODO: Write test
# end


# @testset "Can successfully assemble the comp model parameters." begin
#     TAMode.compTAM(tps, CompParams)
# 
#     @time TAMode.compTAM(tps, CompParams)
# end


# @testset "With no ligand nothing changes." begin
#     pp = copy(CompParams)
#     pp[3] = 0.0
# 
#     solOut = TAMode.compTAM(tps, pp)
# 
#     solOut ./= mean(solOut, dims = 1) .+ eps()
# 
#     @test all(aboutZero.(var(solOut, dims = 1)))
# end
