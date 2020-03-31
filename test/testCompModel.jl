
CompParams = ones(12) * 0.5

@testset "Can successfully assemble the comp model parameters." begin
    TAMode.compTAM(tps, CompParams)

    @time TAMode.compTAM(tps, CompParams)
end


@testset "In the absence of a reaction we just see diffusion." begin
    solInit = zeros(2700)
    solInit[1] = 1000.0

    pp = TAMode.compParamm(CompParams)

    function func(du, u, r, t)
        TAMode.TAMreact(du, u, r, t, reaction = false)
    end

    prob = ODEProblem(func, solInit, maximum(tps), pp)
    solut = solve(prob, AutoTsit5(Rodas5()); saveat = tps).u

    #TODO: Write test
end


@testset "With no ligand nothing changes." begin
    pp = copy(CompParams)
    pp[3] = 0.0

    solOut = TAMode.compTAM(tps, pp)

    solOut ./= mean(solOut, dims = 1) .+ eps()

    @test all(aboutZero.(var(solOut, dims = 1)))
end
