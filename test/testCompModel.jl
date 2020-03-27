
CompParams = ones(12) * 0.5

@testset "Can successfully assemble the comp model parameters." begin
    TAMode.compTAM(tps, CompParams)

    @time TAMode.compTAM(tps, CompParams)
end
