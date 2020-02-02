
CompParams = ones(18) * 0.5

@testset "Can successfully assemble the comp model parameters." begin
    TAMode.calcStim(tps, CompParams, 10.0)

    @time TAMode.calcStim(tps, CompParams, 10.0)
end