using Test
using Profile
using TAMode

tps = [0.1, 1.0, 10.0, 100.0, 1000.0]
params = ones(15) * 0.5


@testset "Can successfully assemble the parameters." begin
    TAMode.runTAM(tps, params, 1.0)
    @time TAMode.runTAM(tps, params, 1.0)

    @profile TAMode.runTAM(tps, params, 1.0)
    @profile TAMode.runTAM(tps, params, 10.0)
    @profile TAMode.runTAM(tps, params, 100.0)
    @profile TAMode.runTAM(tps, params, 1000.0)

    Profile.print(noisefloor=5.0)
end


@testset "Check for detailed balance at steady-state." begin
    rr = TAMode.param(params)

    rr.TAMs[1].expression = 0.0
    rr.TAMs[2].expression = 0.0
    rr.TAMs[3].expression = 0.0
    rr.kDeg = 0.0

    uLong = TAMode.runTAMinit([1000000.0], rr, TAMode.getAutocrine(params))

    dnorm = TAMode.TAM_reacti_dnorm(zeros(55), uLong, rr, 0.0)

    @test dnorm < 0.05
end


@testset "Swapping Ig1 with Ig2 doesn't change anything." begin
    rr = TAMode.param(params)
    rrSwap = TAMode.swapIgs(rr)

    data = TAMode.runTAM(tps, rr, 10.0)
    dataSwap = TAMode.runTAM(tps, rrSwap, 10.0)

    # Test that the phosphorylated receptor matches up
    @test all(data * TAMode.pY .≈ dataSwap * TAMode.pY)
    @test all(data * TAMode.total .≈ dataSwap * TAMode.total)
    @test all(data * TAMode.surface .≈ dataSwap * TAMode.surface)
    @test all(data * TAMode.boundLig .≈ dataSwap * TAMode.boundLig)
end


@testset "Test that if none of the ligand is expressed, we don't end up seeing any." begin
    rr = TAMode.param(params)
    rr.gasCur = 0.0

    data = TAMode.runTAM(tps, rr, 0.0)

    @test all(data * TAMode.pY .≈ 0)
    @test all(data * TAMode.boundLig .≈ 0)
end


@testset "Check if there's no receptor that we don't see any." begin
    for i in 1:3
        rr = TAMode.param(params)
        rr.TAMs[i].expression = 0.0

        data = TAMode.runTAM(tps, rr, 1.0)

        @test all(data * TAMode.recpSpecific[i] .≈ 0.0)
    end
end
