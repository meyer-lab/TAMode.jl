using Test
using Profile
using TAMode
using LinearAlgebra
using ForwardDiff

tps = 10.0 .^ range(-6.0, stop = 4.0, length = 25)
params = ones(15) * 0.5


aboutZero = x -> isapprox(x, 0.0, rtol = 1.0e-5, atol = 1.0e-5)


@testset "Can successfully assemble the parameters." begin
    TAMode.runTAM(tps, params, 10.0)

    @time TAMode.runTAM(tps, params, 10.0)
end


@testset "Check that we can diff though the solver." begin
    g = ForwardDiff.gradient(x -> dot(TAMode.runTAM([10.0], x, 1.0), TAMode.pY), params)
    @test(length(g) == length(params))
end


@testset "Make sure code upholds mass conservation." begin
    tt = TAMode.param(params)
    firstV = TAMode.getAutocrine(tt)

    tt.kDeg = 0
    tt.TAMs[1].expression = 0
    tt.TAMs[2].expression = 0
    tt.TAMs[3].expression = 0
    tt.gasCur *= 1000

    secondV = TAMode.runTAMinit([1000000.0], tt, firstV)

    @test dot(firstV, TAMode.total) - dot(secondV, TAMode.total) < 0.0001
end


@testset "Test that receptor and ligand amounts match expectations." begin
    tt = TAMode.param(params)

    tt.gasCur = 0.0
    tt.TAMs[1].expression = 10.0
    tt.TAMs[2].expression = 10.0
    tt.TAMs[3].expression = 10.0
    tt.kDeg = 1e8
    tt.kRec = 0.0
    tt.internalize = 10.0

    outt = TAMode.getAutocrine(tt)

    # Expect no ligand
    @test outt[13] ≈ 0.0

    for i = 1:3
        @test isapprox(dot(outt, TAMode.recpSpecific[i] .* TAMode.total), 1.0, rtol = 1.0e-5)
    end
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


@testset "Swapping twice leads to the same rates." begin
    rr = TAMode.param(params)
    rrB = TAMode.swapIgs(TAMode.swapIgs(rr))

    for ii = 1:3
        @test all(rr.TAMs[ii].xRev .≈ rrB.TAMs[ii].xRev)
        @test all(rr.TAMs[ii].binding .≈ rrB.TAMs[ii].binding)
        @test all(rr.hetR[ii].xRev .≈ rrB.hetR[ii].xRev)
    end
end


@testset "Swapping Ig1 with Ig2 doesn't change anything." begin
    # TODO: boundLig doesn't work here...
    rr = TAMode.param(params)
    rrSwap = TAMode.swapIgs(rr)

    dataDiff = TAMode.runTAMinit(tps, rr, zeros(55)) .- TAMode.runTAMinit(tps, rrSwap, zeros(55))
    @test all(aboutZero.(dataDiff * TAMode.pY))
    @test all(aboutZero.(dataDiff * TAMode.total))
    @test all(aboutZero.(dataDiff * TAMode.surface))

    for ii = 1:3
        @test all(aboutZero.(dataDiff * TAMode.recpSpecific[ii]))
    end

    autoDiff = TAMode.getAutocrine(rr) .- TAMode.getAutocrine(rrSwap)
    @test aboutZero(dot(autoDiff, TAMode.pY))
    @test aboutZero(dot(autoDiff, TAMode.total))
    @test aboutZero(dot(autoDiff, TAMode.surface))

    for ii = 1:3
        @test aboutZero(dot(autoDiff, TAMode.recpSpecific[ii]))
    end

    dataDiff = TAMode.runTAM(tps, rr, 1.0) .- TAMode.runTAM(tps, rrSwap, 1.0)
    @test all(aboutZero.(dataDiff * TAMode.pY))
    @test all(aboutZero.(dataDiff * TAMode.total))
    @test all(aboutZero.(dataDiff * TAMode.surface))

    for ii = 1:3
        @test all(aboutZero.(dataDiff * TAMode.recpSpecific[ii]))
    end
end


@testset "Test that if none of the ligand is expressed, we don't end up seeing any." begin
    rr = TAMode.param(params)
    rr.gasCur = 0.0

    data = TAMode.runTAM(tps, rr, 0.0)

    @test all(data * TAMode.pY .≈ 0)
    @test all(data * TAMode.boundLig .≈ 0)
end


@testset "Check if there's no receptor that we don't see any." begin
    for i = 1:3
        rr = TAMode.param(params)
        rr.TAMs[i].expression = 0.0

        data = TAMode.runTAM(tps, rr, 1.0)

        @test all(data * TAMode.recpSpecific[i] .≈ 0.0)
    end
end


@testset "Make sure that TAM surface don't explode at long time in reaction code." begin
    tt = TAMode.param(params)
    firstSurf = TAMode.getAutocrine(tt)

    tt.TAMs[1].expression = 0.0
    tt.TAMs[2].expression = 0.0
    tt.TAMs[3].expression = 0.0
    tt.kRec = 0.0
    tt.internalize = 0.0
    tt.pYinternalize = 0.0
    tt.gasCur *= 10.0

    secondSurf = TAMode.runTAMinit([100.0], tt, firstSurf)

    @test_broken isapprox(dot(firstSurf, TAMode.surface), dot(secondSurf, TAMode.surface), rtol = 1.0e-5)
end


@testset "Ensure that system reaches equilibrium." begin
    uLong = TAMode.getAutocrine(params)
    dnorm = zeros(55)
    TAMode.TAM_reacti(dnorm, uLong, params, 0.0)
    @test all(dnorm .< 0.05)
end
