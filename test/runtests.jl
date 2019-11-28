using Test
using Profile
using TAMode
using LinearAlgebra

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


@testset "Make sure code upholds mass conservation." begin
    tt = TAMode.param(params)
    firstV = TAMode.getAutocrine(tt)
    
    tt.kDeg = 0
    tt.TAMs[1].expression = 0
    tt.TAMs[2].expression = 0
    tt.TAMs[3].expression = 0
    tt.gasCur *= 1000

    secondV = TAMode.runTAMinit([1000000.0], tt, firstV)

    @test dot(firstV, TAMode.total) ≈ dot(secondV, TAMode.total)
end

@testset "Test Amount" begin
    
    tt = TAMode.param(params)
    tt.gasCur = 0.0
    tt.TAMs[1].expression = 10.0
    tt.TAMS[2].expression = 10.0
    tt.TAMS[3].expression = 10.0
    tt.kDeg = 1e8
    tt.kRec = 0.0
    tt.internalize = 10.0
    
    outt = TAMode.getAutocrine(tt)
    
    @test outt .* fgMgConv .< 1e-3  #ask if this is used anywhere else
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

@testset "Ensure that system reaches equilibrium." begin  
    uLong = TAMode.getAutocrine(params)
    dnorm = zeros(55)
    TAMode.TAM_reacti(dnorm, uLong, params, 0.0)
    @test all(dnorm .< 0.05)
end


