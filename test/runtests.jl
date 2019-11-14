using Test
using Profile
using TAMode

tps = [0.1, 1.0, 10.0, 100.0]
params = ones(15) * 0.5


@testset "Can successfully assemble the parameters." begin
    TAMode.runTAM(tps, params, 1.0)
    @time TAMode.runTAM(tps, params, 1.0)

    @profile TAMode.runTAM(tps, params, 1.0)

    Profile.print(noisefloor=2.0)
end


@testset "Check for detailed balance at steady-state." begin
    uLong = TAMode.runTAM([1000.0], params, 100.0)

    dnorm = TAMode.TAM_reacti_dnorm(zeros(55), uLong, params, 0.0)

    # TODO: Lower tolerance
    @test dnorm < 2.0
end


@testset "Test that if none of the ligand is expressed, we don't end up seeing any." begin
    
    rr = TAMode.param(params)
    
    rr.autocrine = 0
    
    data = TAMode.runTAM(tps, rr, 0.0)
    
    @test all(data * TAMode.pY .≈ 0)
    @test all(data* TAMode.surfL .≈ 0)
    
end

