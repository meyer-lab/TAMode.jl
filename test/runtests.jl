using Test
using Profile
using TAMode
using LinearAlgebra

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

@testset "Make sure code upholds mass conservation." begin
    #total::Pair{MVector{}} #ladkjfals;dkfj
    
    tt = TAMode.param(params)
    tt.gasCur = tt.autocrine
    
    #CPPUNIT_ASSERT_NO_THROW(tt.initSystem());
    #tt.getAutocrine(params, TAM_reactii, 55)
    
    firstTotal = tt.getTot()
    
    tt.kDeg = 0
    tt.TAMs[1].expression = 0
    tt.TAMs[2].expression = 0
    tt.TAMs[3].expression = 0
    tt.gasCur *= 1000 
    
    #CPPUNIT_ASSERT_NO_THROW(tt.runToT(1000));
    
    secondTotal = tt.getTot()
    
    @test all(firstTotal .≈ secondTotal)
end



