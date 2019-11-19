using Test
using Profile
using TAMode
using LinearAlgebra

<<<<<<< HEAD
tps = [0.1, 1.0, 10.0, 100.0, 1000.0]
=======
tps = [0.1, 1.0, 10.0, 100.0]
>>>>>>> Not done
params = ones(15) * 0.5


@testset "Can successfully assemble the parameters." begin
    TAMode.runTAM(tps, params, 1.0)
    @time TAMode.runTAM(tps, params, 1.0)

    @profile TAMode.runTAM(tps, params, 1.0)
<<<<<<< HEAD
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

    @test dot(firstV, TAMode.total) - dot(secondV, TAMode.total) < 0.0001
=======

    Profile.print(noisefloor=2.0)
>>>>>>> Not done
end


@testset "Check for detailed balance at steady-state." begin
<<<<<<< HEAD
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
=======
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



>>>>>>> Not done
