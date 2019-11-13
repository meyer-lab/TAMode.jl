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
    uLong = TAMode.runTAM([10.0], params, 100.0)

    dnorm = TAMode.TAM_reacti_dnorm(zeros(55), uLong, params, 0.0)

    # TODO: Lower tolerance
    @test dnorm < 2.0
end

@testset "Swapping Ig1 with Ig2 doesn't change anything." begin
    
    rr::Rates = TAMode.param(paramR.getTrafP()) # getTrafP in Distribution.hpp
    rrSwap::Rates = TAMode.swapIgs(rr)
    
    tt::TAM = TAM(rr) # TAM constructor
    data::TAMout = tt.calcStim(tps, 10)
    
    tt::TAM = TAM(rrSwap)
    dataSwap::TAMout = tt.calcStim(tps, 10)
    
    for ii in 1:size(data.pY)
        @test ≈ ((data.pY[ii] - dataSwap.pY[ii])/(data.pY[ii] + dataSwap.pY[ii] + 1e-6), 0, 1e-6)
        @test ≈ ((data.total[ii] - dataSwap.total[ii])/(data.total[ii] + dataSwap.total[ii] + 1e-6), 0, 1e-6)
        @test ≈ ((data.surf[ii] - dataSwap.surf[ii])/(data.surf[ii] + dataSwap.surf[ii] + 1e-6), 0, 1e-6)
        @test ≈ ("surfpY", (data.surfPY[ii] - dataSwap.surfPY[ii])/(data.surfPY[ii] + dataSwap.surfPY[ii] + 1e-6), 0, 1e-6) # assert with message?
    
    end
    
    for ii in 1:size(data.surfL)
        @test ≈ ((data.surfL[ii] - dataSwap.surfL[ii])/(data.surfL[ii] + dataSwap.surfL[ii] + 1e-6), 0, 1e-6)
    end
end

