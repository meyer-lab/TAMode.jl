
@testset "Test if bothLigands function works." begin
    pp = TAMode.Lsparam(fill(0.2, 9))
    ss = TAMode.getAutocrine(pp, TAMode.TAMreactLS, 30)
end


@testset "Make sure code LS upholds mass conservation." begin
    tt = TAMode.Lsparam(fill(0.2, 9))
    firstV = TAMode.getAutocrine(tt, TAMode.TAMreactLS, 30)

    tt.kDeg = 0
    tt.expression = 0
    tt.curL = (100, 100)

    secondV = TAMode.runTAMinit([1000000.0], tt, TAMode.TAMreactLS, firstV)

    @test dot(firstV, TAMode.totalLS) - dot(secondV, TAMode.totalLS) < 0.0001
end


@testset "LS - Make sure that TAM surface don't explode at long time in reaction code." begin
    tt = TAMode.Lsparam(fill(0.2, 9))
    firstSurf = TAMode.getAutocrine(tt, TAMode.TAMreactLS, 30)

    tt.expression = 0.0
    tt.kRec = 0.0
    tt.internalize = 0.0
    tt.pYinternalize = 0.0
    tt.curL *= 10.0

    secondSurf = TAMode.runTAMinit([100.0], tt, TAMode.TAMreactLS, firstSurf)

    @test_broken isapprox(dot(firstSurf, TAMode.surface), dot(secondSurf, TAMode.surface), rtol = 1.0e-5)
end

#failed but no error
@testset "LS - Ensure that system reaches equilibrium." begin
    tt = TAMode.Lsparam(fill(0.2, 9))
    
    uLong = TAMode.getAutocrine(tt, TAMode.TAMreactLS, 30)
    dnorm = zeros(30)
    TAMode.TAMreactLS(dnorm, uLong, tt, 0.0)
    @test all(dnorm .< 0.05)
end
