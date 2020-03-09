
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

@testset "Ensure that system reaches equilibrium(LS)." begin
    pp = TAMode.Lsparam(LSparams)
    uLong = TAMode.getAutocrine(pp, TAMode.TAMreactLS, 30)
    dnorm = zeros(30)
    TAMode.TAMreactLS(dnorm, uLong, pp, 0.0)
    @test all(dnorm .< 0.05)
end

@testset "Make sure that TAM surface don't explode at long time in reaction code (LS)." begin
    tt = TAMode.Lsparam(LSparams)
    firstSurf = TAMode.getAutocrine(tt, TAMode.TAMreactLS, 30)

    tt.expression = 0.0
    tt.kRec = 0.0
    tt.internalize = 0.0
    tt.pYinternalize = 0.0
    tt.curL = (10, 10)

    secondSurf = TAMode.runTAMinit([100.0], tt, TAMode.TAMreactLS, firstSurf)

    @test_broken isapprox(dot(firstSurf, TAMode.surface), dot(secondSurf, TAMode.surface), rtol = 1.0e-5)
end