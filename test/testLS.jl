
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

    secondV = TAMode.runTAMinit([1.0e6], tt, TAMode.TAMreactLS, firstV)

    @test dot(firstV, TAMode.totalLS) - dot(secondV, TAMode.totalLS) < 0.0001
end
