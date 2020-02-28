
@testset "Test if bothLigands function works." begin
    pp = TAMode.Lsparam(fill(0.2, 9))
    ss = TAMode.getAutocrine(params, TAMreactLS, 30)
end


@testset "Make sure code LS upholds mass conservation." begin
    tt = TAMode.Lsparam(fill(0.2, 9))
    firstV = TAMode.getAutocrine(params, TAMreactLS, 30)

    tt.kDeg = 0
    tt.expression = 0
    tt.curL = (100, 100)

    secondV = TAMode.runTAMinitLS([1000000.0], tt, firstV)

    @test dot(firstV, TAMode.totalLS) - dot(secondV, TAMode.totalLS) < 0.0001
end
