
@testset "Test if bothLigands function works." begin
    pp = TAMode.Lsparam(fill(0.2, 9))
    ss = TAMode.getAutocrineLS(pp)
end


@testset "Make sure code LS upholds mass conservation." begin
    tt = TAMode.Lsparam(fill(0.2, 9))
    firstV = TAMode.getAutocrineLS(tt)

    tt.kDeg = 0
    tt.expression = 0
    tt.curL = (100, 1000)

    secondV = TAMode.runTAMinitLS([1000000.0], tt, firstV)

    @test dot(firstV, TAMode.totalLS) - dot(secondV, TAMode.totalLS) < 0.0001
end
