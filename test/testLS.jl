@testset "Two ligand model." begin
    @testset "LS: Test if bothLigands function works." begin
        pp = TAMode.Lsparam(fill(0.2, 9))
        ss = TAMode.getAutocrine(pp)
    end

    @testset "LS: Make sure code upholds mass conservation." begin
        tt = TAMode.Lsparam(fill(0.2, 9))
        firstV = TAMode.getAutocrine(tt)

        tt.kDeg = 0
        tt.expression = 0
        tt.curL = (100, 100)

        secondV = TAMode.runTAMinit([1.0e6], tt, firstV)

        @test dot(firstV, TAMode.totalLS) - dot(secondV, TAMode.totalLS) < 0.0001
    end

    @testset "LS: Ensure that system reaches detailed balance." begin
        rr = TAMode.Lsparam(fill(0.2, 9))
        rr.curL = (0.0, 0.0)

        autoC = TAMode.getAutocrine(rr)

        rr.expression = 0.0
        rr.kDeg = 0.0

        uLong = TAMode.runTAMinit([1.0e7], rr, autoC)

        # Get the Jacobian matrix
        du = zero(uLong)
        J = ForwardDiff.jacobian((y, x) -> TAMode.TAMreact(y, x, rr, 0.0), du, uLong)
        GK = J * diagm(vec(uLong))

        IDXmismatch = findall(GK - transpose(GK) .≠ 0.0)

        @test length(IDXmismatch) == 2
        @test_broken norm(GK - transpose(GK)) < 1.0e-5
    end

    @testset "LS: Make sure that surface total is preserved without trafficking." begin
        tt = TAMode.Lsparam(fill(0.2, 9))
        firstSurf = TAMode.getAutocrine(tt)

        tt.expression = 0.0
        tt.kRec = 0.0
        tt.internalize = 0.0
        tt.pYinternalize = 0.0
        tt.curL = (10.0, 10.0)

        secondSurf = TAMode.runTAMinit([100.0], tt, firstSurf)

        reduce = TAMode.surfaceLS .* TAMode.totalLS
        @test isapprox(dot(firstSurf, reduce), dot(secondSurf, reduce), rtol = 1.0e-5)
    end

    @testset "LS: Ensure that system reaches equilibrium." begin
        tt = TAMode.Lsparam(fill(0.2, 9))

        uLong = TAMode.getAutocrine(tt)
        dnorm = zeros(30)
        TAMode.TAMreact(dnorm, uLong, tt, 0.0)
        @test norm(dnorm) .< 0.001
    end

    @testset "LS: Test that if no Gas6 is present, we don’t see any." begin
        rr = TAMode.Lsparam(fill(0.2, 9))
        rr.curL = (0.0, 0.0)

        data = TAMode.runTAM(tps, rr, (0.0, 10.0))

        @test all(aboutZero.(data * TAMode.GasLS))
    end

    @testset "LS: Test that if no Protein S is present, we don’t see any." begin
        rr = TAMode.Lsparam(fill(0.2, 9))
        rr.curL = (0.0, 0.0)

        data = TAMode.runTAM(tps, rr, (10.0, 0.0))

        @test all(aboutZero.(data * TAMode.PROSLS))
    end

    @testset "LS: Test that if no ligand is present, we don’t see any pY." begin
        rr = TAMode.Lsparam(fill(0.2, 9))
        rr.curL = (0.0, 0.0)

        data = TAMode.runTAM(tps, rr, (0.0, 0.0))

        @test all(aboutZero.(data * TAMode.pYLS))
        @test all(aboutZero.(data * TAMode.PROSLS))
        @test all(aboutZero.(data * TAMode.GasLS))
    end

    @testset "LS: Test for symmetry wrt to both ligands." begin
        rr = TAMode.Lsparam(fill(0.2, 9))

        Gstim = TAMode.runTAM(tps, rr, (10.0, 0.0))
        Pstim = TAMode.runTAM(tps, rr, (0.0, 10.0))
        dataDiff = Gstim .- Pstim

        @test all(aboutZero.(dataDiff * TAMode.pYLS))
        @test all(aboutZero.(dataDiff * TAMode.totalLS))
        @test all(aboutZero.(dataDiff * TAMode.surfaceLS))
        @test all(aboutZero.(Gstim * TAMode.GasLS .- Pstim * TAMode.PROSLS))
    end
end
