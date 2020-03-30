const internalFrac = 0.5
const internalV = 623.0
const fgMgConv = 135.2
const Nspecies = 55


# Mark surface species
# Variables with C on end are only for one compartment
surface = vcat(ones(27), zeros(28))
pYc = vcat(zeros(4), ones(2), zeros(4), ones(2), zeros(4), ones(2), ones(9))
pY = vcat(pYc, pYc * internalFrac, 0)
ligPiece = [0, 1, 1, 2, 1, 2]
hetPiece = [1, 1, 2]
boundLigC = vcat(ligPiece, ligPiece, ligPiece, hetPiece, hetPiece, hetPiece)
boundLig = vcat(boundLigC, boundLigC * internalFrac, 0)
totalPiece = [1, 1, 1, 1, 2, 2]
totalC = vcat(totalPiece, totalPiece, totalPiece, 2 * ones(9))
total = vcat(totalC, totalC * internalFrac, 0)
recpSpecificC = [
    vcat(ones(6), zeros(12), ones(3), zeros(3), ones(3)),  # AXL
    vcat(zeros(6), ones(6), zeros(6), ones(6), zeros(3)),  # MerTK
    vcat(zeros(6), zeros(6), ones(6), zeros(3), ones(6)),
] # Tyro3
recpSpecific = [
    vcat(recpSpecificC[1], recpSpecificC[1] * internalFrac, 0),  # AXL
    vcat(recpSpecificC[2], recpSpecificC[2] * internalFrac, 0),  # MerTK
    vcat(recpSpecificC[3], recpSpecificC[3] * internalFrac, 0),
] # Tyro3


" This makes a rate for a receptor. "
function paramTAMrate(params::AbstractVector{T}) where {T}
    fBnd = 0.6
    AXL = TAMrates{T}([1.2, 0.042, fBnd, fBnd * params[4]], zeros(6), params[1], 1.50) # From Kariolis et al
    MerTK = TAMrates{T}([fBnd, fBnd * params[5], fBnd, fBnd * params[6]], zeros(6), params[2], 0.0)
    Tyro3 = TAMrates{T}([fBnd, fBnd * params[7], fBnd, fBnd * params[8]], zeros(6), params[3], 0.0)

    TAM = TAMsType{T}(AXL, MerTK, Tyro3)

    hetR = hetRates{T}(zeros(10), 0.0, 0.0)
    hetRs = hetRType{T}(hetR, deepcopy(hetR), deepcopy(hetR))

    TAM.Axl.xRev[5] = 0.0144 # From Kariolis et al

    return TAM, hetRs
end


" Setup the parameters for the full TAM receptor model. "
function param(params::Vector{T})::Rates{T} where {T}
    @assert all(params .>= 0.0)
    @assert params[3] < 1.0

    TAM, hetRs = paramTAMrate(view(params, 8:15))

    out = Rates{T}(TAM, params[1], params[2], params[3], params[4], params[5], params[6], params[7], hetRs)

    return detailedBalance(out)
end


" Tracks the hetero-receptor interactions between two receptors. "
function het_module(Rone, Rtwo, dRone, dRtwo, hetR, hetDim, dhetDim, tr, Gas, dLi)
    dRr1 = tr.xFwd * Rone[3] * Rtwo[3] - hetR.xRev[1] * hetDim[3]
    dRr2 = tr.xFwd * Rone[1] * Rtwo[4] - hetR.xRev[2] * hetDim[3]
    dRr3 = tr.xFwd * Rone[2] * Rtwo[2] - hetR.xRev[3] * hetDim[3]
    dRr4 = tr.xFwd * Rone[4] * Rtwo[1] - hetR.xRev[4] * hetDim[3]
    dRr5 = tr.xFwd * Rone[3] * Rtwo[1] - hetR.xRev[5] * hetDim[2]
    dRr6 = tr.xFwd * Rone[1] * Rtwo[2] - hetR.xRev[6] * hetDim[2]
    dRr7 = tr.xFwd * Rone[2] * Rtwo[1] - hetR.xRev[7] * hetDim[1]
    dRr8 = tr.xFwd * Rone[1] * Rtwo[3] - hetR.xRev[8] * hetDim[1]
    dRr9 = hetR.xFwd15 * Gas * hetDim[1] - hetR.xRev[9] * hetDim[3]
    dRr10 = hetR.xFwd16 * Gas * hetDim[2] - hetR.xRev[10] * hetDim[3]

    dRone[1] += -dRr2 - dRr6 - dRr8
    dRone[2] += -dRr3 - dRr7
    dRone[3] += -dRr1 - dRr5
    dRone[4] += -dRr4

    dRtwo[1] += -dRr4 - dRr5 - dRr7
    dRtwo[2] += -dRr3 - dRr6
    dRtwo[3] += -dRr1 - dRr8
    dRtwo[4] += -dRr2

    dhetDim[1] += dRr7 + dRr8 - dRr9 # AMD1
    dhetDim[2] += dRr5 + dRr6 - dRr10 # MAD1
    dhetDim[3] += dRr1 + dRr2 + dRr3 + dRr4 + dRr9 + dRr10 # AMD2

    if !(dLi isa Nothing)
        dLi[1] += -dRr9 - dRr10
    end

    nothing
end


" Tracks the individual receptor-ligand interactions within a given receptor. "
function react_module(R, dR, Gas, dLi, r::TAMrates, tr)
    dRr1 = r.binding[1] * R[1] * Gas - r.binding[2] * R[2]
    dRr2 = r.binding[3] * R[1] * Gas - r.binding[4] * R[3]
    dRr3 = r.binding[3] * R[2] * Gas - r.binding[4] * R[4]
    dRr4 = r.binding[1] * R[3] * Gas - r.binding[2] * R[4]
    dRr5 = tr.xFwd * R[1] * R[2] - r.xRev[1] * R[5]
    dRr6 = tr.xFwd * R[1] * R[3] - r.xRev[2] * R[5]
    dRr7 = tr.xFwd * R[1] * R[4] - r.xRev[3] * R[6]
    dRr8 = tr.xFwd * R[2] * R[2] - r.xRev[4] * R[6]
    dRr9 = tr.xFwd * R[3] * R[3] - r.xRev[5] * R[6]
    dRr10 = r.xFwd6 * Gas * R[5] - r.xRev[6] * R[6]

    dR[1] += -dRr7 - dRr6 - dRr5 - dRr1 - dRr2 # AXL
    dR[2] += -2 * (dRr8) - dRr5 + dRr1 - dRr3 # AXLgas1
    dR[3] += -2 * (dRr9) - dRr6 + dRr2 - dRr4 # AXLgas2
    dR[4] += -dRr7 + dRr3 + dRr4 # AXLgas12
    dR[5] += -dRr10 + dRr6 + dRr5 # AXLdimer1
    dR[6] += dRr10 + dRr9 + dRr8 + dRr7 # AXLdimer2

    if !(dLi isa Nothing)
        dLi[1] += -dRr10 - dRr1 - dRr2 - dRr3 - dRr4
    end

    nothing
end


" Handles trafficking of receptor and ligand. "
function trafFunc(du, u, r, pYcIn)
    for ii = 1:length(pYcIn)
        if pYcIn[ii] == 1
            fElse = 1.0
            intRate = r.pYinternalize
        else
            fElse = r.fElse
            intRate = r.internalize
        end

        du[ii] += -u[ii] * intRate + r.kRec * (1.0 - fElse) * u[ii + length(pYcIn)] * internalFrac # Endocytosis, recycling
        du[ii + length(pYcIn)] += u[ii] * intRate / internalFrac - r.kRec * (1.0 - fElse) * u[ii + length(pYcIn)] - r.kDeg * fElse * u[ii + length(pYcIn)] # ", degradation
    end
end


function compartmentReact(u, du, Gas, dLi, r)
    for (ii, aa) in enumerate((1, 7, 13))
        dR = view(du, aa:(aa + 5))

        react_module(view(u, aa:(aa + 5)), dR, Gas, dLi, r.TAMs[ii], r)

        dR[1] += r.TAMs[ii].expression
    end

    het_module(view(u, 1:6), view(u, 7:12), view(du, 1:6), view(du, 7:12), r.hetR.AM, view(u, 19:21), view(du, 19:21), r, Gas, dLi)
    het_module(view(u, 7:12), view(u, 13:18), view(du, 7:12), view(du, 13:18), r.hetR.MT, view(u, 22:24), view(du, 22:24), r, Gas, dLi)
    het_module(view(u, 1:6), view(u, 13:18), view(du, 1:6), view(du, 13:18), r.hetR.AT, view(u, 25:27), view(du, 25:27), r, Gas, dLi)
end


function TAMreact(du, u, r::Rates, t)
    fill!(du, 0.0)

    compartmentReact(view(u, 1:27), view(du, 1:27), r.gasCur, nothing, r)
    compartmentReact(view(u, 28:54), view(du, 28:54), u[end] / internalV, view(du, Nspecies), r)

    trafFunc(du, u, r, pYc)

    du[end] = -r.kDeg * u[end] # Gas6 degradation

    nothing
end


function detailedBalance(out)
    for T in out.TAMs
        KD1 = T.binding[2] / T.binding[1]
        KD2 = T.binding[4] / T.binding[3]

        if T.binding[4] <= T.binding[2]
            T.xRev[1] = T.binding[4]
            T.xRev[2] = T.xRev[1] * KD1 / KD2
        else
            T.xRev[2] = T.binding[2]
            T.xRev[1] = T.xRev[2] * KD2 / KD1
        end
    end

    for ii = 2:3
        out.TAMs[ii].xFwd6 = max(out.TAMs[ii].binding[1], out.TAMs[ii].binding[3])
        out.TAMs[ii].xRev[6] =
            out.TAMs[1].xRev[6] * out.TAMs[ii].binding[2] * out.TAMs[ii].binding[4] / out.TAMs[1].binding[4] / out.TAMs[1].binding[2]
    end

    for T in out.TAMs
        KD1 = T.binding[2] / T.binding[1]
        KD2 = T.binding[4] / T.binding[3]

        T.xRev[5] = T.xRev[6] * KD1 * T.xRev[1] / KD2 / KD2 / T.xFwd6
        T.xRev[4] = T.xRev[5] * KD2 * KD2 / KD1 / KD1
        T.xRev[3] = T.xRev[4] * KD1 / KD2
    end

    # ligand binding to the one ligand dimer
    for ii = 1:3
        x = [:Axl, :MerTK, :Axl][ii]
        y = [:MerTK, :Tyro3, :Tyro3][ii]

        KD11 = out.TAMs[x].binding[2] / out.TAMs[x].binding[1]
        KD12 = out.TAMs[y].binding[2] / out.TAMs[y].binding[1]
        KD21 = out.TAMs[x].binding[4] / out.TAMs[x].binding[3]
        KD22 = out.TAMs[y].binding[4] / out.TAMs[y].binding[3]

        out.hetR[ii].xFwd15 = max(out.TAMs[y].binding[1], out.TAMs[x].binding[3])
        out.hetR[ii].xFwd16 = max(out.TAMs[x].binding[1], out.TAMs[y].binding[3])

        out.hetR[ii].xRev[1] = out.TAMs[x].xRev[5] * KD12 / KD11
        out.hetR[ii].xRev[2] = out.hetR[ii].xRev[1] * KD21 / KD12
        out.hetR[ii].xRev[3] = KD22 * KD21 / KD11 / KD12 * out.hetR[ii].xRev[1]
        out.hetR[ii].xRev[4] = out.hetR[ii].xRev[1] * KD22 / KD11

        if out.TAMs[x].binding[2] <= out.TAMs[x].binding[4]
            out.hetR[ii].xRev[8] = out.TAMs[x].binding[2]
            out.hetR[ii].xRev[6] = out.hetR[ii].xRev[8] * KD11 / KD21
        else
            out.hetR[ii].xRev[6] = out.TAMs[x].binding[4]
            out.hetR[ii].xRev[8] = out.hetR[ii].xRev[6] * KD21 / KD11
        end

        if out.TAMs[y].binding[2] <= out.TAMs[y].binding[4]
            out.hetR[ii].xRev[5] = out.TAMs[y].binding[2]
            out.hetR[ii].xRev[7] = out.hetR[ii].xRev[5] * KD22 / KD12
        else
            out.hetR[ii].xRev[7] = out.TAMs[y].binding[4]
            out.hetR[ii].xRev[5] = out.hetR[ii].xRev[7] * KD12 / KD22
        end

        out.hetR[ii].xRev[9] = out.hetR[ii].xRev[1] * KD21 / KD12
        out.hetR[ii].xRev[10] = out.hetR[ii].xRev[3] * KD11 / KD22
    end

    return out
end


function swapIgs(out::Rates{TT})::Rates{TT} where {TT <: Real}
    out = detailedBalance(deepcopy(out))

    for T in out.TAMs
        T.binding[1], T.binding[3] = T.binding[3], T.binding[1]
        T.binding[2], T.binding[4] = T.binding[4], T.binding[2]
        T.xRev[4], T.xRev[5] = T.xRev[5], T.xRev[4]
    end

    return detailedBalance(out)
end
