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
recpSpecific = [
    vcat(ones(12), 0, zeros(24), 0.5 * ones(6), zeros(6), 0.5 * ones(6)),  # AXL
    vcat(zeros(12), 0, ones(12), zeros(12), 0.5 * ones(12), zeros(6)),  # MerTK
    vcat(zeros(12), 0, zeros(12), ones(12), zeros(6), 0.5 * ones(12)),
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
function het_module(Rone, Rtwo, dRone, dRtwo, hetR, hetDim, dhetDim, tr::Rates, Gas, dLi, dRr)
    dRr[1] = tr.xFwd * Rone[3] * Rtwo[3] - hetR.xRev[1] * hetDim[3]
    dRr[2] = tr.xFwd * Rone[1] * Rtwo[4] - hetR.xRev[2] * hetDim[3]
    dRr[3] = tr.xFwd * Rone[2] * Rtwo[2] - hetR.xRev[3] * hetDim[3]
    dRr[4] = tr.xFwd * Rone[4] * Rtwo[1] - hetR.xRev[4] * hetDim[3]
    dRr[5] = tr.xFwd * Rone[3] * Rtwo[1] - hetR.xRev[5] * hetDim[2]
    dRr[6] = tr.xFwd * Rone[1] * Rtwo[2] - hetR.xRev[6] * hetDim[2]
    dRr[7] = tr.xFwd * Rone[2] * Rtwo[1] - hetR.xRev[7] * hetDim[1]
    dRr[8] = tr.xFwd * Rone[1] * Rtwo[3] - hetR.xRev[8] * hetDim[1]
    dRr[9] = hetR.xFwd15 * Gas * hetDim[1] - hetR.xRev[9] * hetDim[3]
    dRr[10] = hetR.xFwd16 * Gas * hetDim[2] - hetR.xRev[10] * hetDim[3]

    dRone[1] += -dRr[2] - dRr[6] - dRr[8]
    dRone[2] += -dRr[3] - dRr[7]
    dRone[3] += -dRr[1] - dRr[5]
    dRone[4] += -dRr[4]

    dRtwo[1] += -dRr[4] - dRr[5] - dRr[7]
    dRtwo[2] += -dRr[3] - dRr[6]
    dRtwo[3] += -dRr[1] - dRr[8]
    dRtwo[4] += -dRr[2]

    dhetDim[1] += dRr[7] + dRr[8] - dRr[9] # AMD1
    dhetDim[2] += dRr[5] + dRr[6] - dRr[10] # MAD1
    dhetDim[3] += dRr[1] + dRr[2] + dRr[3] + dRr[4] + dRr[9] + dRr[10] # AMD2

    if !(dLi isa Nothing)
        dLi[1] += -dRr[9] - dRr[10]
    end

    nothing
end


" Tracks the individual receptor-ligand interactions within a given receptor. "
function react_module(R, dR, Gas, dLi, r::TAMrates, tr::Rates, dRr)
    dRr[1] = r.binding[1] * R[1] * Gas - r.binding[2] * R[2]
    dRr[2] = r.binding[3] * R[1] * Gas - r.binding[4] * R[3]
    dRr[3] = r.binding[3] * R[2] * Gas - r.binding[4] * R[4]
    dRr[4] = r.binding[1] * R[3] * Gas - r.binding[2] * R[4]
    dRr[5] = tr.xFwd * R[1] * R[2] - r.xRev[1] * R[5]
    dRr[6] = tr.xFwd * R[1] * R[3] - r.xRev[2] * R[5]
    dRr[7] = tr.xFwd * R[1] * R[4] - r.xRev[3] * R[6]
    dRr[8] = tr.xFwd * R[2] * R[2] - r.xRev[4] * R[6]
    dRr[9] = tr.xFwd * R[3] * R[3] - r.xRev[5] * R[6]
    dRr[10] = r.xFwd6 * Gas * R[5] - r.xRev[6] * R[6]

    dR[1] += -dRr[7] - dRr[6] - dRr[5] - dRr[1] - dRr[2] # AXL
    dR[2] += -2 * (dRr[8]) - dRr[5] + dRr[1] - dRr[3] # AXLgas1
    dR[3] += -2 * (dRr[9]) - dRr[6] + dRr[2] - dRr[4] # AXLgas2
    dR[4] += -dRr[7] + dRr[3] + dRr[4] # AXLgas12
    dR[5] += -dRr[10] + dRr[6] + dRr[5] # AXLdimer1
    dR[6] += dRr[10] + dRr[9] + dRr[8] + dRr[7] # AXLdimer2

    if !(dLi isa Nothing)
        dLi[1] += -dRr[10] - dRr[1] - dRr[2] - dRr[3] - dRr[4]
    end

    nothing
end


" Handles trafficking of receptor and ligand. "
function trafFunc(dextR, dintR, intRate, extR, intR, kRec, kDeg, fElse)
    dextR .+= -extR * intRate + kRec * (1 - fElse) * intR * internalFrac # Endocytosis, recycling
    dintR .+= extR * intRate / internalFrac - kRec * (1 - fElse) * intR - kDeg * fElse * intR # Endocytosis, recycling, degradation
end


function compartmentReact(u, du, Gas, dLi, r, cache)
    for (ii, aa) in enumerate((1, 7, 13))
        dR = view(du, aa:(aa + 5))

        react_module(view(u, aa:(aa + 5)), dR, Gas, dLi, r.TAMs[ii], r, cache)

        dR[1] += r.TAMs[ii].expression
    end

    het_module(view(u, 1:6), view(u, 7:12), view(du, 1:6), view(du, 7:12), r.hetR.AM, view(u, 19:21), view(du, 19:21), r, Gas, dLi, cache)
    het_module(view(u, 7:12), view(u, 13:18), view(du, 7:12), view(du, 13:18), r.hetR.MT, view(u, 22:24), view(du, 22:24), r, Gas, dLi, cache)
    het_module(view(u, 1:12), view(u, 13:18), view(du, 1:12), view(du, 13:18), r.hetR.AT, view(u, 25:27), view(du, 25:27), r, Gas, dLi, cache)
end


function TAMreact(du, u, r::Rates, t)
    fill!(du, 0.0)
    cache = Vector{promote_type(eltype(du), typeof(r.xFwd))}(undef, 10)

    compartmentReact(view(u, 1:27), view(du, 1:27), r.gasCur, nothing, r, cache)
    compartmentReact(view(u, 28:54), view(du, 28:54), u[end] / internalV, view(du, Nspecies), r, cache)

    for ii in 1:27
        if pYc[ii] == 1
            fElse = 1.0
        else
            fElse = r.fElse
        end

        trafFunc(view(du, ii), view(du, ii + 27), r.internalize, view(u, ii), view(u, ii + 27), r.kRec, r.kDeg, fElse)
    end

    du[end] = -r.kDeg * u[end] # Gas6 degradation

    nothing
end


function detailedBalance(out::Union{Rates{TT}, comprates{TT}})::Rates{TT} where {TT <: Real}
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
