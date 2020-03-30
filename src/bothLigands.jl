
function Lsparam(params::Vector)
    fwdDef = 0.6
    @assert all(params .>= 0.0)

    out = Lsrates{eltype(params)}(
        (fwdDef, params[5], fwdDef, params[6]),
        (fwdDef, params[7], fwdDef, params[8]),
        zeros(16),
        0.0,
        0.0,
        5.8e-2,
        2.2e-3, # kRec, kDeg
        0.2,
        0.03,
        0.3, # fElse, internalize, pYinternalize
        params[4], # expression
        params[2:3],
        (0.0, 0.0),
        params[1],
    )

    out.xRev[1] = params[9]

    out.xRev[3] = out.xRev[1] * out.GBinding[2] / out.GBinding[4] * out.GBinding[2] / out.GBinding[4]
    out.xRev[8] = out.xRev[1] / out.GBinding[4] / out.GBinding[4] * out.PBinding[4] * out.PBinding[4]
    out.xRev[9] = out.xRev[1] / out.GBinding[4] / out.GBinding[4] * out.PBinding[2] * out.PBinding[4]
    out.xRev[2] = out.xRev[1] * out.GBinding[2] / out.GBinding[4]
    out.xRev[10] = out.xRev[1] / out.GBinding[4] / out.GBinding[4] * out.PBinding[2] * out.PBinding[2]
    out.xRev[11] = out.xRev[1] / out.GBinding[4] / out.GBinding[4] * out.PBinding[4] * out.GBinding[2]
    out.xRev[12] = out.xRev[1] / out.GBinding[4] / out.GBinding[4] * out.PBinding[2] * out.GBinding[4]
    out.xRev[13] = out.xRev[1] / out.GBinding[4] / out.GBinding[4] * out.PBinding[4] * out.GBinding[4]
    out.xRev[14] = out.xRev[1] / out.GBinding[4] / out.GBinding[4] * out.PBinding[2] * out.GBinding[2]
    out.xFwd27 = max(out.PBinding[1], out.PBinding[3])
    out.xFwd29 = max(out.GBinding[1], out.GBinding[3])
    out.xRev[4] = out.GBinding[4]
    out.xRev[5] = out.GBinding[2]
    out.xRev[6] = out.PBinding[2]
    out.xRev[7] = out.PBinding[4]
    out.xRev[15] = out.xFwd27 * out.xRev[8] / out.xRev[7] * out.PBinding[2] / out.PBinding[1]
    out.xRev[16] = out.xFwd29 * out.xRev[1] / out.xRev[4] * out.GBinding[2] / out.GBinding[1]

    return out
end


function react_module(R, dR, curL, r::Lsrates, dRr)
    dRr[1] = r.GBinding[1] * curL[1] * R[1] - r.GBinding[2] * R[2]
    dRr[2] = r.GBinding[3] * curL[1] * R[1] - r.GBinding[4] * R[3]
    dRr[3] = r.PBinding[1] * curL[2] * R[1] - r.PBinding[2] * R[4]
    dRr[4] = r.PBinding[3] * curL[2] * R[1] - r.PBinding[4] * R[5]
    dRr[5] = r.GBinding[3] * curL[1] * R[2] - r.GBinding[4] * R[6]
    dRr[6] = r.GBinding[1] * curL[1] * R[3] - r.GBinding[2] * R[6]
    dRr[7] = r.PBinding[1] * curL[2] * R[3] - r.PBinding[2] * R[7]
    dRr[8] = r.GBinding[3] * curL[1] * R[4] - r.GBinding[4] * R[7]
    dRr[9] = r.PBinding[3] * curL[2] * R[4] - r.PBinding[4] * R[8]
    dRr[10] = r.PBinding[1] * curL[2] * R[5] - r.PBinding[2] * R[8]
    dRr[11] = r.GBinding[1] * curL[1] * R[5] - r.GBinding[2] * R[9]
    dRr[12] = r.PBinding[3] * curL[2] * R[2] - r.PBinding[4] * R[9]
    dRr[13] = r.xFwd * R[2] * R[2] - r.xRev[1] * R[10]
    dRr[14] = r.xFwd * R[1] * R[6] - r.xRev[2] * R[10]
    dRr[15] = r.xFwd * R[3] * R[3] - r.xRev[3] * R[10]
    dRr[16] = r.xFwd * R[2] * R[1] - r.xRev[4] * R[11]
    dRr[17] = r.xFwd * R[3] * R[1] - r.xRev[5] * R[11]
    dRr[18] = r.xFwd * R[1] * R[5] - r.xRev[6] * R[12]
    dRr[19] = r.xFwd * R[1] * R[4] - r.xRev[7] * R[12]
    dRr[20] = r.xFwd * R[4] * R[4] - r.xRev[8] * R[13]
    dRr[21] = r.xFwd * R[1] * R[8] - r.xRev[9] * R[13]
    dRr[22] = r.xFwd * R[5] * R[5] - r.xRev[10] * R[13]
    dRr[23] = r.xFwd * R[7] * R[1] - r.xRev[11] * R[14]
    dRr[24] = r.xFwd * R[9] * R[1] - r.xRev[12] * R[14]
    dRr[25] = r.xFwd * R[2] * R[4] - r.xRev[13] * R[14]
    dRr[26] = r.xFwd * R[3] * R[5] - r.xRev[14] * R[14]
    dRr[27] = r.xFwd27 * R[12] * curL[2] - r.xRev[15] * R[13]
    dRr[28] = r.xFwd27 * R[11] * curL[2] - r.xRev[15] * R[14]
    dRr[29] = r.xFwd29 * R[12] * curL[1] - r.xRev[16] * R[14]
    dRr[30] = r.xFwd29 * R[11] * curL[1] - r.xRev[16] * R[10]

    dR[1] = -sum(dRr[SVector(1, 2, 3, 4, 14, 16, 17, 18, 19, 21, 23, 24)]) # A
    dR[2] = dRr[1] - sum(dRr[SVector(5, 12, 13, 13, 16, 25)]) # AG (first site)
    dR[3] = dRr[2] - sum(dRr[SVector(6, 7, 15, 15, 17, 26)]) # AG (second site)
    dR[4] = dRr[3] - sum(dRr[SVector(8, 9, 19, 20, 20, 25)]) # AP (first site)
    dR[5] = dRr[4] - sum(dRr[SVector(10, 11, 18, 22, 22, 26)]) # AP (second site)
    dR[6] = dRr[5] + dRr[6] - dRr[14] # AGG
    dR[7] = dRr[7] + dRr[8] - dRr[23] # APG
    dR[8] = dRr[9] + dRr[10] - dRr[21] # APP
    dR[9] = dRr[11] + dRr[12] - dRr[24] # AGP

    dR[10] = sum(dRr[SVector(13, 14, 15, 30)]) # AGGA
    dR[11] = dRr[16] + dRr[17] - dRr[28] - dRr[30] # AGA
    dR[12] = dRr[18] + dRr[19] - dRr[27] - dRr[29] # APA
    dR[13] = sum(dRr[SVector(20, 21, 22, 27)]) # APPA
    dR[14] = sum(dRr[SVector(23, 24, 25, 26, 28, 29)]) # APGA

    if length(dR) == 16
        dR[15] = -sum(dRr[SVector(1, 2, 5, 6, 8, 11, 29, 30)]) - r.kDeg * curL[1] # Gasi
        dR[16] = -sum(dRr[SVector(3, 4, 7, 9, 10, 12, 27, 28)]) - r.kDeg * curL[2] # PROSi
    end

    nothing
end


totalLS = vcat(ones(9), 2 * ones(5), ones(9) * internalFrac, 2 * ones(5) * internalFrac, zeros(2))
surfaceLS = vcat(ones(14), zeros(16))
pYcLS = vcat(zeros(9), ones(5))


function TAMreact(dR, R, tr::Lsrates, t)
    fill!(dR, 0.0)
    cache = Vector{promote_type(eltype(dR), typeof(tr.xFwd))}(undef, 30)

    react_module(view(R, 1:14), view(dR, 1:14), tr.curL, tr, cache)
    react_module(view(R, 15:28), view(dR, 15:30), view(R, 29:30) / internalV, tr, cache)

    dR[1] += tr.expression

    trafFunc(dR, R, tr, pYcLS)

    nothing
end
