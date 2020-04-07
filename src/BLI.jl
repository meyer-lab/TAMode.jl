function importData(cond)
    filepath = joinpath(dirname(pathof(TAMode)), "..", "data", cond)

    df = CSV.read(filepath)
    conc = df[1, 2:end]
    tps = df[5:end, 1]
    measVal = parse.(Float64, Matrix(df[5:end, 2:end]))
    measVal .-= minimum(measVal, dims=2)
    return parse.(Float64, Array(conc)), parse.(Float64, tps), measVal
end


" The actual binding model calculation. This assumes a single site. "
function bindingCalc(tps::Vector, Kon::Real, Kdis::Real, conc::Vector)::Matrix
    Tshift = tps[1] + 599.9
    tBind = tps[tps .< Tshift] .- tps[1]
    tUnbind = tps[tps .> Tshift] .- Tshift
    KD = Kdis / Kon
    conc = reshape(conc, (1, :))

    bind_step = conc ./ (KD .+ conc) .* (1 .- (1 ./ (exp.((Kon .* conc .+ Kdis) .* tBind))))
    bind_end = conc ./ (KD .+ conc) .* (1 .- (1 ./ (exp.((Kon .* conc .+ Kdis) .* 599.9))))
    unbind_step = bind_end .* exp.(-Kdis .* tUnbind)

    return vcat(bind_step, unbind_step)
end


@model BLI(tps, conc, bindData, ::Type{TV} = Vector{Float64}) where {TV} = begin
    Kon ~ LogNormal(9.0, 1.0)
    Kdis ~ LogNormal(-1.0, 1.0)
    Rmax ~ LogNormal(-1.0, 0.1)
    stdev = 0.03

    resid_save = bindingCalc(tps, Kon, Kdis, conc)

    residNorm = norm(bindData - resid_save * Rmax) / stdev
    residNorm ~ Chisq(length(bindData))
end


function sampleModel(pathIn; testt=false)
    conc, tps, bindData = TAMode.importData(pathIn)
    model = BLI(tps, conc, Matrix(bindData))

    if testt
        return sample(model, HMC(0.001, 4), 5)
    end

    return sample(model, NUTS(), 500)
end
