function importData(cond)
    filepath = joinpath(dirname(pathof(TAMode)), "..", "data", cond)

    df = CSV.read(filepath)
    conc = df[1, 2:end]
    tps = df[5:end, 1]
    measVal = df[5:end, 2:end]
    return parse.(Float64, Array(conc)), parse.(Float64, tps), parse.(Float64, Matrix(measVal))
end


" The actual binding model calculation. This assumes a single site. "
function bindingCalc(tps::Vector, Kon::Real, Kdis::Real, conc::Vector)::Matrix
    Tshift = tps[1] + 599.9
    tBind = tps[tps .< Tshift] .- tps[1]
    tUnbind = tps[tps .> Tshift] .- Tshift
    KD = Kdis / Kon
    conc = reshape(conc, (1, :))

    bind_step = conc ./ (KD .+ conc) .* (1 .- (1 ./ (exp.((Kon .* conc .+ Kdis) .* tBind))))
    unbind_step = bind_step[end, :] .* exp.(-Kdis .* tUnbind) # one(conc) is to broadcast shape

    return vcat(bind_step, unbind_step)
end


@model BLI(tps, conc, bindData, ::Type{TV} = Vector{Float64}) where {TV} = begin
    Kon ~ LogNormal(6.0, 0.5)
    Kdis ~ LogNormal(1.0, 1.0)
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
