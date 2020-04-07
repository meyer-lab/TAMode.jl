function importData(cond)
    if cond == :T1
        condPath = "T1-010820.csv"
    elseif cond == :T2
        condPath = "T2-010820.csv"
    elseif cond == :TFL
        condPath = "TFL-010820.csv"
    end

    filepath = joinpath(dirname(pathof(TAMode)), "..", "data", condPath)

    df = CSV.read(filepath)
    conc = df[1, 2:end]
    tps = df[5:end, 1]
    measVal = df[5:end, 2:end]
    return parse.(Float64, Array(conc)), parse.(Float64, tps), parse.(Float64, Matrix(measVal))
end


" The actual binding model calculation. This assumes a single site. "
function bindingCalc(tps::Vector, Kon::Real, Kdis::Real, conc::Real)
    Tshift = tps[1] + 599.9
    tBind = tps[tps .< Tshift] .- tps[1]
    tUnbind = tps[tps .> Tshift] .- Tshift

    bind_step = conc ./ (Kdis ./ Kon + conc) .* (1 .- (1 ./ (exp.((Kon .* conc .+ Kdis) .* tBind))))
    unbind_step = bind_step[end] .* exp.(-Kdis .* tUnbind)
    theor_bind = vcat(bind_step[:], unbind_step)

    return theor_bind
end


@model BLI(tps, conc, bindData, ::Type{TV} = Vector{Float64}) where {TV} = begin
    Kon ~ LogNormal(6.0, 0.5)
    Kdis ~ LogNormal(1.0, 1.0)
    Rmax ~ LogNormal(-1.0, 0.1)
    stdev = 0.03

    resid_save = Matrix{typeof(Kon)}(undef, length(tps), length(conc))
    for i = 1:length(conc)
        resid_save[:, i] .= bindingCalc(tps, Kon, Kdis, conc[i])
    end

    residNorm = norm(bindData - resid_save * Rmax) / stdev
    residNorm ~ Chisq(length(bindData))
end


function buildModel(pathIn)
    conc, tps, bindData = TAMode.importData(pathIn)

    return BLI(tps, conc, Matrix(bindData))
end


function sampleModel(pathIn; testt=false)
    model = buildModel(pathIn)

    if testt
        return sample(model, HMC(0.001, 4), 5)
    end

    return sample(model, NUTS(), 500)
end
