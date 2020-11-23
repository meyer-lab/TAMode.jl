mutable struct TAMrates{T}
    binding::MVector{4, T} # fwd/rev binding rate for Ig1, then Ig2
    xRev::MVector{6, T} # xRev 1, 2, 3, 4, 5, 6
    expression::T # AXL expression rate.
    xFwd6::T
end

mutable struct hetRates{T}
    xRev::MVector{10, T} # xRev 7, 8, 9, 10, 11, 12, 13, 14, 15, 16
    xFwd15::T
    xFwd16::T
end


TAMsType{T} = LabelledArrays.@SLVector TAMrates{T} (:Axl, :MerTK, :Tyro3)
hetRType{T} = LabelledArrays.@SLVector hetRates{T} (:AM, :MT, :AT)


mutable struct Rates{T}
    TAMs::TAMsType{T}
    internalize::T # Non-pY species internalization rate.
    pYinternalize::T # pY species internalization rate.
    fElse::T # Recycling fraction for non-D2 species.
    kRec::T # Recycling rate.
    kDeg::T # Degradation rate.
    xFwd::T
    gasCur::T
    hetR::hetRType{T}
end


mutable struct comprates{T}
    TAMs::TAMsType{T}
    diff::T # Receptor diffusivity
    gasPart::T # Partitioning of ligand into PtdSer region
    gasCur::T
    xFwd::T
    hetR::hetRType{T}
end

mutable struct Lsrates{T}
    GBinding::MVector{4, T} # Gas6 fwd/rev rate for Ig1, then Ig2
    PBinding::MVector{4, T} # ProS fwd/rev rate for Ig1, then Ig2
    xRev::MVector{16, T} # Different receptor-receptor unbinding rates

    xFwd27::T
    xFwd29::T

    kRec::T # Recycling rate.
    kDeg::T # Degradation rate.
    fElse::T # Recycling fraction for non-D2 species.
    internalize::T # Non-pY species internalization rate.
    pYinternalize::T # pY species internalization rate.
    expression::T # Receptor expression rate.
    curL::MVector{2, T}
    xFwd::T
end

function flatten(d::Union{TAMrates{T}, hetRates{T}, Rates{T}, comprates{T}, Lsrates{T}, TAMsType{T}, hetRType{T}}) where{T}
    arr::ArrayPartition{T}
    if d isa TAMrates{T}
        arr = ArrayPartition(d.binding, d.xRev, d.expression, d.xFwd6)
    elseif d isa hetRates{T}
        arr = ArrayPartition(d.xRev, d.xFwd15, d.xFwd16)
    elseif d isa Rates{T}
        arr = ArrayPartition(flatten(d.TAMs), d.internalize, d.pYinternalize, d.fElse, d.kRec, d.kDeg, d.xFwd, d.gasCur, flatten(d.hetR))
    elseif d isa comprates{T}
        arr = ArrayPartition(flatten(d.TAMs), d.diff, d.gasPart, d.gasCur, d.xFwd, flatten(d.hetR))
    elseif d isa Lsrates{T}
        arr = ArrayPartition(d.GBinding, d.PBinding, d.xRev, d.xFwd27, d.xFwd29, d.kRec, d.kDeg, d.fElse, d.internalize, d.pYinternalize, d.expression, d.curL, d.xFwd)
    elseif d isa TAMsType{T}
        arr = ArrayPartition(flatten(d.Axl), flatten(d.MerTK), flatten(d.Tyro3))
    elseif d isa hetRates{T}
        arr = ArrayPartition(flatten(d.AM), flatten(d.MT), flatten(d.TM))
    return arr
end

#works for TAMrates, hetRates, LsRates
function unflatten(a::ArrayPartition{T}, type::DataType) where{T}
    return type(a.x...)
end
