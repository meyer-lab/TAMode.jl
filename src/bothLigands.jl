

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
	internalFrac::T # Ratio of endosomal to surface membrane
	internalV::T # Endosomal volume
	expression::T # Receptor expression rate.
	autocrine::MVector{2, T}
	curL::MVector{2, T}
	xFwd::T
end

# @brief Callback function from CVode to handle ODE solving
# @param[in] R  Double vector for current state of model
# @param dR  Double vector for returning derivative of current state
# @return Always returns 0 for success


function TAM_react(const R::T, const dR::T)
   
    react_module(R, dR, tr.curL) # We give the same internal ligand address as below because it will be overwritten
    react_module(R+iR, dR+iR, (R[15+iR]/tr.internalV, R[16+iR]/tr.internalV))

    dR[1] += tr.expression

    trafFunc(@view dR[1:8], @view dR[(1:8)+(1:8)*R], tr.internalize, R[1:8], R[1:8)+(1:8)*R], tr.kRec, tr.kDeg, tr.fElse, tr.internalFrac)
    
    trafFunc(@view dR[9:13], @view dR[(9:13)+(9:13)*R], tr.pYinternalize, R[9:13], R[(9:13)+(9:13)*R], tr.kRec, tr.kDeg, tr.fElse, tr.internalFrac)
        
    return 0
end