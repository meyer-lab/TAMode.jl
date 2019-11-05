

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
