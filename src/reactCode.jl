

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


mutable struct Rates{T}
	AXL::TAMrates{T}
	MerTK::TAMrates{T}
	Tyro3::TAMrates{T}
	internalize::T # Non-pY species internalization rate.
	pYinternalize::T # pY species internalization rate.
	fElse::T # Recycling fraction for non-D2 species.
	kRec::T # Recycling rate.
	kDeg::T # Degradation rate.
	xFwd::T
	internalFrac::T
	internalV::T
	gasCur::T
	
	AM::hetRates{T}
	AT::hetRates{T}
	MT::hetRates{T}
end


function param(params)
	@assert all(params .>= 0.0)
	@assert params[3] < 1.0
	fBnd = 0.06

	AXL = TAMrates{eltype(params)}([1.2, 0.042, fBnd, fBnd*params[11]], zeros(6), params[8], 1.50) # From Kariolis et al
	MerTK = TAMrates{eltype(params)}([fBnd, fBnd*params[12], fBnd, fBnd*params[13]], zeros(6), params[9], 0.0)
	Tyro3 = TAMrates{eltype(params)}([fBnd, fBnd*params[14], fBnd, fBnd*params[15]], zeros(6), params[10], 0.0)

	hetR = hetRates{eltype(params)}(zeros(10), 0.0, 0.0)

	out = Rates{eltype(params)}(AXL, MerTK, Tyro3, params[1], params[2], params[3],
								params[4], params[5], params[6], 0.5, 623.0, params[7],
								hetR, deepcopy(hetR), deepcopy(hetR))


	out.AXL.xRev[5] = 0.0144 # From Kariolis et al

	# ==== Detailed balance

	# TODO: Add.

	return out
end


function het_module(Rone, Rtwo, dRone, dRtwo, hetR, hetDim, dhetDim, tr, Gas, dLi)
	dRr = SVector{10}([tr.xFwd*Rone[3]*Rtwo[3] - hetR.xRev[1]*hetDim[3],
					   tr.xFwd*Rone[1]*Rtwo[4] - hetR.xRev[2]*hetDim[3],
					   tr.xFwd*Rone[2]*Rtwo[2] - hetR.xRev[3]*hetDim[3],
					   tr.xFwd*Rone[4]*Rtwo[1] - hetR.xRev[4]*hetDim[3],
					   tr.xFwd*Rone[3]*Rtwo[1] - hetR.xRev[5]*hetDim[2],
					   tr.xFwd*Rone[1]*Rtwo[2] - hetR.xRev[6]*hetDim[2],
					   tr.xFwd*Rone[2]*Rtwo[1] - hetR.xRev[7]*hetDim[1],
					   tr.xFwd*Rone[1]*Rtwo[3] - hetR.xRev[8]*hetDim[1],
					   hetR.xFwd15 * Gas * hetDim[1] - hetR.xRev[9] * hetDim[3],
					   hetR.xFwd16 * Gas * hetDim[2] - hetR.xRev[10] * hetDim[3]])

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
	
	if isnothing(dLi) == false
		dLi[1] += -dRr[9] - dRr[10]
	end
end


function trafFunc(dextR, dintR, intRate::Float64, extR, intR, kRec::Float64, kDeg::Float64, fElse::Float64, internalFrac::Float64)
	dextR[:] .+= -extR*intRate + kRec*(1-fElse)*intR*internalFrac # Endocytosis, recycling
	dintR[:] .+= extR*intRate/internalFrac - kRec*(1-fElse)*intR - kDeg*fElse*intR # Endocytosis, recycling, degradation
end


function heteroTAM(Rone, Rtwo, dRone, dRtwo, hetR, hetDim, dhetDim, tr, Li, dLi)
	het_module(Rone, Rtwo, dRone, dRtwo, hetR, hetDim, dhetDim, tr, tr.gasCur, nothing)
	het_module(view(Rone, 7:10), view(Rtwo, 7:10), view(dRone, 7:10), view(dRtwo, 7:10), hetR, view(hetDim, 4:6), view(dhetDim, 4:6), tr, Li / tr.internalV, dLi);

	trafFunc(view(dhetDim, 1:3), view(dhetDim, 4:6), tr.pYinternalize, hetDim[1:3], hetDim[4:6], tr.kRec, tr.kDeg, 1.0, tr.internalFrac)
end

function TAM_reactii(R, Li, dR, dLi, r::TAMrates, tr::Rates)
	# react_module(R, dR, &temp, tr->gasCur, r, tr);
	# react_module(R+6, dR+6, dLi, (*Li)/tr->internalV, r, tr);
	
	dR[1] += r.expression
	
	trafFunc(view(dR, 1:4), view(dR, 7:10), tr.internalize, R[1:4], R[7:10], tr.kRec, tr.kDeg, tr.fElse, tr.internalFrac)
	trafFunc(view(dR, 4:5), view(dR, 10:11), tr.pYinternalize, R[4:5], R[10:11], tr.kRec, tr.kDeg, 1.0, tr.internalFrac)
end

function TAM_reacti(dxdt_d, x_d, params, t)
	fill!(dxdt_d, 0.0)
	r = param(params)

	TAM_reactii(view(x_d, 1:12), x_d[13], view(dxdt_d, 1:12), view(dxdt_d, 13), r.AXL, r)
	TAM_reactii(view(x_d, 14:25), x_d[13], view(dxdt_d, 14:25), view(dxdt_d, 13), r.MerTK, r)
	TAM_reactii(view(x_d, 26:37), x_d[13], view(dxdt_d, 26:37), view(dxdt_d, 13), r.Tyro3, r)
	
	heteroTAM(x_d[1:12],  x_d[14:25], view(dxdt_d, 1:12),  view(dxdt_d, 14:25), r.AM, x_d[38:43], view(dxdt_d, 38:43), r, x_d[13], view(dxdt_d, 13))
	heteroTAM(x_d[14:25], x_d[26:37], view(dxdt_d, 14:25), view(dxdt_d, 26:37), r.MT, x_d[44:49], view(dxdt_d, 44:49), r, x_d[13], view(dxdt_d, 13))
	heteroTAM(x_d[1:12],  x_d[26:37], view(dxdt_d, 1:12),  view(dxdt_d, 26:37), r.AT, x_d[50:55], view(dxdt_d, 50:55), r, x_d[13], view(dxdt_d, 13))
	
	dxdt_d[13] = -r.kDeg*x_d[13] # Gas6 degradation
end

function internalSurfpYCalc(state)
	state[5]+
end





