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
    #=AXL::TAMrates{T}
	MerTK::TAMrates{T}
	Tyro3::TAMrates{T}=#
    TAMs = @SLVector TAMrates{T} (:Axl,:MerTK,:Tyro3)
    TAMs::TAMs
	
    internalize::T # Non-pY species internalization rate.
	pYinternalize::T # pY species internalization rate.
	fElse::T # Recycling fraction for non-D2 species.
	kRec::T # Recycling rate.
	kDeg::T # Degradation rate.
	xFwd::T
	internalFrac::T
	internalV::T
	gasCur::T
	
	#=AM::hetRates{T}
	AT::hetRates{T}
	MT::hetRates{T}=#
    hetR = @SLVector TAMrates{T} (:AM,:AT,:MT)
    hetR::hetR
end


" Setup the parameters for the full TAM receptor model. "
function param(params::Vector)
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


" Tracks the hetero-receptor interactions between two receptors. "
function het_module(Rone, Rtwo, dRone, dRtwo, hetR, hetDim, dhetDim, tr::Rates, Gas, dLi)
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

	return norm(dRr)
end


" Tracks the individual receptor-ligand interactions within a given receptor. "
function react_module(R, dR, dLi, Gas, r::TAMrates, tr::Rates)
	dRr = SVector{10}([r.binding[1]*R[1]*Gas - r.binding[2]*R[2],
                       r.binding[3]*R[1]*Gas - r.binding[4]*R[3],
                       r.binding[3]*R[2]*Gas - r.binding[4]*R[4],
                       r.binding[1]*R[3]*Gas - r.binding[2]*R[4],
                       tr.xFwd*R[1]*R[2] - r.xRev[1]*R[5],
                       tr.xFwd*R[1]*R[3] - r.xRev[2]*R[5],
                       tr.xFwd*R[1]*R[4] - r.xRev[3]*R[6],
                       tr.xFwd*R[2]*R[2] - r.xRev[4]*R[6],
                       tr.xFwd*R[3]*R[3] - r.xRev[5]*R[6],
                       r.xFwd6*Gas*R[5] - r.xRev[6]*R[6]])
    
	dR[1] += - dRr[7] - dRr[6] - dRr[5] - dRr[1] - dRr[2] # AXL
	dR[2] += -2*(dRr[8]) - dRr[5] + dRr[1] - dRr[3] # AXLgas1
	dR[3] += -2*(dRr[9]) - dRr[6] + dRr[2] - dRr[4] # AXLgas2
	dR[4] += -dRr[7] + dRr[3] + dRr[4] # AXLgas12
	dR[5] += -dRr[10] + dRr[6] + dRr[5] # AXLdimer1
	dR[6] += dRr[10] + dRr[9] + dRr[8] + dRr[7] # AXLdimer2

	if isnothing(dLi) == false
		dLi[1] += -dRr[10] - dRr[1] - dRr[2] - dRr[3] - dRr[4]
	end
    
    return norm(dRr)
end


" Handles trafficking of receptor and ligand. "
function trafFunc(dextR, dintR, intRate::Float64, extR, intR, kRec::Float64, kDeg::Float64, fElse::Float64, internalFrac::Float64)
	dextR[:] .+= -extR*intRate + kRec*(1-fElse)*intR*internalFrac # Endocytosis, recycling
	dintR[:] .+= extR*intRate/internalFrac - kRec*(1-fElse)*intR - kDeg*fElse*intR # Endocytosis, recycling, degradation
end


" Handles hetero-receptor interactions. "
function heteroTAM(Rone, Rtwo, dRone, dRtwo, hetR, hetDim, dhetDim, tr, Li, dLi)
	dnorm = het_module(Rone, Rtwo, dRone, dRtwo, hetR, hetDim, dhetDim, tr, tr.gasCur, nothing)
	dnorm += het_module(view(Rone, 7:10), view(Rtwo, 7:10), view(dRone, 7:10), view(dRtwo, 7:10), hetR, view(hetDim, 4:6), view(dhetDim, 4:6), tr, Li / tr.internalV, dLi);

	trafFunc(view(dhetDim, 1:3), view(dhetDim, 4:6), tr.pYinternalize, hetDim[1:3], hetDim[4:6], tr.kRec, tr.kDeg, 1.0, tr.internalFrac)

	return dnorm
end


function TAM_reactii(R, Li, dR, dLi, r::TAMrates, tr::Rates)
	dnorm = react_module(R, dR, nothing, tr.gasCur, r, tr)
	dnorm += react_module(view(R, 1:6), view(dR, 1:6), dLi, Li/tr.internalV, r, tr)

	dR[1] += r.expression

	trafFunc(view(dR, 1:4), view(dR, 7:10), tr.internalize, R[1:4], R[7:10], tr.kRec, tr.kDeg, tr.fElse, tr.internalFrac)
	trafFunc(view(dR, 4:5), view(dR, 10:11), tr.pYinternalize, R[4:5], R[10:11], tr.kRec, tr.kDeg, 1.0, tr.internalFrac)

	return dnorm
end


function TAM_reacti_dnorm(dxdt_d, x_d, params, t)
	fill!(dxdt_d, 0.0)
	r = param(params)

	dnorm = TAM_reactii(view(x_d, 1:12), x_d[13], view(dxdt_d, 1:12), view(dxdt_d, 13), r.AXL, r)
	dnorm += TAM_reactii(view(x_d, 14:25), x_d[13], view(dxdt_d, 14:25), view(dxdt_d, 13), r.MerTK, r)
	dnorm += TAM_reactii(view(x_d, 26:37), x_d[13], view(dxdt_d, 26:37), view(dxdt_d, 13), r.Tyro3, r)

	dnorm += heteroTAM(x_d[1:12],  x_d[14:25], view(dxdt_d, 1:12),  view(dxdt_d, 14:25), r.AM, x_d[38:43], view(dxdt_d, 38:43), r, x_d[13], view(dxdt_d, 13))
	dnorm += heteroTAM(x_d[14:25], x_d[26:37], view(dxdt_d, 14:25), view(dxdt_d, 26:37), r.MT, x_d[44:49], view(dxdt_d, 44:49), r, x_d[13], view(dxdt_d, 13))
	dnorm += heteroTAM(x_d[1:12],  x_d[26:37], view(dxdt_d, 1:12),  view(dxdt_d, 26:37), r.AT, x_d[50:55], view(dxdt_d, 50:55), r, x_d[13], view(dxdt_d, 13))

	dxdt_d[13] = -r.kDeg*x_d[13] # Gas6 degradation

	return dnorm
end


function TAM_reacti(dxdt_d, x_d, params, t)
	TAM_reacti_dnorm(dxdt_d, x_d, params, t)
end


function detailedBalance (out::Rates)
   
    for T in out.TAMs 
        KD1 = T.binding[2]/T.binding[1]
        KD2 = T.binding[3]/T.binding[2]
        
        if T.binding[2] <= T.binding[2]
            T.xRev[1] = T.binding[4]
            T.xRev[2] = T.xRev[1]*KD1/KD2
        else
            T.xRev[2] = T.binding[2]
            T.xRev[1] = T.xRev[2]*KD2/KD1 
        end
    end  
    
    for ii = 2:3
        out.TAMs[ii].xFwd6 = max(out.TAMs[ii].binding[1], out.TAMs[ii].binding[3])
        out.TAMs[ii].xRev[5] = out.TAMs[1].xRev[6]*out.TAMs[ii].binding[2]*out.TAMs[ii].binding[4]/out.TAMs[1].binding[4]/out.TAMs[1].binding[2];
    end
    
    for T in out.TAMs
        KD1 = T.binding[2]/T.binding[1]
        KD2 = T.binding[4]/T.binding[3]
            
        T.xRev[5] = T.xRev[6]*KD1*T.xRev[1]/KD2/KD2/T.xFwd6;
        T.xRev[4] = T.xRev[5]*KD2*KD2/KD1/KD1;
        T.xRev[3] = T.xRev[4]*KD1/KD2;
    end
    
    #=for ii = 1:3
            hetRates * const A = &out.hetR[ii];
            const unsigned int * const Ai = hetRi[ii];
            
            // ligand binding to the one ligand dimer
            KD11 = out.TAMs[Ai[0]].binding[1] / out.TAMs[hetRi[ii][0]].binding[0];
            KD12 = out.TAMs[Ai[1]].binding[1] / out.TAMs[hetRi[ii][1]].binding[0];
            KD21 = out.TAMs[Ai[0]].binding[3] / out.TAMs[hetRi[ii][0]].binding[2];
            KD22 = out.TAMs[Ai[1]].binding[3] / out.TAMs[hetRi[ii][1]].binding[2];
            
            A->xFwd15 = std::fmax(out.TAMs[Ai[1]].binding[0],out.TAMs[Ai[0]].binding[2]);
            A->xFwd16 = std::fmax(out.TAMs[Ai[0]].binding[0],out.TAMs[Ai[1]].binding[2]);
            
            A->xRev[0] = out.TAMs[hetRi[ii][0]].xRev[4]*KD12/KD11;
            A->xRev[1] = A->xRev[0]*KD21/KD12;
            A->xRev[2] = KD22*KD21/KD11/KD12*A->xRev[0];
            A->xRev[3] = A->xRev[0]*KD22/KD11;
            
            if (out.TAMs[Ai[0]].binding[1] <= out.TAMs[Ai[0]].binding[3]) {
                A->xRev[7] = out.TAMs[Ai[0]].binding[1];
                A->xRev[5] = A->xRev[7]*KD11/KD21;
            } else {
                A->xRev[5] = out.TAMs[Ai[0]].binding[3];
                A->xRev[7] = A->xRev[5]*KD21/KD11;
            }
            
            if (out.TAMs[Ai[1]].binding[1] <= out.TAMs[Ai[1]].binding[3]) {
                A->xRev[4] = out.TAMs[Ai[1]].binding[1];
                A->xRev[6] = A->xRev[4]*KD22/KD12;
            } else {
                A->xRev[6] = out.TAMs[Ai[1]].binding[3];
                A->xRev[4] = A->xRev[6]*KD12/KD22;
            }
            
            A->xRev[8] = A->xRev[0]*KD21/KD12;
            A->xRev[9] = A->xRev[2]*KD11/KD22;
        }
        
        return out;=#
end
