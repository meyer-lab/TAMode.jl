

function het_module(Rone, Rtwo, dRone, dRtwo, hetR, hetDim, dhetDim, tr, Gas, dLi)
	dRr = zeros(eltype(Rone), 10)

	dRr[1] = tr.xFwd*Rone[3]*Rtwo[3] - hetR.xRev[1]*hetDim[3]
    dRr[2] = tr.xFwd*Rone[1]*Rtwo[4] - hetR.xRev[2]*hetDim[3]
    dRr[3] = tr.xFwd*Rone[2]*Rtwo[2] - hetR.xRev[3]*hetDim[3]
    dRr[4] = tr.xFwd*Rone[4]*Rtwo[1] - hetR.xRev[4]*hetDim[3]
    dRr[5] = tr.xFwd*Rone[3]*Rtwo[1] - hetR.xRev[5]*hetDim[2]
    dRr[6] = tr.xFwd*Rone[1]*Rtwo[2] - hetR.xRev[6]*hetDim[2]
    dRr[7] = tr.xFwd*Rone[2]*Rtwo[1] - hetR.xRev[7]*hetDim[1]
    dRr[8] = tr.xFwd*Rone[1]*Rtwo[3] - hetR.xRev[8]*hetDim[1]
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
    
    dLi[1] += -dRr[9] - dRr[10]

	return norm(dRr)
end