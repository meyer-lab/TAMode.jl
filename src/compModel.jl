
mutable struct comprates{T} #do we want this to be T and what does {T} mean
	rr::Rates #do i need values for all the rates and does this want to be in the reactCode file or a new one
	fraction::T #Fraction of cell surface covered with PtdSer
	partIn::T #Partitioning rate into PtdSer regions
	gasPart::T #Partitioning of ligand into PtdSer region
end


function react_module(R, dR, curL, r)
	dRr = SVector{30}([r.GBinding[1]*curL[1]*R[1] - r.GBinding[2]*R[2],
					   r.GBinding[3]*curL[1]*R[1] - r.GBinding[4]*R[3],
					   r.PBinding[1]*curL[2]*R[1] - r.PBinding[2]*R[4],
					   r.PBinding[3]*curL[2]*R[1] - r.PBinding[4]*R[5],
					   r.GBinding[3]*curL[1]*R[2] - r.GBinding[4]*R[6],
					   r.GBinding[1]*curL[1]*R[3] - r.GBinding[2]*R[6],
					   r.PBinding[1]*curL[2]*R[3] - r.PBinding[2]*R[7],
					   r.GBinding[3]*curL[1]*R[4] - r.GBinding[4]*R[7],
					   r.PBinding[3]*curL[2]*R[4] - r.PBinding[4]*R[8],
					   r.PBinding[1]*curL[2]*R[5] - r.PBinding[2]*R[8],
					   r.GBinding[1]*curL[1]*R[5] - r.GBinding[2]*R[9],
					   r.PBinding[3]*curL[2]*R[2] - r.PBinding[4]*R[9],

					   r.xFwd*R[2]*R[2] - r.xRev[1]*R[10],
					   r.xFwd*R[1]*R[6] - r.xRev[2]*R[10],
					   r.xFwd*R[3]*R[3] - r.xRev[3]*R[10],
					   r.xFwd*R[2]*R[1] - r.xRev[4]*R[11],
					   r.xFwd*R[3]*R[1] - r.xRev[5]*R[11],
					   r.xFwd*R[1]*R[5] - r.xRev[6]*R[12],
					   r.xFwd*R[1]*R[4] - r.xRev[7]*R[12],
					   r.xFwd*R[4]*R[4] - r.xRev[8]*R[13],
					   r.xFwd*R[1]*R[8] - r.xRev[9]*R[13],
					   r.xFwd*R[5]*R[5] - r.xRev[10]*R[13],
					   r.xFwd*R[7]*R[1] - r.xRev[11]*R[14],
					   r.xFwd*R[9]*R[1] - r.xRev[12]*R[14],
					   r.xFwd*R[2]*R[4] - r.xRev[13]*R[14],
					   r.xFwd*R[3]*R[5] - r.xRev[14]*R[14],
					   r.xFwd27*R[12]*curL[2] - r.xRev[15]*R[13],
					   r.xFwd27*R[11]*curL[2] - r.xRev[15]*R[14],
					   r.xFwd29*R[12]*curL[1] - r.xRev[16]*R[14],
					   r.xFwd29*R[11]*curL[1] - r.xRev[16]*R[10]])

	dR[1] = -sum(dRr[1, 2, 3, 4, 14, 16, 17, 18, 19, 21, 23, 24]) # A
	dR[1] = dRr[1] - sum(dRr[5, 12, 13, 13, 16, 25]) # AG (first site)
	dR[2] = dRr[2] - sum(dRr[6, 7, 15, 15, 17, 26]) # AG (second site)
	dR[3] = dRr[3] - sum(dRr[8, 9, 19, 20, 20, 25]) # AP (first site)
	dR[4] = dRr[4] - sum(dRr[10, 11, 18, 22, 22, 26]) # AP (second site)
	dR[5] = sumIdx(dRr.data(), {4, 5, -13}); // AGG
	dR[6] = sumIdx(dRr.data(), {6, 7, -22}); // APG
	dR[7] = sumIdx(dRr.data(), {8, 9, -20}); // APP
	dR[8] = sumIdx(dRr.data(), {10, 11, -23}); // AGP
	
	dR[9] = sumIdx(dRr.data(), {12, 13, 14, 29}); // AGGA
	dR[10] = sumIdx(dRr.data(), {15, 16, -27, -29}); // AGA
	dR[11] = sumIdx(dRr.data(), {17, 18, -26, -28}); // APA
	dR[12] = sumIdx(dRr.data(), {19, 20, 21, 26}); // APPA
	dR[13] = sumIdx(dRr.data(), {22, 23, 24, 25, 27, 28}); // APGA
	
	dR[14] = -sumIdx(dRr.data(), {0, 1, 4, 5, 7, 10, 28, 29}) - r->kDeg*R[14]; // Gasi
	dR[15] = -sumIdx(dRr.data(), {2, 3, 6, 8, 9, 11, 26, 27}) - r->kDeg*R[15]; // PROSi

	return norm(dRr)
end
