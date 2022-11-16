/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Martin Brehm.

    ---------------------------------------------------------------------------

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*****************************************************************************/


// This must always be the first include directive
#include "config.h"

#include "order.h"
#include "order_chain.h"
#include "order_vector.h"
#include "globalvar.h"
#include "largeinteger.h"


const char *GetRevisionInfo_order(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_order() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



void COrderEngine::REC_DetectChains(int depth, int pos, std::vector<int> &touched) {

	int z, z2, i, nc, nh;
	CMolBond *mb;
	CElement *el;
	std::vector<int> nbh;


	el = ((CAtom*)g_oaAtoms[g_waAtomRealElement[pos]])->m_pElement;

	if (mystricmp(el->m_sLabel,"C") != 0)
		return;

//	for (z=0;z<depth;z++)
//		mprintf("  ");
//	mprintf(">>> REC %d >>>\n",pos);

	nh = 0;
	nc = 0;

	for (z=0;z<m_pTempSM->m_oaBonds.GetSize();z++) {

		mb = (CMolBond*)m_pTempSM->m_oaBonds[z];

		if (mb->m_iAtomOffset[0] == pos)
			i = mb->m_iAtomOffset[1];
		else if (mb->m_iAtomOffset[1] == pos)
			i = mb->m_iAtomOffset[0];
		else
			continue;


		el = ((CAtom*)g_oaAtoms[g_waAtomRealElement[i]])->m_pElement;

		if (mystricmp(el->m_sLabel,"C") == 0) {
			nc++;
			if (touched[i] == 0)
				nbh.push_back(i);
			goto _accept;
		} else {
			for (z2=0;z2<(int)m_saChainAllowedSub.size();z2++) {
				if (mystricmp(el->m_sLabel,m_saChainAllowedSub[z2].c_str()) == 0) {
					nh++;
					goto _accept;
				}
			}
		}
		nc++;

_accept:;
	}

//	for (z=0;z<depth+1;z++)
//		mprintf("  ");
//	mprintf("nc=%d, nh=%d\n",nc,nh);

	// Only accept sp3 carbons
	if (nh + nc != 4)
		return;

	// Reject branched chains
	if (nc > 2)
		return;

	touched[pos] = 2;

	for (z=0;z<(int)nbh.size();z++)
		REC_DetectChains(depth+1,nbh[z],touched);

//	for (z=0;z<depth;z++)
//		mprintf("  ");
//	mprintf("<<< REC %d <<<\n",pos);
}


void COrderEngine::DetectChains() {

	int z, z2, z3, z4, z5, ic, n, ti, ti2, mai, tn, ntot;
	CMolecule *m;
	std::vector<int> tia, tlist;
	COrderChain *oc;
	LargeInteger li;
	CxString buf;
	char *p, *q, cbuf[16];


	mprintf("    Detecting alkyl chains...\n\n");

	m_iChainsTotal = 0;

	AskString("    Which elements to allow as carbon bond partners in alkyl chains? [H,D,F] ",&buf,"H,D,F");

	p = buf.GetWritePointer();

	while (*p != 0) {
		while (*p == ' ')
			p++;
		q = p;
		while ((*q != 0) && (*q != ',') && (*q != ' '))
			q++;
		memcpy(cbuf,p,q-p);
		cbuf[q-p] = 0;
		m_saChainAllowedSub.push_back((const char*)cbuf);
		while ((*q == ' ') || (*q == ','))
			q++;
		if (*q == 0)
			break;
		p = q;
	}
	mprintf("\n");

	//for (z=0;z<(int)m_saChainAllowedSub.size();z++)
	//	mprintf("\"%s\"\n",m_saChainAllowedSub[z].c_str());

	tia.resize(g_iGesAtomCount);

	m_oaChains.resize(g_oaMolecules.GetSize());

	ntot = 0;

	for (z=0;z<g_oaMolecules.GetSize();z++) {

		m = (CMolecule*)g_oaMolecules[z];
		m_pTempSM = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];

		mprintf(WHITE,"      %s ...\n",m->m_sName);

		ic = -1;
		for (z2=0;z2<m->m_baAtomIndex.GetSize()-1;z2++) {
			if (mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_pElement->m_sLabel,"C") == 0) {
				ic = z2;
				break;
			}
		}
		if (ic == -1) {
			mprintf("        Does not contain carbon atoms.\n\n");
			continue;
		}

		for (z2=0;z2<g_iGesAtomCount;z2++)
			tia[z2] = 0;

		for (z2=0;z2<g_iGesAtomCount;z2++) {

			if (g_laAtomSMIndex[z2] != m->m_laSingleMolIndex[0])
				continue;

			if (tia[z2] != 0)
				continue;

			REC_DetectChains(0,z2,tia);

			n = 0;
			tlist.clear();
			for (z3=0;z3<g_iGesAtomCount;z3++) {

				if (tia[z3] == 2) {
	/*	//			if (b)
		//				mprintf(", ");
		//			else
		//				mprintf("        ");
					for (z4=0;z4<((CxIntArray*)m_pTempSM->m_oaAtomOffset[ic])->GetSize();z4++) {
						if (((CxIntArray*)m_pTempSM->m_oaAtomOffset[ic])->GetAt(z4) == z3) {
		//					mprintf("C%d",z4+1);
							goto _done;
						}
					}
					eprintf("Internal error: Not a carbon atom.\n");
					abort();
_done:*/
					//b = true;
					n++;
					tia[z3] = 1;
					tlist.push_back(z3);
				}
			}
		//	if (b)
		//		mprintf(" (%d)\n",n);

			if (tlist.size() < 2) // Process further if at least an ethyl group
				continue;

			oc = new COrderChain();
			m_iChainsTotal++;
			oc->m_iMolecule = z;
			oc->m_iCarbonIndex = ic;
			m_oaChains[z].push_back(oc);
			ntot++;

			ti = -1;
			tn = -1;
			for (z3=0;z3<(int)tlist.size();z3++) {

				n = 0;
				for (z4=0;z4<m_pTempSM->m_oaBonds.GetSize();z4++) {
					if (((CMolBond*)m_pTempSM->m_oaBonds[z4])->m_iAtomOffset[0] == tlist[z3]) {
						for (z5=0;z5<(int)tlist.size();z5++) {
							if (((CMolBond*)m_pTempSM->m_oaBonds[z4])->m_iAtomOffset[1] == tlist[z5]) {
								n++;
								goto _bf;
							}
						}
					} else if (((CMolBond*)m_pTempSM->m_oaBonds[z4])->m_iAtomOffset[1] == tlist[z3]) {
						for (z5=0;z5<(int)tlist.size();z5++) {
							if (((CMolBond*)m_pTempSM->m_oaBonds[z4])->m_iAtomOffset[0] == tlist[z5]) {
								n++;
								goto _bf;
							}
						}
					}
_bf:;
				}

				if (n == 1) {
					mai = -1;
					for (z4=0;z4<m_pTempSM->m_oaMolAtoms.GetSize();z4++) {
						if (((CMolAtom*)m_pTempSM->m_oaMolAtoms[z4])->m_iOffset == tlist[z3]) {
							mai = z4;
							break;
						}
					}
					if (mai == -1) {
						eprintf("Internal error: CMolAtom not found.\n");
						abort();
					}
					if (ti != -1) {
						if (((CMolAtom*)m_pTempSM->m_oaMolAtoms[mai])->m_liAtomCode > li) {
							ti = tlist[z3];
							li = ((CMolAtom*)m_pTempSM->m_oaMolAtoms[mai])->m_liAtomCode;
							tn = ((CMolAtom*)m_pTempSM->m_oaMolAtoms[mai])->m_iNumber;
						}
					} else {
						ti = tlist[z3];
						li = ((CMolAtom*)m_pTempSM->m_oaMolAtoms[mai])->m_liAtomCode;
						tn = ((CMolAtom*)m_pTempSM->m_oaMolAtoms[mai])->m_iNumber;
					}
				} else if (n != 2) {
					eprintf("Internal error: Chain atom has > 2 bonds: %d.\n",n);
					abort();
				}
			}

			//mprintf("Start is C%d.\n",tn+1);

			oc->m_iaAtoms.push_back(tn);

			ti2 = -1;
			for (z3=0;z3<(int)tlist.size()-1;z3++) {

				for (z4=0;z4<m_pTempSM->m_oaBonds.GetSize();z4++) {
					if (((CMolBond*)m_pTempSM->m_oaBonds[z4])->m_iAtomOffset[0] == ti) {
						for (z5=0;z5<(int)tlist.size();z5++) {
							if (((CMolBond*)m_pTempSM->m_oaBonds[z4])->m_iAtomOffset[1] == tlist[z5]) {
								if (tlist[z5] != ti2) {
									ti2 = ti;
									ti = tlist[z5];
									goto _cd;
								}
							}
						}
					} else if (((CMolBond*)m_pTempSM->m_oaBonds[z4])->m_iAtomOffset[1] == ti) {
						for (z5=0;z5<(int)tlist.size();z5++) {
							if (((CMolBond*)m_pTempSM->m_oaBonds[z4])->m_iAtomOffset[0] == tlist[z5]) {
								if (tlist[z5] != ti2) {
									ti2 = ti;
									ti = tlist[z5];
									goto _cd;
								}
							}
						}
					}
				}
_cd:
				for (z4=0;z4<m_pTempSM->m_oaMolAtoms.GetSize();z4++) {
					if (((CMolAtom*)m_pTempSM->m_oaMolAtoms[z4])->m_iOffset == ti) {
						oc->m_iaAtoms.push_back(((CMolAtom*)m_pTempSM->m_oaMolAtoms[z4])->m_iNumber);
						goto _cd2;
					}
				}
				eprintf("Internal error: CMolAtom not found (2).\n");
				abort();
_cd2:;
				//mprintf(" C%d",oc->m_iaAtoms[oc->m_iaAtoms.size()-1]+1);
			}
			//mprintf("\n");
		}

		// Sort chains descending by length
		for (z2=0;z2<(int)m_oaChains[z].size()-1;z2++) {
			ti = -1;
			ti2 = -1;
			for (z3=z2;z3<(int)m_oaChains[z].size();z3++) {
				if ((int)m_oaChains[z][z3]->m_iaAtoms.size() > ti) {
					ti = (int)m_oaChains[z][z3]->m_iaAtoms.size();
					ti2 = z3;
				}
			}
			if (ti2 != z2) {
				oc = m_oaChains[z][z2];
				m_oaChains[z][z2] = m_oaChains[z][ti2];
				m_oaChains[z][ti2] = oc;
			}
		}

		for (z2=0;z2<(int)m_oaChains[z].size();z2++) {
			oc = (COrderChain*)m_oaChains[z][z2];
			mprintf("       %2d.) ",z2+1);
			for (z3=0;z3<(int)oc->m_iaAtoms.size();z3++) {
				if (z3 != 0)
					mprintf("-");
				mprintf("C%d",oc->m_iaAtoms[z3]+1);
			}
			mprintf("  (%lu)\n",(unsigned long)oc->m_iaAtoms.size());
		}

		mprintf("\n");
	}

	mprintf("    Chain detection completed. Detected %d alkyl chains.\n\n\n",ntot);
}


bool COrderEngine::Parse(CTimeStep *ts) {

	COrderAnalysis *v;
	CxString buf;

	mprintf(WHITE,"\n    >>> Order Parameters >>>\n\n");

	DetectChains();

	while (true) {

		mprintf(WHITE,"    *** Analysis %lu ***\n\n",(unsigned long)m_oaAnalyses.size()+1);

		mprintf("    Available analyses:\n");
		mprintf("      vec    - Vector-based Order Parameters\n");
		mprintf("      chain  - Alkyl Chain Order Parameters\n");
		mprintf("\n");

_again:
		AskString_ND("    Which analysis to perform? ",&buf);

		if (strchr((const char*)buf,',') != NULL) {
			eprintf("\n    Please enter only one analysis at a time.\n\n");
			goto _again;
		}

		if (mystricmp((const char*)buf,"vec") == 0)
			v = new COrderAnalysis_Vector(this);
		else if (mystricmp((const char*)buf,"chain") == 0)
			v = new COrderAnalysis_Chain(this);
		else {
			eprintf("\n    Unknown analysis: \"%s\".\n\n",(const char*)buf);
			goto _again;
		}

		v->m_iNumber = (int)m_oaAnalyses.size();
		m_oaAnalyses.push_back(v);

		if (!v->Parse(ts))
			return false;

		if (!AskYesNo("    Add another order parameter analysis (y/n)? [no] ",false))
			break;

		mprintf("\n");
	}

	mprintf("\n    Added %lu analyses.\n",(unsigned long)m_oaAnalyses.size());

	mprintf(WHITE,"\n    <<< Order Parameters <<<\n\n");

	return true;
}


void COrderEngine::ProcessStep(CTimeStep *ts) {

	int z;

	for (z=0;z<(int)m_oaAnalyses.size();z++)
		m_oaAnalyses[z]->ProcessStep(ts);
}


void COrderEngine::Finish() {

	int z;

	mprintf(WHITE,"\n    Finishing Order Parameter Analyses...\n");

	for (z=0;z<(int)m_oaAnalyses.size();z++)
		m_oaAnalyses[z]->Finish();
}



