/***********************************************************************************

    LibBQB - File Format and Compression Algorithms for Trajectories of
             Volumetric Data and Atom Positions

    https://brehm-research.de/bqb

    Free software, licensed under GNU LGPL v3

    Copyright (c) Martin Brehm and Martin Thomas,
                  Martin Luther University Halle-Wittenberg, Germany,
                  2016 - 2021.

    Please cite:  M. Brehm, M. Thomas: "An Efficient Lossless Compression Algorithm
                  for Trajectories of Atom Positions and Volumetric Data",
                  J. Chem. Inf. Model. 2018, 58 (10), pp 2092-2107.

    --------------------------------------------------------------------------------

    LibBQB is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

***********************************************************************************/



// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_extrapolator.h"
#include "bqb_math.h"
#include "bqb_linalg.h"
#include "bqb_interface.h"
#include <set>


const char *GetRevisionInfo_bqb_extrapolator(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_extrapolator() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



void CBQBExtrapolator::Reset() {

	int z;

	m_iCubeFramePos = 0;
	m_iCubeFrameCount = 0;
	m_iCubeFrameVisibleCount = 0;
	for (z=0;z<(int)m_oaCubeFrames.size();z++)
		m_oaCubeFrames[z] = NULL;
	for (z=0;z<(int)m_oaAtomFrames.size();z++)
		m_oaAtomFrames[z] = NULL;
}


void CBQBExtrapolator::InitializeXYZ(CBQBParameterSet_Position *p, bool silent) {

	InitializeXYZ(
		p->GetExtraTRange(),
		p->GetExtraTOrder(),
		p->GetExtraTimeExpo(),
		silent
	);
}


const CBQBCubeFrame* CBQBExtrapolator::GetCubeFrame(int depth) const {

	int i;

	if (depth >= m_iCubeFrameVisibleCount) {
		m_IF.eprintf("CBQBExtrapolator::GetCubeFrame(): Error: %d >= %d (%d).\n",
			depth,m_iCubeFrameVisibleCount,m_iCubeFrameCount);
		abort();
	}
	i = m_iCubeFramePos - depth;
	if (i < 0)
		i += m_iTRange;
	return m_oaCubeFrames[i];
}


const CBQBAtomSet* CBQBExtrapolator::GetAtomFrame(int depth) const {

	int i;

	if (depth >= m_iCubeFrameVisibleCount) {
		m_IF.eprintf("CBQBExtrapolator::GetAtomFrame(): Error: %d >= %d (%d).\n",
			depth,m_iCubeFrameVisibleCount,m_iCubeFrameCount);
		abort();
	}
	i = m_iCubeFramePos - depth;
	if (i < 0)
		i += m_iTRange;
/*	if (((void*)m_oaAtomFrames[i]) < (void*)0xFFFF) {
		m_IF.printf("@ GetAtomFrame(%d): %d, %08X\n",depth,i,m_oaAtomFrames[i]);
		m_IF.printf("Pos %d, Count %d, VisCo %d.\n",m_iCubeFramePos,m_iCubeFrameCount,m_iCubeFrameVisibleCount);
	}*/
	return m_oaAtomFrames[i];
}


double CBQBExtrapolator::ExtrapolateKnown(int pattern, int index) const {

	double tf;
	unsigned int z;

	const CBQBExtraPattern &p = m_oaPatterns[pattern];
	tf = 0;
	for (z=0;z<p.m_iaSpatialIndex.size();z++) {

		const int ti = index + p.m_iaSpatialIndex[z];

		#ifdef CHECK_EXTRAPOLATOR
			if (ti < 0)
				m_IF.eprintf("CBQBExtrapolator::ExtrapolateKnown(): Pattern %d, Term %u/%lu: Fail1 %d + %d = %d\n",
					pattern,z,(unsigned long)p.m_iaSpatialIndex.size(),index,p.m_iaSpatialIndex[z],ti);
			if (ti >= m_iRes012/*m_iRes[0]*m_iRes[1]*m_iRes[2]*/)
				m_IF.eprintf("CBQBExtrapolator::ExtrapolateKnown(): Pattern %d, Term %u/%lu: Fail2 %d + %d = %d\n",
					pattern,z,(unsigned long)p.m_iaSpatialIndex.size(),index,p.m_iaSpatialIndex[z],ti);
			if ((p.m_iaSpatialIndex[z] >= 0) && (p.m_iaTemporalIndex[z] == 0))
				m_IF.eprintf("CBQBExtrapolator::ExtrapolateKnown(): Pattern %d, Term %u/%lu: Fail3.\n",
					pattern,z,(unsigned long)p.m_iaSpatialIndex.size());
		#endif

		const CBQBCubeFrame *cfr = GetCubeFrame(p.m_iaTemporalIndex[z]);
	//	tf += p.m_faCoeff[z] * cfr->m_faBin[index+p.m_iaSpatialIndex[z]];
	//	tf += p.m_faCoeff[z] * cfr->m_iaMantis[ti] * pow10(cfr->m_iaExpo[ti]);
		tf += p.m_faCoeff[z] * cfr->m_faBin[ti];
	}
	return tf;
}


double CBQBExtrapolator::ExtrapolateKnownCorr(const std::vector<double> &fa, int pattern, int index) const {

	double tf;
	unsigned int z;

	const CBQBExtraPattern &p = m_oaPatterns[pattern];
	tf = 0;
	for (z=0;z<p.m_iaSpatialIndex.size();z++) {
		const int ti = index + p.m_iaSpatialIndex[z];

		#ifdef CHECK_EXTRAPOLATOR
			if (ti < 0)
				m_IF.eprintf("CBQBExtrapolator::ExtrapolateKnownCorr(): Pattern %d, Term %u/%lu: Fail1 %d + %d = %d\n",
					pattern,z,(unsigned long)p.m_iaSpatialIndex.size(),index,p.m_iaSpatialIndex[z],ti);
			if (ti >= m_iRes012)
				m_IF.eprintf("CBQBExtrapolator::ExtrapolateKnownCorr(): Pattern %d, Term %u/%lu: Fail2 %d + %d = %d\n",
					pattern,z,(unsigned long)p.m_iaSpatialIndex.size(),index,p.m_iaSpatialIndex[z],ti);
			if (p.m_iaSpatialIndex[z] >= 0)
				m_IF.eprintf("CBQBExtrapolator::ExtrapolateKnownCorr(): Pattern %d, Term %u/%lu: Fail3.\n",
					pattern,z,(unsigned long)p.m_iaSpatialIndex.size());
		#endif

		tf += p.m_faCoeff[z] * fa[ti];
	}
	return tf;
}

	
void CBQBExtrapolator::Initialize(
		int resx,
		int resy,
		int resz,
		int srangex,
		int srangey,
		int srangez,
		int trange,
		int sorder,
		int torder,
		int offsetx,
		int offsety,
		int offsetz,
		bool crosss,
		bool crosst,
		bool wrap,
		bool crossranges,
		bool crossranget,
		double distexpo,
		double timeexpo,
		bool silent
	) {

	int ix, iy, iz, it, rx, ry, rz, rt, dp, i, j2, morder;
	int jx, jy, jz, jt, kx, ky, kz, kt;
	int dx, dy, dz, dt;
	int maxdp, maxcoeff;
	double tf;
	CBQBDMatrixMN minp, mout;
	std::vector<int> iax, iay, iaz, iat, iaind;
	std::set<int> setx, sety, setz, sett;
	std::set<int>::iterator setit;
	CBQBLinAlg la(m_IF);
	CBQBMath bqbmath(m_IF);


	if ((offsetx < 0) || (offsetx >= srangex)) {
		m_IF.eprintf("CBQBExtrapolator::Initialize(): Error: Invalid value for offsetx (found %d, allowed 0 .. %d).\n",
			offsetx,srangex-1);
		abort();
	}
	if ((offsety < 0) || (offsety >= srangey)) {
		m_IF.eprintf("CBQBExtrapolator::Initialize(): Error: Invalid value for offsety (found %d, allowed 0 .. %d).\n",
			offsety,srangey-1);
		abort();
	}
	if ((offsetz < 0) || (offsetz >= srangez)) {
		m_IF.eprintf("CBQBExtrapolator::Initialize(): Error: Invalid value for offsetz (found %d, allowed 0 .. %d).\n",
			offsetz,srangez-1);
		abort();
	}

	if (trange < 1) {
		m_IF.eprintf("CBQBExtrapolator::Initialize(): Error: Parameter trange needs to be >= 1.\n");
		abort();
	}

	m_fDistExpoAdd = 1.0;
	m_fTimeExpoAdd = 1.0;

	if (sorder < 0) { // Order -1 means no extrapolation, which is less than order 0 (constant extrapolation)
		sorder = 0;
		srangex = 1;
		srangey = 1;
		srangez = 1;
	}
	
	if (torder < 0) {
		torder = 0;
		trange = 1;
	}

	m_iSRange[0] = srangex;
	m_iSRange[1] = srangey;
	m_iSRange[2] = srangez;
	m_iTRange = trange;
	m_iSOffset[0] = offsetx;
	m_iSOffset[1] = offsety;
	m_iSOffset[2] = offsetz;
	m_iSOrder = sorder;
	m_iTOrder = torder;

	m_bCrossS = crosss;
	m_bCrossRangeS = crossranges;
	m_bWrap = wrap;
	m_bCrossT = crosst;
	m_bCrossRangeT = crossranget;
	m_fDistExpo = distexpo;
	m_fTimeExpo = timeexpo;

	m_iRes[0] = resx;
	m_iRes[1] = resy;
	m_iRes[2] = resz;

	m_iSRange01 = m_iSRange[0] * m_iSRange[1];
	m_iSRange012 = m_iSRange[0] * m_iSRange[1] * m_iSRange[2];
	m_iTempVal[0] = m_iRes[0] - m_iSRange[0] + m_iSOffset[0];
	m_iTempVal[1] = m_iRes[1] - m_iSRange[1] + m_iSOffset[1];
	m_iTempVal[2] = m_iRes[2] - m_iSRange[2] + m_iSOffset[2];
	m_iRes012 = m_iRes[0] * m_iRes[1] * m_iRes[2];
	m_iRes12 = m_iRes[1] * m_iRes[2];

	maxdp = 0;
	maxcoeff = 0;


	m_oaCubeFrames.resize(m_iTRange);
	m_oaAtomFrames.resize(m_iTRange);
	m_iCubeFramePos = 0;
	m_iCubeFrameCount = 0;
	m_iCubeFrameVisibleCount = 0;
	for (ix=0;ix<m_iTRange;ix++) {
		m_oaCubeFrames[ix] = NULL;
		m_oaAtomFrames[ix] = NULL;
	}

	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("Initializing Extrapolator...\n");

	for (it=0;it<m_iTRange;it++) {

		for (iz=0;iz<m_iSRange[2];iz++) {

			for (iy=0;iy<m_iSRange[1];iy++) {

				for (ix=0;ix<m_iSRange[0];ix++) {

					if (m_IF.IsPL(BQB_PL_DEBUG)) {
						m_IF.printf("\n");
						m_IF.printf("*** Pattern %lu: Preparing Matrix for t=%d x=%d y=%d z=%d ...\n",
							(unsigned long)m_oaPatterns.size(),it,ix,iy,iz);
					}

					m_oaPatterns.push_back(CBQBExtraPattern());

					dp = 0;

					iax.clear();
					iay.clear();
					iaz.clear();
					iat.clear();
					iaind.clear();

					setx.clear();
					sety.clear();
					setz.clear();
					sett.clear();

					for (jt=0;jt<m_iTRange-it;jt++) {

						for (jx=-m_iSOffset[0];jx<m_iSRange[0]-m_iSOffset[0];jx++) {

							if (!wrap) {
								if (ix > m_iSOffset[0]) {
									if (jx >= m_iSRange[0]-ix)
										continue;
								} else /*if (ix <= m_iSOffset[0])*/ {
									if (jx < -m_iSOffset[0]+ix)
										continue;
								}
							}

							if ((ix <= m_iSOffset[0]) && (jx < -m_iSOffset[0]+ix) && wrap)
								kx = jx+resx;
							else if ((ix > m_iSOffset[0]) && (jx > m_iSRange[0]-1-ix) && wrap)
								kx = jx-resx;
							else
								kx = jx;

							for (jy=-m_iSOffset[1];jy<m_iSRange[1]-m_iSOffset[1];jy++) {

								if (!wrap) {
									if (iy > m_iSOffset[1]) {
										if (jy >= m_iSRange[1]-iy)
											continue;
									} else /*if (iy <= m_iSOffset[1])*/ {
										if (jy < -m_iSOffset[1]+iy)
											continue;
									}
								}

								if ((iy <= m_iSOffset[1]) && (jy < -m_iSOffset[1]+iy) && wrap)
									ky = jy+resy;
								else if ((iy > m_iSOffset[1]) && (jy > m_iSRange[1]-1-iy) && wrap)
									ky = jy-resy;
								else
									ky = jy;

								for (jz=-m_iSOffset[2];jz<m_iSRange[2]-m_iSOffset[2];jz++) {

									if (!wrap) {
										if (iz > m_iSOffset[2]) {
											if (jz >= m_iSRange[2]-iz)
												continue;
										} else /*if (iz <= m_iSOffset[2])*/ {
											if (jz < -m_iSOffset[2]+iz)
												continue;
										}
									}

									if ((iz <= m_iSOffset[2]) && (jz < -m_iSOffset[2]+iz) && wrap)
										kz = jz+resz;
									else if ((iz > m_iSOffset[2]) && (jz > m_iSRange[2]-1-iz) && wrap)
										kz = jz-resz;
									else
										kz = jz;

									if (!crossranges && (((jx != 0) && (jy != 0)) || ((jx != 0) && (jz != 0)) || ((jy != 0) && (jz != 0))))
										continue;

									if (!crossranget && (jt != 0) && ((jx != 0) || (jy != 0) || (jz != 0)))
										continue;

									if ((jt == 0) && ((long)kx*resy*resz+(long)ky*resz+(long)kz >= 0))
										continue;

									if (m_IF.IsPL(BQB_PL_DEBUG)) {
										m_IF.printf("%3d:  ( ",dp+1);
										m_IF.printf("%2d",jt);
										m_IF.printf(" | ");
										m_IF.printf("%2d:%4d",jx,kx);
										m_IF.printf(" | ");
										m_IF.printf("%2d:%4d",jy,ky);
										m_IF.printf(" | ");
										m_IF.printf("%2d:%4d",jz,kz);
										m_IF.printf(" ) ");
										m_IF.printf("DIdx %10ld",(long)kx*resy*resz+(long)ky*resz+(long)kz);
										m_IF.printf("\n");
									}

									iax.push_back(jx);
									iay.push_back(jy);
									iaz.push_back(jz);
									iat.push_back(jt);
									setx.insert(jx);
									sety.insert(jy);
									setz.insert(jz);
									sett.insert(jt);
									iaind.push_back((long)kx*resy*resz+(long)ky*resz+(long)kz);
									dp++;
								}
							}
						}
					}

					if (dp == 0) {
						if (m_IF.IsPL(BQB_PL_DEBUG))
							m_IF.printf("  Skipping.\n");
						continue;
					}


					rx = (int)setx.size();
					ry = (int)sety.size();
					rz = (int)setz.size();
					rt = (int)sett.size();

					if (m_IF.IsPL(BQB_PL_DEBUG)) {

						m_IF.printf("  %lu different X values:",(unsigned long)setx.size());
						for (setit=setx.begin();setit!=setx.end();++setit)
							m_IF.printf(" %d",*setit);
						m_IF.printf("\n");

						m_IF.printf("  %lu different Y values:",(unsigned long)sety.size());
						for (setit=sety.begin();setit!=sety.end();++setit)
							m_IF.printf(" %d",*setit);
						m_IF.printf("\n");

						m_IF.printf("  %lu different Z values:",(unsigned long)setz.size());
						for (setit=setz.begin();setit!=setz.end();++setit)
							m_IF.printf(" %d",*setit);
						m_IF.printf("\n");

						m_IF.printf("  %lu different T values:",(unsigned long)sett.size());
						for (setit=sett.begin();setit!=sett.end();++setit)
							m_IF.printf(" %d",*setit);
						m_IF.printf("\n");
					}

					morder = MIN(sorder,MAX3(rx-1,ry-1,rz-1));
_norder:
					if (rx-1 < morder) dx = rx-1; else dx = morder;
					if (ry-1 < morder) dy = ry-1; else dy = morder;
					if (rz-1 < morder) dz = rz-1; else dz = morder;
					if (rt-1 < torder) dt = rt-1; else dt = torder;
					if (m_IF.IsPL(BQB_PL_DEBUG))
						m_IF.printf("  Using MaxOrder=%d, degree X=%d, Y=%d, Z=%d, T=%d\n",morder,dx,dy,dz,dt);
					i = 0;
					for (jx=0;jx<=dx;jx++) {
						for (jy=0;jy<=dy;jy++) {
							if (jx+jy > morder)
								break;
							for (jz=0;jz<=dz;jz++) {
								if (jx+jy+jz > morder)
									break;
								if (!crosss && (((jx != 0) && (jy != 0)) || ((jx != 0) && (jz != 0)) || ((jy != 0) && (jz != 0))))
									continue;
								for (jt=0;jt<=dt;jt++) {
									if (!crosst && (jt != 0) && ((jx != 0) || (jy != 0) || (jz != 0)))
										break;
									i++;
								}
							}
						}
					}
					if (i > dp) {
						morder--;
						goto _norder;
					}

					if (m_IF.IsPL(BQB_PL_DEBUG)) {

						m_IF.printf("  These are %d coefficients.\n",i);
						m_IF.printf("  Constructing %dx%d matrix...\n",dp,i);
					}

					if (i > maxcoeff)
						maxcoeff = i;

					minp = CBQBDMatrixMN(dp,i);

					if (m_IF.IsPL(BQB_PL_DEBUG))
						m_IF.printf("    ");
					i = 0;
					for (jt=0;jt<=dt;jt++) {
						for (jx=0;jx<=dx;jx++) {
							for (jy=0;jy<=dy;jy++) {
								if (jx+jy > morder)
									break;
								for (jz=0;jz<=dz;jz++) {
									if (jx+jy+jz > morder)
										break;
									if (!crosss && (((jx != 0) && (jy != 0)) || ((jx != 0) && (jz != 0)) || ((jy != 0) && (jz != 0))))
										continue;
									if (!crosst && (jt != 0) && ((jx != 0) || (jy != 0) || (jz != 0)))
										break;

									if (m_IF.IsPL(BQB_PL_DEBUG)) {

										if ((jx==0) && (jy==0) && (jz==0) && (jt==0))
											m_IF.printf("  1");
										if (jx > 1)
											m_IF.printf("  x^%d",jx);
										else if (jx == 1)
											m_IF.printf("  x");
										if (jy > 1)
											m_IF.printf("  y^%d",jy);
										else if (jy == 1)
											m_IF.printf("  y");
										if (jz > 1)
											m_IF.printf("  z^%d",jz);
										else if (jz == 1)
											m_IF.printf("  z");
										if (jt > 1)
											m_IF.printf("  t^%d",jt);
										else if (jt == 1)
											m_IF.printf("  t");
									}

									for (kt=0;kt<dp;kt++)
										minp(kt,i) = (double)(pow(iat[kt],jt) * pow(iax[kt],jx) * pow(iay[kt],jy) * pow(iaz[kt],jz)) *
													1.0/pow((double)iax[kt]*iax[kt]+iay[kt]*iay[kt]+iaz[kt]*iaz[kt]+m_fDistExpoAdd,0.5*distexpo) *
													1.0/pow((double)iat[kt]+m_fTimeExpoAdd,timeexpo);

									i++;
								}
							}
						}
					}
					if (m_IF.IsPL(BQB_PL_DEBUG))
						m_IF.printf("\n");

					mout = la.BQB_ComputePseudoInverse(minp);

					j2 = 0;
					for (kt=0;kt<dp;kt++)
						if (fabs(mout(0,kt)) > 1.0E-13)
							j2++;

					if (m_IF.IsPL(BQB_PL_DEBUG))
						m_IF.printf("  Have %d zero coefficients, %d data points remain.\n",dp-j2,j2);

					if (j2 > maxdp)
						maxdp = j2;

					m_oaPatterns.back().m_iaSpatialIndex.resize(j2);
					m_oaPatterns.back().m_iaTemporalIndex.resize(j2);
					m_oaPatterns.back().m_faCoeff.resize(j2);

					if (m_IF.IsPL(BQB_PL_DEBUG)) {
						m_IF.printf("  Matrix:\n");
						m_IF.printf("  ");
					}
					i = 0;
					for (kt=0;kt<dp;kt++) {
						if (fabs(mout(0,kt)) <= 1.0E-13)
							continue;
						m_oaPatterns.back().m_iaSpatialIndex[i] = iaind[kt];
						m_oaPatterns.back().m_iaTemporalIndex[i] = iat[kt];
						if (m_IF.IsPL(BQB_PL_DEBUG))
							m_IF.printf("%2d|%2d|%2d|%2d  ",iat[kt],iax[kt],iay[kt],iaz[kt]);
						i++;
					}
					if (m_IF.IsPL(BQB_PL_DEBUG))
						m_IF.printf("\n");

					if (m_IF.IsPL(BQB_PL_DEBUG)) {

						if ((minp.GetCols() < 20) && (minp.GetRows() < 20))
							bqbmath.DumpMatrix(minp);
						else
							m_IF.printf("[too large]\n");

						m_IF.printf("  Pseudo Inverse:\n");

						if ((minp.GetCols() < 20) && (minp.GetRows() < 20))
							bqbmath.DumpMatrix(mout);
						else
							m_IF.printf("[too large]\n");
					}

					i = 0;
					tf = 0;
					for (kt=0;kt<dp;kt++) {
						if (fabs(mout(0,kt)) <= 1.0E-13)
							continue;
						m_oaPatterns.back().m_faCoeff[i] = mout(0,kt) *
							1.0/pow((double)iax[kt]*iax[kt]+iay[kt]*iay[kt]+iaz[kt]*iaz[kt]+m_fDistExpoAdd,0.5*distexpo) *
							1.0/pow((double)iat[kt]+m_fTimeExpoAdd,timeexpo);
						tf += m_oaPatterns.back().m_faCoeff[i];

				//		m_oaPatterns.back().m_faCoeff[i] = floor(m_oaPatterns.back().m_faCoeff[i]*10000.0+0.5) / 10000.0;

						i++;
					}
					if (m_IF.IsPL(BQB_PL_DEBUG))
						m_IF.printf("  Vector: ");
					for (kt=0;kt<(int)m_oaPatterns.back().m_faCoeff.size();kt++) {
				//		m_oaPatterns.back().m_faCoeff[kt] *= corrfactor / tf;
						m_oaPatterns.back().m_faCoeff[kt] /= tf;
						if (m_IF.IsPL(BQB_PL_DEBUG))
							m_IF.printf("  %.16f",m_oaPatterns.back().m_faCoeff[kt]);
					}
					if (m_IF.IsPL(BQB_PL_DEBUG))
						m_IF.printf("\n");

				}
			}
		}
	}

	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		m_IF.printf("Extrapolator initialized.\n");
		m_IF.printf("%d data points, %d coefficients, %lu patterns.\n",maxdp,maxcoeff,(unsigned long)m_oaPatterns.size());
		m_IF.printf("\n");
	} else if (!silent)
		m_IF.printf("        %d data points, %d coefficients, %lu patterns.\n",maxdp,maxcoeff,(unsigned long)m_oaPatterns.size());
}


void CBQBExtrapolator::PushCubeFrame(const CBQBCubeFrame *frame) {

	if (m_iCubeFrameCount < m_iTRange)
		m_iCubeFrameCount++;
	m_iCubeFrameVisibleCount = m_iCubeFrameCount;

	m_iCubeFramePos++;
	if (m_iCubeFramePos >= m_iTRange)
		m_iCubeFramePos -= m_iTRange;

	m_oaCubeFrames[m_iCubeFramePos] = frame;
}


void CBQBExtrapolator::PushAtomFrame(const CBQBAtomSet *frame) {

	//m_IF.printf("PushAtomFrame before: count %d, visco %d, pos %d\n",m_iCubeFrameCount,m_iCubeFrameVisibleCount,m_iCubeFramePos);

	if (m_iCubeFrameCount < m_iTRange)
		m_iCubeFrameCount++;
	m_iCubeFrameVisibleCount = m_iCubeFrameCount;

	m_iCubeFramePos++;
	if (m_iCubeFramePos >= m_iTRange)
		m_iCubeFramePos -= m_iTRange;

	m_oaAtomFrames[m_iCubeFramePos] = frame;

	//m_IF.printf("PushAtomFrame after: count %d, visco %d, pos %d\n",m_iCubeFrameCount,m_iCubeFrameVisibleCount,m_iCubeFramePos);
}


/*
#include "df.h"


// Momentan beste Settings:
// travis compress cube -check no -cextra -cextrange 4 -cextorder 3 -cexsrange 5 -cexoffset 2 -cexsorder 1
//        -cexscross -cexdistexpo 2.5 -cexscrossrange -cexpredcorr -cexcorrfac 0.97 -cexwrap /home/brehm/cubefill/first10.cube new.bqb


void TestExtrapolator() {

	CBQBExtrapolator extra, extra2;
	CBQBCubeFrame cfr;
	int ix, iy, iz;
	double tf, tf2;
	FILE *a;
	CDF dfabs, dfrel;


	extra.Initialize(
		160,    // resx,
		160,    // resy,
		160,    // resz,
		3,      // srangex,
		3,      // srangey,
		3,      // srangez,
		5,      // trange,
		1,      // sorder,
		4,      // torder,
		1,      // offsetx,
		1,      // offsety,
		1,      // offsetz,
		false,   // crosss,
		false,  // crosst,
		false,   // wrap,
		false,   // crossranges,
		true,  // crossranget
		0.0,      // distexpo
		0.0,     // timeexpo
		1.0,  // corrfactor
		true,    // verbose
		false
		);

	return;

	extra.Initialize(
		160,    // resx,
		160,    // resy,
		160,    // resz,
		5,      // srangex,
		5,      // srangey,
		5,      // srangez,
		1,      // trange,
		2,      // sorder,
		0,      // torder,
		2,      // offsetx,
		2,      // offsety,
		2,      // offsetz,
		true,   // crosss,
		false,  // crosst,
		true,   // wrap,
		true,   // crossranges,
		false,  // crossranget
		6.0,      // distexpo
		0.0,    // timeexpo
		1.0, // corrfactor
		true,    // verbose
		false
		);

	a = fopen("first10.cube","rt");

	if (!cfr.ReadFrame(a,12,5,6,true)) {
		m_IF.printf("Error reading cube frame.\n");
		return;
	}

	fclose(a);

	extra.PushCubeFrame(&cfr);

	dfabs.m_iResolution = 100;
	dfabs.m_fMinVal = 0;
	dfabs.m_fMaxVal = 1.0;
	dfabs.Create();
	dfabs.SetLabelX("Absolute Deviation");
	dfabs.SetLabelY("Occurrence");

	dfrel.m_iResolution = 100;
	dfrel.m_fMinVal = 1.0;
	dfrel.m_fMaxVal = 2.0;
	dfrel.Create();
	dfrel.SetLabelX("Relative Deviation");
	dfrel.SetLabelY("Occurrence");

	m_IF.printf("Extrapolating:");

	std::vector<double> tda;

	tda.resize(160*160*160);

	for (ix=0;ix<160;ix++) {
		if ((ix%2) == 0)
			m_IF.printf(".");
		for (iy=0;iy<160;iy++) {
			for (iz=0;iz<160;iz++) {
		//		m_IF.printf("%d|%d|%d --> %d/%lu\n",ix,iy,iz,extra.FindPattern(0,ix,iy,iz),(unsigned long)extra.m_oaPatterns.size());
				tf = extra.Extrapolate(ix,iy,iz);
				tf2 = cfr.m_faBin[ix*160*160+iy*160+iz];
				dfabs.AddToBin(fabs(tf - tf2));
				if (tf != 0) {
					if (tf > tf2)
						dfrel.AddToBin(tf/tf2);
					else
						dfrel.AddToBin(tf2/tf);
				}
			}
		}
	}

	m_IF.printf("Done.\n");

	dfabs.NormBinIntegral(1000.0);
	dfabs.Write("","absf.csv","",false);

	dfrel.NormBinIntegral(1000.0);
	dfrel.Write("","relf.csv","",false);
}
*/



