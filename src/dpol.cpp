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

#include "dpol.h"
#include "timestep.h"
#include "maintools.h"
#include "conversion.h"


const char *GetRevisionInfo_dpol(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_dpol() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}





void ReRa_Neu2() {

	//std::vector<std::vector<double> > tempaa;
	std::vector<std::vector<std::vector<double> > > dip;
	std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > polre, polim;
	std::vector<int> smidx, tia;
	std::vector<std::vector<double> > acfre_iso, acfim_iso, acfre_aniso, acfim_aniso;
	std::vector<std::vector<double> > specre, specim;
	std::vector<double> tempin, tempin2, tempout, ev, /*afac,*/ tda;
	int z0, z, z2, m, n, acfdepth, o, rtpzps, bomdzps, bomdstart;
	int rtpwfn, mode, rmode, speclen, rty, procmode;
	int zm, za, zb, zt, zr, xfrom, xto, yfrom, yto;
	FILE *a;
	char buf[4096];
	bool crosscorr;
	CAutoCorrelation acf;
	CCrossCorrelation ccf;
	double tf=0, mwn, fac, tfs, bomdstride, rtpstride, estrength, contwn1, contwn2, contev1, contev2, rtpwfnparm;
	double freq, laserfreq, temperature;
	CFFT fft;
	CxString xbuf, path[3];
	CTimeStep *ts;
	CMolecule *pm;
	CSingleMolecule *sm;
	CBQBFile *bf;
	std::vector<int> afdir;
	const char *name = NULL;



	mprintf("\n");
	mprintf(WHITE,"    #########################################################\n");
	mprintf(WHITE,"    ####        Resonance Raman Prediction Engine        ####\n");
	mprintf(WHITE,"    ####          (c) Martin Brehm, 2019 - 2022          ####\n");
	mprintf(WHITE,"    ####            https://brehm-research.de            ####\n");
	mprintf(WHITE,"    #########################################################\n");
	mprintf("\n\n");

	for (z=0;z<g_oaMolecules.GetSize();z++) {
		((CMolecule*)g_oaMolecules[z])->m_iDipoleMode = 7;
		ParseAtom("#2", z, ((CMolecule *)g_oaMolecules[z])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[z])->m_iDipoleCenterIndex);
	}

	g_bVoroIntegrateCharge = true;
	g_bVoroIntegrateDipoleMoment = true;
	g_bVoroIntegrateQuadrupoleMoment = true;
	g_bVoroIntegrateTotalCurrent = false;
	g_bVoroIntegrateMagneticMoment = false;
	g_bCubeTimeDev = false;
	g_bUseVelocities = false;

	for (z=0;z<g_oaMolecules.GetSize();z++) {

		pm = (CMolecule*)g_oaMolecules[z];

		if (AskYesNo("    Observe molecules of type %s for this resonance Raman spectrum (y/n)? [yes] ",true,pm->m_sName)) {

			if (pm->m_laSingleMolIndex.GetSize() > 1) {
_msagain:
				AskString("      Which molecules of type %s to observe (e.g. 1,5-7, allowed range 1..%d)? [all] ",&xbuf,"*",pm->m_sName,pm->m_laSingleMolIndex.GetSize());
				if (mystricmp((const char*)xbuf,"*") == 0) {
					for (z2=0;z2<pm->m_laSingleMolIndex.GetSize();z2++)
						smidx.push_back(pm->m_laSingleMolIndex[z2]);
				} else {
					tia.clear();
					if (!ParseIntList((const char*)xbuf,tia))
						goto _msagain;
					for (z2=0;z2<(int)tia.size();z2++) {
						if ((tia[z2] < 1) || (tia[z2] > pm->m_laSingleMolIndex.GetSize())) {
							eprintf("Error: Molecule number %d out of allowed range.\n",tia[z2]);
							goto _msagain;
						}
						smidx.push_back(tia[z2]);
					}
				}
			} else
				smidx.push_back( pm->m_laSingleMolIndex[0] );
		}
	}

	mprintf("\n    Observing %lu molecules in total.\n\n",(unsigned long)smidx.size());

	if (smidx.size() > 1) {
		crosscorr = AskYesNo("    Include intermolecular cross-correlations (y/n)? [no] ",false);
		mprintf("\n");
	} else
		crosscorr = false;

	mprintf("    What you should have prepared as input:\n");
	mprintf("      (*) Performed RTP simulations for the step response under three different external field directions\n");
	mprintf("      (*) Wrote the total electron density along each RTP run to a volumetric trajectory\n");
	mprintf("      (*) Performed a Voronoi integration for each volumetric trajectory, *including* the SCF frame with field on.\n\n");

	mprintf("    The module expects all the \"properties.emp\" files from the Voronoi integrations of the RTP runs\n");
	mprintf("    in directories with the following path and name scheme:\n\n");

	mprintf("      /path/to/dpol_nnnnnn_d/properties.emp\n\n");

	mprintf("    /path/to is the prefix path to the subdirectories which can be chosen below.\n\n");
	mprintf("    nnnnnn is the BOMD snapshot number in six-digit format with leading zeroes, starting from 0 and consecutively.\n\n");
	mprintf("    d is the external field direction (x,y,z).\n\n");

	mprintf("\n");

	mode = AskRangeInteger("    Choose mode (0=XX, 1=YY, 2=ZZ, 3=XYZ): [3] ",0,3,3);

	mprintf("\n");

	if (mode == 3) {
		AskString("    Enter prefix path of RTP X sub-directories (absolute or relative): [.] ",&path[0],".");
		AskString("    Enter prefix path of RTP Y sub-directories (absolute or relative): [.] ",&path[1],".");
		AskString("    Enter prefix path of RTP Z sub-directories (absolute or relative): [.] ",&path[2],".");
	} else
		AskString("    Enter prefix path of RTP sub-directories (absolute or relative): [.] ",&path[0],".");

	if (mode == 3)
		rmode = AskRangeInteger("    Choose resonance Raman spectrum type (0=iso, 1=aniso, 2=unpol, 3=para, 4=ortho): [2] ",0,4,2);
	else
		rmode = 0;

	procmode = AskRangeInteger("    Choose processing mode (0=real, 1=complex, 2=abs.sqr., 3=abs.): [0] ",0,3,0);

	AskString("    Enter name of properties file in each directory: [properties.emp] ",&xbuf,"properties.emp");
	estrength = AskFloat("    Enter electric field strength that was applied (in a.u.): [5.0e-4] ",5.0e-4);
	temperature = AskFloat("    Enter temperature to compute spectrum for (in K): [350.0] ",350.0);

	mprintf("\n");

	bomdstride = AskFloat("    Enter BOMD snapshot time increment (in fs): [2.5] ",2.5);
	bomdstart = AskUnsignedInteger("    Start with which BOMD snapshot? [1] ",1) - 1;
	n = AskUnsignedInteger("    How many BOMD snapshots to evaluate from this position? [1600] ",1600);
	if (n < 256)
		acfdepth = n;
	else
		acfdepth = 256;
	acfdepth = AskUnsignedInteger("    Enter ACF depth along BOMD trajectory (in steps): [%d] ",acfdepth,acfdepth);
	bomdzps = AskUnsignedInteger("    Enter ACF length after zero padding for BOMD trajectory: [1024] ",1024);

	mprintf("\n");

	rtpstride = AskFloat("    Enter RTP trajectory time increment (in fs): [0.0625] ",0.0625);
	m = AskUnsignedInteger("    How many data points from the RTP time series to use? [256] ",256);
	rtpzps = AskUnsignedInteger("    ACF length after zero padding for RTP time series: [1024] ",1024);
	rtpwfn = AskUnsignedInteger("    Select RTP window function (0=cos^2, 1=exp, 2=triangular, 3=Gauss): [0] ",0);
	if (rtpwfn != 0)
		rtpwfnparm = AskFloat("    Enter parameter for window function: [4.0] ",4.0);
	else
		rtpwfnparm = 0;

	mprintf("\n");

	contwn1 = AskFloat("    Write contour plot starting from which wave number (in cm^-1): [200] ",200);
	contwn2 = AskFloat("    Write contour plot up to which wave number (in cm^-1): [1800] ",1800);
	contev1 = AskFloat("    Write contour plot starting from which laser energy (in eV): [1.0] ",1.0);
	contev2 = AskFloat("    Write contour plot up to which laser energy (in eV): [10.0] ",10.0);

	mprintf("\n");

	if (crosscorr) {

		eprintf("Cross-Correlations not yet implemented.\n");
		abort();

	} else {

		mprintf("    This will require %s of RAM. Allocating...\n",FormatBytes((double)smidx.size()*9*n*(rtpzps/2)*sizeof(double)));

		polre.resize( smidx.size() );
		polim.resize( smidx.size() );
		dip.resize( smidx.size() );

		for (zm=0;zm<(int)smidx.size();zm++) {

			polre[zm].resize(3);
			polim[zm].resize(3);
			dip[zm].resize(3);

			for (zb=0;zb<3;zb++)
				dip[zm][zb].resize(m);

			for (za=0;za<3;za++) {

				polre[zm][za].resize(3);
				polim[zm][za].resize(3);

				for (zb=0;zb<3;zb++) {

					polre[zm][za][zb].resize(n);
					polim[zm][za][zb].resize(n);

					for (zt=0;zt<n;zt++) {

						polre[zm][za][zb][zt].resize(rtpzps/2);
						polim[zm][za][zb][zt].resize(rtpzps/2);

						for (zr=0;zr<rtpzps/2;zr++) {

							polre[zm][za][zb][zt][zr] = 0;
							polim[zm][za][zb][zt][zr] = 0;
						}
					}
				}
			}
		}
	}

	mprintf("\n");

	fft.PrepareFFT_C2C(rtpzps);

	// From Debye to e*pm^2/V
	fac = DIP_DEBYE2CM / estrength / EFIELD_AU2VPM * 1.0E12 * 2.0 * Pi / CONST_ELEMENTARY_CHARGE;


	switch(mode) {
		case 0: // XX
			name = "xx";
			afdir.push_back(0);
			break;
		case 1: // YY
			name = "yy";
			afdir.push_back(1);
			break;
		case 2: // ZZ
			name = "zz";
			afdir.push_back(2);
			break;
		case 3: // Trace
			name = "trace";
			afdir.push_back(0);
			afdir.push_back(1);
			afdir.push_back(2);
			break;
		default:
			eprintf("Error: Unknown mode.\n");
			abort();
	}


	for (z0=0;z0<(int)afdir.size();z0++) {

		mprintf(WHITE,"    *** Now Field Direction %c ***\n\n",'X'+afdir[z0]);

		mprintf("      Reading data and performing FFT...\n");
		mprintf(WHITE,"        [");
		tfs = n/80.0;

		for (zt=0;zt<n;zt++) {

			if (fmod(zt,tfs) < 1.0)
				mprintf(WHITE,"#");

			sprintf(buf,"%s/dpol_%06d_%c/%s",(const char*)path[z0],zt+bomdstart+1,'x'+afdir[z0],(const char*)xbuf);

			bf = new CBQBFile(*g_pBQBInterface);
			ts = new CTimeStep();

			if (!bf->OpenRead(buf)) {
				eprintf("\nError: Could not open file \"%s\".\n",buf);
				abort();
			}

			zr = 0;
			while (ts->ReadBQB(bf,true)) {

				if (zr == 0) { // Coordinates do not change during RTP...
					ts->UniteMolecules(false);
					ts->CalcCenters();
				}

				ts->CalcDipoles(false);

				for (zm=0;zm<(int)smidx.size();zm++) {

					sm = (CSingleMolecule*)g_oaSingleMolecules[smidx[zm]];

					for (zb=0;zb<3;zb++)
						dip[zm][zb][zr] = sm->m_vDipole[zb];
				}

				zr++;

				if (zr >= m)
					break;
			}

			bf->Close();

			if (zr != m) {
				mprintf("\nError: Only %d instead of %d RTP samples in file %d.",zr,m,z+bomdstart);
				abort();
			}

			delete ts;
			delete bf;

			for (zm=0;zm<(int)smidx.size();zm++) {

				for (zb=0;zb<3;zb++) {

					for (zr=0;zr<m;zr++) {

						switch(rtpwfn) {

							case 0: // Cos^2
								tf = SQR(cos( (double)zr/m*Pi/2.0 ));
								break;

							case 1: // Exp
								tf = exp(-(double)zr/512.0*rtpwfnparm);
								break;

							case 2: // Triangular
								if (zr == 0)
									tf = 1.0;
								else
									tf = SQR(sin((double)zr/512.0*rtpwfnparm*Pi)/((double)zr/512.0*rtpwfnparm*Pi));
								break;

							case 3: // Gauss
								tf = exp(-SQR((double)zr/512.0*rtpwfnparm));
								break;
						}

						fft.m_pInput[zr*2] = (dip[zm][zb][zr] - dip[zm][zb][0]) * tf;
						fft.m_pInput[zr*2+1] = 0;
					}
					for (;zr<rtpzps;zr++) {
						fft.m_pInput[zr*2]   = 0;
						fft.m_pInput[zr*2+1] = 0;
					}

					fft.DoFFT();

					for (zr=0;zr<rtpzps/2;zr++) {
						polim[zm][afdir[z0]][zb][zt][zr] = fft.m_pOutput[zr*2]   * (double)zr/rtpzps * fac;
						polre[zm][afdir[z0]][zb][zt][zr] = fft.m_pOutput[zr*2+1] * (double)zr/rtpzps * fac;
					}

				} // FOR zb

			} // FOR zm

		} // FOR zt

		mprintf(WHITE,"] Done.\n\n");

		mprintf("      Read and transformed %d BOMD snapshots with %d RTP samples each.\n\n",n,m);

	} // FOR z0


	ev.resize(rtpzps/2);

	for (z=0;z<rtpzps/2;z++)
		ev[z] = (double)z/rtpzps * CONST_PLANCK / CONST_ELEMENTARY_CHARGE * 1.0E15 / rtpstride;

	m = rtpzps/2;


	mprintf("    Writing average polarizabilities...\n");

	sprintf(buf,"polre_%s.csv",name);
	a = fopen(buf,"wt");
	for (zr=0;zr<m;zr++) {
		fprintf(a,"%f",ev[zr]);
		for (zt=0;zt<n;zt++) {
			tf = 0;
			switch(mode) {
				case 0:
					for (zm=0;zm<(int)smidx.size();zm++)
						tf += polre[zm][0][0][zt][zr];
					break;
				case 1:
					for (zm=0;zm<(int)smidx.size();zm++)
						tf += polre[zm][1][1][zt][zr];
					break;
				case 2:
					for (zm=0;zm<(int)smidx.size();zm++)
						tf += polre[zm][2][2][zt][zr];
					break;
				case 3:
					for (zm=0;zm<(int)smidx.size();zm++) {
						tf += polre[zm][0][0][zt][zr] / 3.0;
						tf += polre[zm][1][1][zt][zr] / 3.0;
						tf += polre[zm][2][2][zt][zr] / 3.0;
					}
					break;
			}
			tf /= smidx.size();
			fprintf(a,";  %f",tf);
		}
		fprintf(a,"\n");
	}
	fclose(a);

	sprintf(buf,"polim_%s.csv",name);
	a = fopen(buf,"wt");
	for (zr=0;zr<m;zr++) {
		fprintf(a,"%f",ev[zr]);
		for (zt=0;zt<n;zt++) {
			tf = 0;
			switch(mode) {
				case 0:
					for (zm=0;zm<(int)smidx.size();zm++)
						tf += polim[zm][0][0][zt][zr];
					break;
				case 1:
					for (zm=0;zm<(int)smidx.size();zm++)
						tf += polim[zm][1][1][zt][zr];
					break;
				case 2:
					for (zm=0;zm<(int)smidx.size();zm++)
						tf += polim[zm][2][2][zt][zr];
					break;
				case 3:
					for (zm=0;zm<(int)smidx.size();zm++) {
						tf += polim[zm][0][0][zt][zr] / 3.0;
						tf += polim[zm][1][1][zt][zr] / 3.0;
						tf += polim[zm][2][2][zt][zr] / 3.0;
					}
					break;
			}
			tf /= smidx.size();
			fprintf(a,";  %f",tf);
		}
		fprintf(a,"\n");
	}
	fclose(a);

	mprintf("\n");

	mprintf("    Writing transposed average polarizabilities...\n");

	sprintf(buf,"transpolre_%s.csv",name);
	a = fopen(buf,"wt");
	for (zt=0;zt<n;zt++) {
		fprintf(a,"%f",(double)zt*bomdstride);
		for (zr=0;zr<m;zr++) {
			tf = 0;
			switch(mode) {
				case 0:
					for (zm=0;zm<(int)smidx.size();zm++)
						tf += polre[zm][0][0][zt][zr];
					break;
				case 1:
					for (zm=0;zm<(int)smidx.size();zm++)
						tf += polre[zm][1][1][zt][zr];
					break;
				case 2:
					for (zm=0;zm<(int)smidx.size();zm++)
						tf += polre[zm][2][2][zt][zr];
					break;
				case 3:
					for (zm=0;zm<(int)smidx.size();zm++) {
						tf += polre[zm][0][0][zt][zr] / 3.0;
						tf += polre[zm][1][1][zt][zr] / 3.0;
						tf += polre[zm][2][2][zt][zr] / 3.0;
					}
					break;
			}
			tf /= smidx.size();
			fprintf(a,";  %f",tf);
		}
		fprintf(a,"\n");
	}
	fclose(a);

	sprintf(buf,"transpolim_%s.csv",name);
	a = fopen(buf,"wt");
	for (zt=0;zt<n;zt++) {
		fprintf(a,"%f",(double)zt*bomdstride);
		for (zr=0;zr<m;zr++) {
			tf = 0;
			switch(mode) {
				case 0:
					for (zm=0;zm<(int)smidx.size();zm++)
						tf += polim[zm][0][0][zt][zr];
					break;
				case 1:
					for (zm=0;zm<(int)smidx.size();zm++)
						tf += polim[zm][1][1][zt][zr];
					break;
				case 2:
					for (zm=0;zm<(int)smidx.size();zm++)
						tf += polim[zm][2][2][zt][zr];
					break;
				case 3:
					for (zm=0;zm<(int)smidx.size();zm++) {
						tf += polim[zm][0][0][zt][zr] / 3.0;
						tf += polim[zm][1][1][zt][zr] / 3.0;
						tf += polim[zm][2][2][zt][zr] / 3.0;
					}
					break;
			}
			tf /= smidx.size();
			fprintf(a,";  %f",tf);
		}
		fprintf(a,"\n");
	}
	fclose(a);


	mprintf("\n    Computing Derivatives...\n");

	n -= 2;

	for (zm=0;zm<(int)smidx.size();zm++) {
		for (za=0;za<3;za++) {
			for (zb=0;zb<3;zb++) {
				for (zt=0;zt<n;zt++) {
					for (zr=0;zr<m;zr++) {
						polre[zm][za][zb][zt][zr] = (polre[zm][za][zb][zt+2][zr] - polre[zm][za][zb][zt][zr]) / 2.0;
						polim[zm][za][zb][zt][zr] = (polim[zm][za][zb][zt+2][zr] - polim[zm][za][zb][zt][zr]) / 2.0;
					}
				}
			}
		}
	}

	mprintf("\n");

	mprintf("    Autocorrelation...\n");

	mprintf(WHITE,"      [");

	acfre_iso.resize(m);
	acfre_aniso.resize(m);

	acf.Init(n,acfdepth,true);

	tempin.resize( n );
	tempout.resize( acfdepth );

	if (procmode == 1) { // Complex processing

		acfim_iso.resize(m);
		acfim_aniso.resize(m);
		ccf.Init(n,acfdepth,true);
		tempin2.resize( n );

	} else if (procmode == 2) { // Abs. Sqr.

		for (zm=0;zm<(int)smidx.size();zm++)
			for (za=0;za<3;za++)
				for (zb=0;zb<3;zb++)
					for (zt=0;zt<n;zt++)
						for (zr=0;zr<m;zr++)
							polre[zm][za][zb][zt][zr] = SQR(polre[zm][za][zb][zt][zr]) + SQR(polim[zm][za][zb][zt][zr]);

	} else if (procmode == 3) { // Abs.

		for (zm=0;zm<(int)smidx.size();zm++)
			for (za=0;za<3;za++)
				for (zb=0;zb<3;zb++)
					for (zt=0;zt<n;zt++)
						for (zr=0;zr<m;zr++)
							polre[zm][za][zb][zt][zr] = sqrt( SQR(polre[zm][za][zb][zt][zr]) + SQR(polim[zm][za][zb][zt][zr]) );
	}

	tfs = m / 80.0;

	for (zr=0;zr<m;zr++) {

		if (fmod(zr,tfs) < 1.0)
			mprintf(WHITE,"#");

		acfre_iso[zr].resize(acfdepth);
		acfre_aniso[zr].resize(acfdepth);

		for (zt=0;zt<acfdepth;zt++) {
			acfre_iso[zr][zt] = 0;
			acfre_aniso[zr][zt] = 0;
		}

		if (procmode == 1) {
			acfim_iso[zr].resize(acfdepth);
			acfim_aniso[zr].resize(acfdepth);
			for (zt=0;zt<acfdepth;zt++) {
				acfim_iso[zr][zt] = 0;
				acfim_aniso[zr][zt] = 0;
			}
		}

		for (zm=0;zm<(int)smidx.size();zm++) {

			// Isotropic Part
			for (zt=0;zt<n;zt++)
				tempin[zt] = (polre[zm][0][0][zt][zr] + polre[zm][1][1][zt][zr] + polre[zm][2][2][zt][zr]) / 3.0;
			acf.AutoCorrelate( tempin, tempout );
			for (zt=0;zt<acfdepth;zt++)
				acfre_iso[zr][zt] += tempout[zt];

			if (procmode == 1) {

				for (zt=0;zt<n;zt++)
					tempin2[zt] = (polim[zm][0][0][zt][zr] + polim[zm][1][1][zt][zr] + polim[zm][2][2][zt][zr]) / 3.0;
				acf.AutoCorrelate( tempin2, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfre_iso[zr][zt] += tempout[zt];
				ccf.CrossCorrelate( tempin, tempin2, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfim_iso[zr][zt] += tempout[zt];
				ccf.CrossCorrelate( tempin2, tempin, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfim_iso[zr][zt] -= tempout[zt];
			}


			// First Anisotropic Part
			for (zt=0;zt<n;zt++)
				tempin[zt] = polre[zm][0][0][zt][zr] - polre[zm][1][1][zt][zr];
			acf.AutoCorrelate( tempin, tempout );
			for (zt=0;zt<acfdepth;zt++)
				acfre_aniso[zr][zt] += tempout[zt] / 2.0;

			if (procmode == 1) {

				for (zt=0;zt<n;zt++)
					tempin2[zt] = polim[zm][0][0][zt][zr] - polim[zm][1][1][zt][zr];
				acf.AutoCorrelate( tempin2, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfre_aniso[zr][zt] += tempout[zt] / 2.0;
				ccf.CrossCorrelate( tempin, tempin2, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfim_aniso[zr][zt] += tempout[zt] / 2.0;
				ccf.CrossCorrelate( tempin2, tempin, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfim_aniso[zr][zt] -= tempout[zt] / 2.0;
			}

			// Second Anisotropic Part
			for (zt=0;zt<n;zt++)
				tempin[zt] = polre[zm][1][1][zt][zr] - polre[zm][2][2][zt][zr];
			acf.AutoCorrelate( tempin, tempout );
			for (zt=0;zt<acfdepth;zt++)
				acfre_aniso[zr][zt] += tempout[zt] / 2.0;

			if (procmode == 1) {

				for (zt=0;zt<n;zt++)
					tempin2[zt] = polim[zm][1][1][zt][zr] - polim[zm][2][2][zt][zr];
				acf.AutoCorrelate( tempin2, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfre_aniso[zr][zt] += tempout[zt] / 2.0;
				ccf.CrossCorrelate( tempin, tempin2, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfim_aniso[zr][zt] += tempout[zt] / 2.0;
				ccf.CrossCorrelate( tempin2, tempin, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfim_aniso[zr][zt] -= tempout[zt] / 2.0;
			}

			// Third Anisotropic Part
			for (zt=0;zt<n;zt++)
				tempin[zt] = polre[zm][2][2][zt][zr] - polre[zm][0][0][zt][zr];
			acf.AutoCorrelate( tempin, tempout );
			for (zt=0;zt<acfdepth;zt++)
				acfre_aniso[zr][zt] += tempout[zt] / 2.0;

			if (procmode == 1) {

				for (zt=0;zt<n;zt++)
					tempin2[zt] = polim[zm][2][2][zt][zr] - polim[zm][0][0][zt][zr];
				acf.AutoCorrelate( tempin2, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfre_aniso[zr][zt] += tempout[zt] / 2.0;
				ccf.CrossCorrelate( tempin, tempin2, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfim_aniso[zr][zt] += tempout[zt] / 2.0;
				ccf.CrossCorrelate( tempin2, tempin, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfim_aniso[zr][zt] -= tempout[zt] / 2.0;
			}

			// Fourth Anisotropic Part
			for (zt=0;zt<n;zt++)
				tempin[zt] = polre[zm][0][1][zt][zr];
			acf.AutoCorrelate( tempin, tempout );
			for (zt=0;zt<acfdepth;zt++)
				acfre_aniso[zr][zt] += tempout[zt] * 3.0;

			if (procmode == 1) {

				for (zt=0;zt<n;zt++)
					tempin2[zt] = polim[zm][0][1][zt][zr];
				acf.AutoCorrelate( tempin2, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfre_aniso[zr][zt] += tempout[zt] * 3.0;
				ccf.CrossCorrelate( tempin, tempin2, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfim_aniso[zr][zt] += tempout[zt] * 3.0;
				ccf.CrossCorrelate( tempin2, tempin, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfim_aniso[zr][zt] -= tempout[zt] * 3.0;
			}

			// Fifth Anisotropic Part
			for (zt=0;zt<n;zt++)
				tempin[zt] = polre[zm][1][2][zt][zr];
			acf.AutoCorrelate( tempin, tempout );
			for (zt=0;zt<acfdepth;zt++)
				acfre_aniso[zr][zt] += tempout[zt] * 3.0;

			if (procmode == 1) {

				for (zt=0;zt<n;zt++)
					tempin2[zt] = polim[zm][1][2][zt][zr];
				acf.AutoCorrelate( tempin2, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfre_aniso[zr][zt] += tempout[zt] * 3.0;
				ccf.CrossCorrelate( tempin, tempin2, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfim_aniso[zr][zt] += tempout[zt] * 3.0;
				ccf.CrossCorrelate( tempin2, tempin, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfim_aniso[zr][zt] -= tempout[zt] * 3.0;
			}

			// Sixth Anisotropic Part
			for (zt=0;zt<n;zt++)
				tempin[zt] = polre[zm][2][0][zt][zr];
			acf.AutoCorrelate( tempin, tempout );
			for (zt=0;zt<acfdepth;zt++)
				acfre_aniso[zr][zt] += tempout[zt] * 3.0;

			if (procmode == 1) {

				for (zt=0;zt<n;zt++)
					tempin2[zt] = polim[zm][2][0][zt][zr];
				acf.AutoCorrelate( tempin2, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfre_aniso[zr][zt] += tempout[zt] * 3.0;
				ccf.CrossCorrelate( tempin, tempin2, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfim_aniso[zr][zt] += tempout[zt] * 3.0;
				ccf.CrossCorrelate( tempin2, tempin, tempout );
				for (zt=0;zt<acfdepth;zt++)
					acfim_aniso[zr][zt] -= tempout[zt] * 3.0;
			}


		}

	}

	mprintf(WHITE,"] Done.\n\n");

	mprintf("    Writing ACF outputs...\n");

	sprintf(buf,"acfre_iso_%s.csv",name);
	a = fopen(buf,"wt");
	for (zt=0;zt<acfdepth;zt++) {
		fprintf(a,"%f",zt*bomdstride);
		for (zr=0;zr<m;zr++)
			fprintf(a,";  %f",acfre_iso[zr][zt]);
		fprintf(a,"\n");
	}
	fclose(a);

	sprintf(buf,"acfre_aniso_%s.csv",name);
	a = fopen(buf,"wt");
	for (zt=0;zt<acfdepth;zt++) {
		fprintf(a,"%f",zt*bomdstride);
		for (zr=0;zr<m;zr++)
			fprintf(a,";  %f",acfre_aniso[zr][zt]);
		fprintf(a,"\n");
	}
	fclose(a);

	if (procmode == 1) {

		sprintf(buf,"acfim_iso_%s.csv",name);
		a = fopen(buf,"wt");
		for (zt=0;zt<acfdepth;zt++) {
			fprintf(a,"%f",zt*bomdstride);
			for (zr=0;zr<m;zr++)
				fprintf(a,";  %f",acfim_iso[zr][zt]);
			fprintf(a,"\n");
		}
		fclose(a);

		sprintf(buf,"acfim_aniso_%s.csv",name);
		a = fopen(buf,"wt");
		for (zt=0;zt<acfdepth;zt++) {
			fprintf(a,"%f",zt*bomdstride);
			for (zr=0;zr<m;zr++)
				fprintf(a,";  %f",acfim_aniso[zr][zt]);
			fprintf(a,"\n");
		}
		fclose(a);
	}

	mprintf("    Computing spectra...\n");


	mwn = 5000.0;
	speclen = (int)(mwn / (33356.41 / bomdstride / 2.0) * bomdzps);

	mprintf("    Max. WN of %.2f cm^-1 --> speclen = %d / %d points\n",mwn,speclen,2*bomdzps);

	specre.resize(m);
	if (procmode == 1)
		specim.resize(m);
	fft.PrepareFFT_C2C(2*bomdzps);

	mprintf(WHITE,"      [");
	tfs = m / 80.0;

	for (zr=0;zr<m;zr++) {

		if (fmod(zr,tfs) < 1.0)
			mprintf(WHITE,"#");

		acfre_iso[zr].resize( 2*bomdzps );
		acfre_aniso[zr].resize( 2*bomdzps );

		for (zt=0;zt<acfdepth;zt++) {
			acfre_iso[zr][zt] *= SQR(cos( (double)zt/acfdepth * Pi/2.0 ));
			acfre_aniso[zr][zt] *= SQR(cos( (double)zt/acfdepth * Pi/2.0 ));
		}
		for (;zt<bomdzps;zt++) {
			acfre_iso[zr][zt] = 0;
			acfre_aniso[zr][zt] = 0;
		}
		o = bomdzps;
		for (zt=1;zt<o;zt++) {
			acfre_iso[zr][o+zt] = acfre_iso[zr][o-zt];
			acfre_aniso[zr][o+zt] = acfre_aniso[zr][o-zt];
		}
		acfre_iso[zr][o] = 0.0;
		acfre_aniso[zr][o] = 0.0;

		specre[zr].resize(speclen);

		if (procmode == 1) {

			acfim_iso[zr].resize( 2*bomdzps );
			acfim_aniso[zr].resize( 2*bomdzps );

			for (zt=0;zt<acfdepth;zt++) {
				acfim_iso[zr][zt] *= SQR(cos( (double)zt/acfdepth * Pi/2.0 ));
				acfim_aniso[zr][zt] *= SQR(cos( (double)zt/acfdepth * Pi/2.0 ));
			}
			for (;zt<bomdzps;zt++) {
				acfim_iso[zr][zt] = 0;
				acfim_aniso[zr][zt] = 0;
			}
			o = bomdzps;
			for (zt=1;zt<o;zt++) {
				acfim_iso[zr][o+zt] = acfim_iso[zr][o-zt];
				acfim_aniso[zr][o+zt] = acfim_aniso[zr][o-zt];
			}
			acfim_iso[zr][o] = 0.0;
			acfim_aniso[zr][o] = 0.0;

			specim[zr].resize(speclen);
		}


		switch(rmode) {

			case 0: // Iso
				for (zt=0;zt<2*bomdzps;zt++) {
					fft.m_pInput[zt*2] = acfre_iso[zr][zt];
					if (procmode == 1)
						fft.m_pInput[zt*2+1] = acfim_iso[zr][zt];
					else
						fft.m_pInput[zt*2+1] = 0;
				}
				break;

			case 1: // Aniso
				for (zt=0;zt<2*bomdzps;zt++) {
					fft.m_pInput[zt*2] = acfre_aniso[zr][zt];
					if (procmode == 1)
						fft.m_pInput[zt*2+1] = acfim_aniso[zr][zt];
					else
						fft.m_pInput[zt*2+1] = 0;
				}
				break;

			case 2: // Unpol
				for (zt=0;zt<2*bomdzps;zt++) {
					fft.m_pInput[zt*2] = acfre_iso[zr][zt] + 7.0/45.0 * acfre_aniso[zr][zt];
					if (procmode == 1)
						fft.m_pInput[zt*2+1] = acfim_iso[zr][zt] + 7.0/45.0 * acfim_aniso[zr][zt];
					else
						fft.m_pInput[zt*2+1] = 0;
				}
				break;

			case 3: // Para
				for (zt=0;zt<2*bomdzps;zt++) {
					fft.m_pInput[zt*2] = acfre_iso[zr][zt] + 4.0/45.0 * acfre_aniso[zr][zt];
					if (procmode == 1)
						fft.m_pInput[zt*2+1] = acfim_iso[zr][zt] + 4.0/45.0 * acfim_aniso[zr][zt];
					else
						fft.m_pInput[zt*2+1] = 0;
				}
				break;

			case 4: // Ortho
				for (zt=0;zt<2*bomdzps;zt++) {
					fft.m_pInput[zt*2] = acfre_aniso[zr][zt] / 15.0;
					if (procmode == 1)
						fft.m_pInput[zt*2+1] = acfim_aniso[zr][zt] / 15.0;
					else
						fft.m_pInput[zt*2+1] = 0;
				}
				break;
		}

		laserfreq = ev[zr] * 8065.5; // in cm^-1

		fft.DoFFT();

		for (zt=0;zt<speclen;zt++) {

			freq = (double)zt/speclen*mwn; // in cm^-1

			if ((zr == 0) || (zt == 0))
				tf = 0;
			else
				tf = 1.0E40
					* 100.0
					* 2.0
					* CONST_PLANCK
					/ (8.0 * CONST_BOLTZMANN * CONST_EPSILON0 * CONST_EPSILON0)
					* 1.0E6
					* CONST_ELEMENTARY_CHARGE * CONST_ELEMENTARY_CHARGE
					* 1.0E-48
					/ 1.0E-15
					* pow4(laserfreq - freq)
					/ freq
					/ (1.0 - exp(-1.438777 * freq / temperature))
					* bomdstride; // Output in 1e-30*m^2*cm

			specre[zr][zt] = tf * fft.m_pOutput[zt*2];

			if (procmode == 1)
				specim[zr][zt] = tf * fft.m_pOutput[zt*2+1];
		}

		// Finite differences correction
		fac = mwn/speclen * bomdstride * 1.883652e-4;

		for (zt=1;zt<speclen;zt++)
			specre[zr][zt] *= pow2(fac * zt / sin(fac * zt)); // Divide by sinc function to correct finite difference derivation

		if (procmode == 1)
			for (zt=1;zt<speclen;zt++)
				specim[zr][zt] *= pow2(fac * zt / sin(fac * zt));
	}

	mprintf(WHITE,"] Done.\n\n");

	mprintf("    Writing spectra...\n");

	sprintf(buf,"specre_%s.csv",name);
	a = fopen(buf,"wt");
	for (zt=0;zt<speclen;zt++) {
		fprintf(a,"%f",(double)zt/speclen*mwn);
		for (zr=0;zr<m;zr++)
			fprintf(a,";  %f",specre[zr][zt]);
		fprintf(a,"\n");
	}
	fclose(a);

	if (procmode == 1) {

		sprintf(buf,"specim_%s.csv",name);
		a = fopen(buf,"wt");
		for (zt=0;zt<speclen;zt++) {
			fprintf(a,"%f",(double)zt/speclen*mwn);
			for (zr=0;zr<m;zr++)
				fprintf(a,";  %f",specim[zr][zt]);
			fprintf(a,"\n");
		}
		fclose(a);

		sprintf(buf,"specabs_%s.csv",name);
		a = fopen(buf,"wt");
		for (zt=0;zt<speclen;zt++) {
			fprintf(a,"%f",(double)zt/speclen*mwn);
			for (zr=0;zr<m;zr++)
				fprintf(a,";  %f",sqrt(SQR(specre[zr][zt])+SQR(specim[zr][zt])));
			fprintf(a,"\n");
		}
		fclose(a);
	}

	mprintf("\n");


	C2DF df;

	xfrom = (int)(contwn1/mwn*speclen);
	xto = (int)(contwn2/mwn*speclen);

	if (xto > speclen)
		xto = speclen;

	df.m_iRes[0] = xto-xfrom+1;

	mprintf("    X axis range from %7.2f to %7.2f cm^-1 is data point %4d to %4d.\n",contwn1,contwn2,xfrom,xto);

	yfrom = -1;
	yto = -1;
	for (z=0;z<m;z++) {
		if (yfrom == -1)
			if (ev[z] >= contev1)
				yfrom = z;
		if (yto == -1) {
			if (ev[z] >= contev2) {
				yto = z;
				break;
			}
		}
	}
	df.m_iRes[1] =  yto - yfrom + 1;

	mprintf("    Y axis range from %7.2f to %7.2f eV    is data point %4d to %4d.\n\n",contev1,contev2,yfrom,yto);

	mprintf("    Using a resolution of %d x %d for the contour plots.\n\n",df.m_iRes[0],df.m_iRes[1]);

	mprintf("    Writing ReRa contour plot...\n");

	df.m_fMinVal[0] = contwn1;
	//df.m_fMaxVal[0] = ((double)df.m_iRes[0]-1.0)/speclen*mwn;
	df.m_fMaxVal[0] = contwn2;
	df.m_fMinVal[1] = contev1;
	//df.m_fMaxVal[1] = ev[df.m_iRes[1]];
	df.m_fMaxVal[1] = contev2;
	df.Create();
	df.SetLabelX("Wavenumber / cm^-1");
	df.SetLabelY("Laser Energy / eV");
	df.m_iPlotType = 2;
	df.m_iSmoothGrade = 0;
	df.m_iInterpolationOrder = 0;
	df.m_bDrawMesh = false;

	for (z=0;z<df.m_iRes[0];z++)
		for (z2=0;z2<df.m_iRes[1];z2++)
			df.m_pBin[z2*df.m_iRes[0]+z] = specre[z2+yfrom][z+xfrom];

	df.m_fPlotExp = 0.5;

	sprintf(buf,"specre_contour_%s.nb",name);
	df.WriteMathematicaNb("",buf,"",false);

	sprintf(buf,"specre_contour_%s",name);
	df.WriteGnuplotInput("",buf,"",false);



	for (z=0;z<df.m_iRes[0];z++) {
		for (z2=0;z2<df.m_iRes[1];z2++) {
			if (specre[z2+yfrom][z+xfrom] < 1000.0)
				df.m_pBin[z2*df.m_iRes[0]+z] = 0.0;
			else
				df.m_pBin[z2*df.m_iRes[0]+z] = log10(specre[z2+yfrom][z+xfrom])-3.0;
		}
	}

	df.m_fPlotExp = 1.0;

	sprintf(buf,"specre_contour_log_%s.nb",name);
	df.WriteMathematicaNb("",buf,"",false);

	sprintf(buf,"specre_contour_log_%s",name);
	df.WriteGnuplotInput("",buf,"",false);



	for (z2=0;z2<df.m_iRes[1];z2++) {
		tf = 0;
		for (z=0;z<df.m_iRes[0];z++)
			if (specre[z2+yfrom][z+xfrom] > 0)
				tf += specre[z2+yfrom][z+xfrom];
		for (z=0;z<df.m_iRes[0];z++)
			if (specre[z2+yfrom][z+xfrom] > 0)
				df.m_pBin[z2*df.m_iRes[0]+z] = specre[z2+yfrom][z+xfrom] * 1000.0/tf;
			else
				df.m_pBin[z2*df.m_iRes[0]+z] = 0;
	}

	df.m_fPlotExp = 0.7;

	sprintf(buf,"specre_contour_rownorm_%s.nb",name);
	df.WriteMathematicaNb("",buf,"",false);

	sprintf(buf,"specre_contour_rownorm_%s",name);
	df.WriteGnuplotInput("",buf,"",false);



	for (z2=0;z2<df.m_iRes[1];z2++) {
		tf = 0;
		for (z=0;z<df.m_iRes[0];z++)
			if (tf < specre[z2+yfrom][z+xfrom])
				tf = specre[z2+yfrom][z+xfrom];
		for (z=0;z<df.m_iRes[0];z++)
			if (specre[z2+yfrom][z+xfrom] > 0)
				df.m_pBin[z2*df.m_iRes[0]+z] = specre[z2+yfrom][z+xfrom] * 1000.0/tf;
			else
				df.m_pBin[z2*df.m_iRes[0]+z] = 0;
	}

	df.m_fPlotExp = 0.7;

	sprintf(buf,"specre_contour_rowmax_%s.nb",name);
	df.WriteMathematicaNb("",buf,"",false);

	sprintf(buf,"specre_contour_rowmax_%s",name);
	df.WriteGnuplotInput("",buf,"",false);


	sprintf(buf,"specre_%s_normpeak.csv",name);
	a = fopen(buf,"wt");
	tda.resize(m);
	for (zr=0;zr<m;zr++) {
		tf = 0;
		for (zt=xfrom;zt<=xto;zt++)
			if (specre[zr][zt] > tf)
				tf = specre[zr][zt];
		if (tf < 1.0e-6)
			tda[zr] = 1.0e-6;
		else
			tda[zr] = tf;
	}
	for (zt=0;zt<speclen;zt++) {
		fprintf(a,"%f",(double)zt/speclen*mwn);
		for (zr=0;zr<m;zr++)
			fprintf(a,";  %f",specre[zr][zt]/tda[zr]*1000.0);
		fprintf(a,"\n");
	}
	fclose(a);


	mprintf("\n");
	mprintf(WHITE,"    #########################################################\n");
	mprintf(WHITE,"    ####    Resonance Raman Prediction Engine leaving    ####\n");
	mprintf(WHITE,"    #########################################################\n");
	mprintf("\n");
}







