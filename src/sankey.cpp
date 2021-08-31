/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

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

#include "sankey.h"
#include "tools.h"
#include "svgwriter.h"
#include "xstring.h"


const char *GetRevisionInfo_sankey(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_sankey() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



CSankeyDiagramEngine::CSankeyDiagramEngine() {

	m_bSmooth = true;
	m_bSmoothHSV = true;
	m_fPadding = 100.0;
}



void CSankeyDiagramEngine::BuildSankeyDiagram( const char *s ) {

	int z, z2, i1, i2, ti;
	int count;
	CxString buf;
	const char *defaultcolors[7];
	char *p, *q;
	unsigned long rs;
	std::vector<std::string> sadesc;
	std::vector<std::vector<double> > favalues;
	std::vector<int> iasource, iasink;
	std::vector<double> fasum;
	char cbuf[256];
	FILE *a;


	defaultcolors[0] = "FF6020";
	defaultcolors[1] = "00FF00";
	defaultcolors[2] = "3030FF";
	defaultcolors[3] = "FFFF00";
	defaultcolors[4] = "C000E0";
	defaultcolors[5] = "A0A0A0";
	defaultcolors[6] = "00FFFF";


	mprintf("\n");
	mprintf(WHITE,"    #############################################\n");
	mprintf(WHITE,"    ####        Sankey Diagram Engine        ####\n");
	mprintf(WHITE,"    ####    (c) Martin Brehm, 2020 - 2021    ####\n");
	mprintf(WHITE,"    ####      https://brehm-research.de      ####\n");
	mprintf(WHITE,"    #############################################\n");
	mprintf("\n");

	if (AskYesNo("    Create a linear (left-to-right) Sankey diagram (y) or a circular Sankey diagram (n)? [yes] ",true)) {

		m_iType = 0;

		mprintf("\n");

		if (s != NULL) {
			buf = s;
			goto _hafil;
		}

		if (AskYesNo("    Read data from a text file (y) or enter it manually / use random data (n)? [no] ",false)) {

_f1again:
			mprintf("\n");

			AskString_ND("    Enter the name of the data file: ",&buf);

_hafil:
			mprintf("\n    Reading file \"%s\" ...\n",(const char*)buf);

			a = fopen((const char*)buf,"rt");
			if (a == NULL) {
				eprintf("Error: Could not open \"%s\" for reading.\n",(const char*)buf);
				if (s != NULL)
					abort();
				goto _f1again;
			}
			while (!feof(a)) {
				(void)!fgets(cbuf,256,a);
				if (feof(a)) {
					eprintf("Error: Unexpected end of file (A).\n");
					abort();
				}
				if ((cbuf[0] != '#') && (cbuf[0] != 0) && (cbuf[0] != '\r') && (cbuf[0] != '\n'))
					break;
			}

			if (!IsValidInteger( cbuf )) {
				eprintf("Error: Invalid integer number: \"%s\".\n",cbuf);
				abort();
			}

			count = atoi(cbuf);

			favalues.resize( count );
			fasum.resize( count );

			for (z=0;z<count;z++) {

				while (!feof(a)) {
					(void)!fgets(cbuf,256,a);
					if (feof(a)) {
						eprintf("Error: Unexpected end of file (B, z=%d).\n",z);
						abort();
					}
					if ((cbuf[0] != '#') && (cbuf[0] != 0) && (cbuf[0] != '\r') && (cbuf[0] != '\n'))
						break;
				}
				p = cbuf;
				while (*p == ' ')
					p++;
				q = p;
				while (*q != 0)
					q++;
				q--;
				while ((*q == ' ') || (*q == '\r') || (*q == '\n'))
					q--;
				q++;
				*q = 0;
				sadesc.push_back( p );

				favalues[z].resize( count );
				fasum[z] = 0;

				for (z2=0;z2<count;z2++) {

					while (!feof(a)) {
						(void)!fgets(cbuf,256,a);
						if (feof(a)) {
							eprintf("Error: Unexpected end of file (C, z=%d, z2=%d).\n",z,z2);
							abort();
						}
						if ((cbuf[0] != '#') && (cbuf[0] != 0) && (cbuf[0] != '\r') && (cbuf[0] != '\n'))
							break;
					}

					if (!IsValidFloat( cbuf )) {
						eprintf("Error: Invalid floating point number: \"%s\".\n",cbuf);
						abort();
					}

					favalues[z][z2] = atof( cbuf );
					fasum[z] += favalues[z][z2];
				}
			}
			fclose(a);

			mprintf("\n");
			mprintf("    Read %d groups from file:\n\n",count);

			for (z=0;z<count;z++) {
				mprintf(WHITE,"      (*) Group %2d:",z+1);
				mprintf("  Outgoing sum %7.4f,  %s\n",fasum[z],sadesc[z].c_str());
			}

			mprintf("\n");

_ila1:
			AskString_ND( "    Which of the groups to use for the left (source) side (e.g. \"1,3-5\")? ",&buf);
			if (!ParseIntList( (const char*)buf, iasource )) {
				eprintf("    Error: Could not parse integer list.\n");
				goto _ila1;
			}
			for (z=0;z<(int)iasource.size();z++) {
				if ((iasource[z] < 1) || (iasource[z] > count)) {
					eprintf("Error: Group index %d out of allowed range (1 .. %d).\n",iasource[z],count);
					goto _ila1;
				}
			}

_ila2:
			AskString_ND( "    Which of the groups to use for the right (sink) side (e.g. \"1,3-5\")? ",&buf);
			if (!ParseIntList( (const char*)buf, iasink )) {
				eprintf("    Error: Could not parse integer list.\n");
				goto _ila2;
			}
			for (z=0;z<(int)iasink.size();z++) {
				if ((iasink[z] < 1) || (iasink[z] > count)) {
					eprintf("Error: Group index %d out of allowed range (1 .. %d).\n",iasink[z],count);
					goto _ila2;
				}
			}

			for (z=0;z<(int)iasource.size();z++)
				iasource[z]--;

			for (z=0;z<(int)iasink.size();z++)
				iasink[z]--;

			m_iSourceCount = (int)iasource.size();
			m_iSinkCount = (int)iasink.size();

			m_faMatrixSource.resize( m_iSourceCount*m_iSinkCount );
			m_faMatrixSink.resize( m_iSourceCount*m_iSinkCount );
			m_iaGroupPermutationSource.resize(m_iSourceCount);
			m_iaGroupPermutationSink.resize(m_iSinkCount);


			for (z=0;z<m_iSinkCount;z++) {
				for (z2=0;z2<m_iSourceCount;z2++) {
					m_faMatrixSource[ z2*m_iSinkCount + z ] = favalues[ iasource[z2] ][ iasink[z] ];
					m_faMatrixSink[ z2*m_iSinkCount + z ] = favalues[ iasink[z] ][ iasource[z2] ];
				}
			}

			for (z=0;z<m_iSourceCount;z++)
				m_iaGroupPermutationSource[z] = z;

			for (z=0;z<m_iSinkCount;z++)
				m_iaGroupPermutationSink[z] = z;

		} else {

			mprintf("\n");

			m_iSourceCount = AskUnsignedInteger("    Enter source (left) group count of Sankey diagram: [5] ",5);
			m_iSinkCount = AskUnsignedInteger("    Enter sink (right) group count of Sankey diagram: [5] ",5);

			m_faMatrixSource.resize( m_iSourceCount*m_iSinkCount );
			m_faMatrixSink.resize( m_iSourceCount*m_iSinkCount );
			m_iaGroupPermutationSource.resize(m_iSourceCount);
			m_iaGroupPermutationSink.resize(m_iSinkCount);

			mprintf("\n");

			if (AskYesNo("    Create random diagram (y/n)? [no] ",false)) {

				mprintf("\n");

				m_fZeroThreshold = 0.1;

				ti = AskInteger("    Enter random seed for the \"random\" data (-1 = \"truly random\"): [-1] ",-1);

				mprintf("\n");

				if (ti < 0) {
					rs = (unsigned long)time(NULL);
					mprintf("\n\n    Using %lu as random seed.\n",rs);
					srand(rs);
					mprintf(WHITE,"    Hint: ");
					mprintf("Try it multiple times! The diagram will look different each time.\n\n");
				} else
					srand(ti);

				for (z=0;z<m_iSourceCount*m_iSinkCount;z++) {
					m_faMatrixSource[z] = (700 + rand()%5000) / 5000.0;
					m_faMatrixSink[z] = (700 + rand()%5000) / 5000.0;
				}

				if (!AskYesNo("    Populate all connections (y) or only some (n)? [no] ",false)) {
					for (z=0;z<(m_iSourceCount*m_iSinkCount/2);z++) {
						i1 = rand()%m_iSourceCount;
						i2 = rand()%m_iSinkCount;
						m_faMatrixSource[ i1*m_iSinkCount + i2 ] = 0;
						m_faMatrixSink[ i1*m_iSinkCount + i2 ] = 0;
					}
				}

				for (z=0;z<m_iSourceCount;z++)
					m_iaGroupPermutationSource[z] = z;

				for (z=0;z<m_iSinkCount;z++)
					m_iaGroupPermutationSink[z] = z;
			}
		}

		mprintf("\n");
		m_fZeroThreshold = AskFloat("    Enter zero threshold for entries: [0.01] ",0.01);

		mprintf("\n");
		m_iaColorsSource.resize(m_iSourceCount);
		for (z=0;z<m_iSourceCount;z++) {
_colagainsource:
			AskString("    Enter color for source group %2d (in RGB-Hex, 6 digits): [%s] ",&buf,defaultcolors[z%7],z+1,defaultcolors[z%7]);
			if (buf.GetLength() != 6) {
				eprintf("Error: Please enter exactly 6 hexadecimal digits.\n");
				goto _colagainsource;
			}
			m_iaColorsSource[z] = strtol( buf.GetWritePointer(), &p, 16 );
			if (*p != 0) {
				eprintf("Error: Invalid character \"%c\" found. Enter only characters 0 .. F.\n",*p);
				goto _colagainsource;
			}
		}

		mprintf("\n");

		m_iaColorsSink.resize(m_iSinkCount);
		for (z=0;z<m_iSinkCount;z++) {
_colagainsink:
			AskString("    Enter color for sink group %2d (in RGB-Hex, 6 digits): [%s] ",&buf,defaultcolors[z%7],z+1,defaultcolors[z%7]);
			if (buf.GetLength() != 6) {
				eprintf("Error: Please enter exactly 6 hexadecimal digits.\n");
				goto _colagainsink;
			}
			m_iaColorsSink[z] = strtol( buf.GetWritePointer(), &p, 16 );
			if (*p != 0) {
				eprintf("Error: Invalid character \"%c\" found. Enter only characters 0 .. F.\n",*p);
				goto _colagainsink;
			}
		}

		mprintf("\n");
		m_bCorrectWidth = AskYesNo("    Perform width correction of the connections (y/n)? [yes] ",true);

		mprintf("\n");
		m_bCorrectRatio = AskYesNo("    Keep correct ratio between sink and source (y/n)? [yes] ",true);


/*****************************************************************************************************************************************/
	} else { // Circular Sankey Diagram
/*****************************************************************************************************************************************/


		m_iType = 1;

		mprintf("\n");

		if (s != NULL) {
			buf = s;
			goto _hafic;
		}

		if (AskYesNo("    Read data from a text file (y) or enter it manually / use random data (n)? [no] ",false)) {

_f2again:
			mprintf("\n");

			AskString_ND("    Enter the name of the data file: ",&buf);

_hafic:
			mprintf("\n    Reading file \"%s\" ...\n",(const char*)buf);

			a = fopen((const char*)buf,"rt");
			if (a == NULL) {
				eprintf("Error: Could not open \"%s\" for reading.\n",(const char*)buf);
				if (s != NULL)
					abort();
				goto _f2again;
			}
			while (!feof(a)) {
				(void)!fgets(cbuf,256,a);
				if (feof(a)) {
					eprintf("Error: Unexpected end of file (A).\n");
					abort();
				}
				if ((cbuf[0] != '#') && (cbuf[0] != 0) && (cbuf[0] != '\r') && (cbuf[0] != '\n'))
					break;
			}

			if (!IsValidInteger( cbuf )) {
				eprintf("Error: Invalid integer number: \"%s\".\n",cbuf);
				abort();
			}

			count = atoi(cbuf);

			favalues.resize( count );
			fasum.resize( count );

			for (z=0;z<count;z++) {

				while (!feof(a)) {
					(void)!fgets(cbuf,256,a);
					if (feof(a)) {
						eprintf("Error: Unexpected end of file (B, z=%d).\n",z);
						abort();
					}
					if ((cbuf[0] != '#') && (cbuf[0] != 0) && (cbuf[0] != '\r') && (cbuf[0] != '\n'))
						break;
				}
				p = cbuf;
				while (*p == ' ')
					p++;
				q = p;
				while (*q != 0)
					q++;
				q--;
				while ((*q == ' ') || (*q == '\r') || (*q == '\n'))
					q--;
				q++;
				*q = 0;
				sadesc.push_back( p );

				favalues[z].resize( count );
				fasum[z] = 0;

				for (z2=0;z2<count;z2++) {

					while (!feof(a)) {
						(void)!fgets(cbuf,256,a);
						if (feof(a)) {
							eprintf("Error: Unexpected end of file (C, z=%d, z2=%d).\n",z,z2);
							abort();
						}
						if ((cbuf[0] != '#') && (cbuf[0] != 0) && (cbuf[0] != '\r') && (cbuf[0] != '\n'))
							break;
					}

					if (!IsValidFloat( cbuf )) {
						eprintf("Error: Invalid floating point number: \"%s\".\n",cbuf);
						abort();
					}

					favalues[z][z2] = atof( cbuf );
					fasum[z] += favalues[z][z2];
				}
			}
			fclose(a);

			mprintf("\n");
			mprintf("    Read %d groups from file:\n\n",count);

			for (z=0;z<count;z++) {
				mprintf(WHITE,"      (*) Group %2d:",z+1);
				mprintf("  Outgoing sum %7.4f,  %s\n",fasum[z],sadesc[z].c_str());
			}

			m_iNodeCount = count;

			m_faMatrix.resize( m_iNodeCount*m_iNodeCount );
			m_iaGroupPermutation.resize(m_iNodeCount);

			for (z=0;z<m_iNodeCount*m_iNodeCount;z++)
				m_faMatrix[z] = 0;

			for (z=0;z<m_iNodeCount;z++)
				for (z2=0;z2<m_iNodeCount;z2++)
					m_faMatrix[ z2*m_iNodeCount + z ] = favalues[z][z2];

		} else {

			mprintf("\n");

			m_iNodeCount = AskUnsignedInteger("    Enter group count of Sankey diagram: [5] ",5);

			m_faMatrix.resize( m_iNodeCount*m_iNodeCount );
			m_iaGroupPermutation.resize(m_iNodeCount);

			mprintf("\n");

			if (AskYesNo("    Create random diagram (y/n)? [no] ",false)) {

				mprintf("\n");

				m_fZeroThreshold = 0.1;

				ti = AskInteger("    Enter random seed for the \"random\" data (-1 = \"truly random\"): [-1] ",-1);

				mprintf("\n");

				if (ti < 0) {
					rs = (unsigned long)time(NULL);
					mprintf("\n\n    Using %lu as random seed.\n",rs);
					srand(rs);
					mprintf(WHITE,"    Hint: ");
					mprintf("Try it multiple times! The diagram will look different each time.\n\n");
				} else
					srand(ti);

				for (z=0;z<m_iNodeCount*m_iNodeCount;z++)
					m_faMatrix[z] = (700 + rand()%5000) / 5000.0;

				for (z=0;z<m_iNodeCount;z++)
					m_faMatrix[ z*m_iNodeCount + z ] = 0;

				if (!AskYesNo("    Populate all connections (y) or only some (n)? [no] ",false)) {
					for (z=0;z<(m_iNodeCount*(m_iNodeCount-1)/4);z++) {
						do {
							i1 = rand()%m_iNodeCount;
							i2 = rand()%m_iNodeCount;
						} while (i1 == i2);
						m_faMatrix[ i1*m_iNodeCount + i2 ] = 0;
						m_faMatrix[ i2*m_iNodeCount + i1 ] = 0;
					}
				}

				for (z=0;z<m_iNodeCount;z++)
					m_iaGroupPermutation[z] = z;

			} else {

				mprintf("\n");

				for (z=0;z<m_iNodeCount*m_iNodeCount;z++)
					m_faMatrix[z] = 0;

				for (z=0;z<m_iNodeCount;z++) {
					mprintf(WHITE,"    Group %2d:\n",z+1);
					for (z2=0;z2<m_iNodeCount;z2++) {
						if (z == z2)
							continue;
						m_faMatrix[ z2*m_iNodeCount + z ] = AskFloat("      Enter connection weight to group %2d: [0.0] ",0,z2+1);
					}
					mprintf("\n");
				}
			}
		}

		mprintf("\n");
		m_fZeroThreshold = AskFloat("    Enter zero threshold for entries: [0.01] ",0.01);

		mprintf("\n");
		m_iaColors.resize(m_iNodeCount);
		for (z=0;z<m_iNodeCount;z++) {
_colagain:
			AskString("    Enter color for group %2d (in RGB-Hex, 6 digits): [%s] ",&buf,defaultcolors[z%7],z+1,defaultcolors[z%7]);
			if (buf.GetLength() != 6) {
				eprintf("Error: Please enter exactly 6 hexadecimal digits.\n");
				goto _colagain;
			}
			m_iaColors[z] = strtol( buf.GetWritePointer(), &p, 16 );
			if (*p != 0) {
				eprintf("Error: Invalid character \"%c\" found. Enter only characters 0 .. F.\n",*p);
				goto _colagain;
			}
		}

		mprintf("\n");
		if (AskYesNo("    Permute the groups (y/n)? [no] ",false)) {
			mprintf("\n");
			for (z=0;z<m_iNodeCount;z++)
				m_iaGroupPermutation[z] = AskRangeInteger("      Which group to use in position %d (1-%d)? [%d] ",1,m_iNodeCount,z+1,z+1,m_iNodeCount,z+1) - 1;
		} else
			for (z=0;z<m_iNodeCount;z++)
				m_iaGroupPermutation[z] = z;
	}

	if (m_iType == 0) {

		mprintf("\n");
		mprintf(WHITE,"    *** Source Matrix ***\n\n");
		for (z=0;z<m_iSinkCount;z++) {
			mprintf("      ");
			for (z2=0;z2<m_iSourceCount;z2++)
				mprintf("  %10.6f", m_faMatrixSource[ z2*m_iSinkCount + z ] );
			mprintf("\n");
		}

		mprintf("\n\n");
		mprintf(WHITE,"    *** Sink Matrix ***\n\n");
		for (z=0;z<m_iSinkCount;z++) {
			mprintf("      ");
			for (z2=0;z2<m_iSourceCount;z2++)
				mprintf("  %10.6f", m_faMatrixSink[ z2*m_iSinkCount + z ] );
			mprintf("\n");
		}

	} else {

		mprintf("\n");
		mprintf(WHITE,"    *** Data Matrix ***\n\n");
		for (z=0;z<m_iNodeCount;z++) {
			mprintf("      ");
			for (z2=0;z2<m_iNodeCount;z2++)
				mprintf("  %10.6f", m_faMatrix[ z2*m_iNodeCount + z ] );
			mprintf("\n");
		}
	}

	mprintf("\n\n");

	if (m_iType == 0) // Linear
		m_fPadding = AskFloat("    Enter padding percentage: [50.0] ",50.0);
	else // Circular
		m_fPadding = AskFloat("    Enter padding percentage: [100.0] ",100.0);

	if (AskYesNo("    Use color gradients for connections (y/n)? [yes] ",true)) {
		m_bSmooth = true;
		m_bSmoothHSV = AskYesNo("    Interpolate colors in HSV (y) or RGB (n) space? [yes] ",true);
	} else
		m_bSmooth = false;

	mprintf("\n");

	if (m_iType == 1) {
		mprintf("    Writing Sankey diagram to \"sankey_circular.svg\"...\n");
		WriteSVG( "sankey_circular.svg" );
	} else {
		mprintf("    Writing Sankey diagram to \"sankey_linear.svg\"...\n");
		WriteSVG( "sankey_linear.svg" );
	}

	mprintf("    Use a web browser or a vector drawing program to open the diagram!\n\n");

	mprintf("\n");
	mprintf(WHITE,"    #############################################\n");
	mprintf(WHITE,"    ####    Sankey Diagram Engine leaving    ####\n");
	mprintf(WHITE,"    #############################################\n");
	mprintf("\n");
}



void ArcPart( CSVGObject_Path *p, double mainrad, double a1, double a2 ) {

	double r2, d;
	int f;


	r2 = mainrad * tan( fabs((a2-a1)/2.0) );

	d = a2-a1;
	if (d > Pi)
		d -= 2.0*Pi;
	if (d < -Pi)
		d += 2.0*Pi;
	if (d > 0)
		f = 0;
	else
		f = 1;
	p->AddArcTo( fabs(r2), fabs(r2), 0, 0, f, 500.0+mainrad*cos(a2), 500.0+mainrad*sin(a2) );
}



unsigned long ColorMix(unsigned long col1, unsigned long col2) {

	unsigned char r, g, b;


	r = (unsigned char)((col1>>17)+(col2>>17));
	g = (unsigned char)(((col1&0xFFFF)>>9)+((col2&0xFFFF)>>9));
	b = (unsigned char)(((col1&0xFF)>>1)+((col2&0xFF)>>1));
	return (r<<16)+(g<<8)+b;
}



unsigned long ColorInterpolate(unsigned long col1, unsigned long col2, double f, bool hsv) {

	unsigned char r1, g1, b1, r2, g2, b2, r, g, b;
	double th, fr, fg, fb, h1, s1, v1, h2, s2, v2;


	r1 = (unsigned char)(col1>>16);
	g1 = (unsigned char)((col1&0xFFFF)>>8);
	b1 = (unsigned char)(col1&0xFF);
	r2 = (unsigned char)(col2>>16);
	g2 = (unsigned char)((col2&0xFFFF)>>8);
	b2 = (unsigned char)(col2&0xFF);

	if (hsv) {
		RGB2HSV( r1/255.0, g1/255.0, b1/255.0, h1, s1, v1 );
		RGB2HSV( r2/255.0, g2/255.0, b2/255.0, h2, s2, v2 );
		if (fabs(h2-h1) > 0.5) {
			if (h2 > h1)
				th = (1.0-f)*(h1+1.0)+f*h2;
			else
				th = (1.0-f)*h1+f*(h2+1.0);
			if (th > 1.0)
				th -= 1.0;
		} else
			th = (1.0-f)*h1+f*h2;
		HSV2RGB( th, (1.0-f)*s1+f*s2, (1.0-f)*v1+f*v2, fr, fg, fb );
		r = (unsigned char)floor(fr*255.0);
		g = (unsigned char)floor(fg*255.0);
		b = (unsigned char)floor(fb*255.0);
	} else {
		r = (unsigned char)((1.0-f)*r1+f*r2);
		g = (unsigned char)((1.0-f)*g1+f*g2);
		b = (unsigned char)((1.0-f)*b1+f*b2);
	}
	return (r<<16)+(g<<8)+b;
}



void CircleCenter( double rad, double x1, double y1, double x2, double y2, double &cx, double &cy ) {

	double q, sqt;
	double tx, ty;


	q = sqrt( pow2( x2-x1 ) + pow2( y2-y1 ) );

	sqt = sqrt( pow2(rad) - pow2(q)/4.0 );

	cx = (x1+x2)/2.0 + sqt * (y1 - y2) / q;
	cy = (y1+y2)/2.0 + sqt * (x2 - x1) / q;

	tx = (x1+x2)/2.0 - sqt * (y1 - y2) / q;
	ty = (y1+y2)/2.0 - sqt * (x2 - x1) / q;

	if (sqrt( pow2(tx-500.0) + pow2(ty-500.0) ) > sqrt( pow2(cx-500.0) + pow2(cy-500.0) )) {
		cx = tx;
		cy = ty;
	}
}



void DrawSankeyConnection( CSVGWriter *svg, double mainrad, double a1from, double a1to, double a2from, double a2to, unsigned long col1, unsigned long col2, bool fill, bool smooth, bool hsv ) {

	CSVGObject_Path *svgp;
	double trad1, trad2;
	double a1d, a2d;
	double tcen1x, tcen1y, tcen2x, tcen2y;
	double tcart1fromx, tcart1fromy, tcart1tox, tcart1toy, tcart2fromx, tcart2fromy, tcart2tox, tcart2toy;
	double ta1from, ta1to, ta2from, ta2to;
	double ta1, ta2, ainc1, ainc2;
	double tf, dda;
	int z, steps;


	steps = 200;

	if (fill && smooth) {

		dda = 1.4;

		tcart1fromx = 500.0 + mainrad * cos(a1from);
		tcart1fromy = 500.0 + mainrad * sin(a1from);
		tcart1tox   = 500.0 + mainrad * cos(a1to);
		tcart1toy   = 500.0 + mainrad * sin(a1to);
		tcart2fromx = 500.0 + mainrad * cos(a2from);
		tcart2fromy = 500.0 + mainrad * sin(a2from);
		tcart2tox   = 500.0 + mainrad * cos(a2to);
		tcart2toy   = 500.0 + mainrad * sin(a2to);

		a1d = a1from-a2to;
		if (a1d < Pi)
			a1d += 2.0*Pi;
		if (a1d > Pi)
			a1d -= 2.0*Pi;

		a2d = a2from-a1to;
		if (a2d < Pi)
			a2d += 2.0*Pi;
		if (a2d > Pi)
			a2d -= 2.0*Pi;

		trad1 = mainrad * tan( fabs(a1d/2.0) );
		trad2 = mainrad * tan( fabs(a2d/2.0) );

		CircleCenter( trad1, tcart1fromx, tcart1fromy, tcart2tox, tcart2toy, tcen1x, tcen1y );
		CircleCenter( trad2, tcart2fromx, tcart2fromy, tcart1tox, tcart1toy, tcen2x, tcen2y );

		ta1from = atan2( tcart1fromy-tcen1y, tcart1fromx-tcen1x );
		ta1to   = atan2( tcart2toy-tcen1y,   tcart2tox-tcen1x );
		ta2from = atan2( tcart2fromy-tcen2y, tcart2fromx-tcen2x );
		ta2to   = atan2( tcart1toy-tcen2y,   tcart1tox-tcen2x );

		tf = ta1to - ta1from;
		if (tf > Pi)
			tf -= 2.0*Pi;
		if (tf < -Pi)
			tf += 2.0*Pi;
		ainc1 = tf / steps;

		tf = ta2from - ta2to;
		if (tf > Pi)
			tf -= 2.0*Pi;
		if (tf < -Pi)
			tf += 2.0*Pi;
		ainc2 = tf / steps;

		svg->SetFill(true);
		svg->SetStroke(false);
		svg->SetStrokeColor( 0x000000 );
		svg->SetFillOpacity(1.0);

		ta1 = ta1from;
		ta2 = ta2to;
		for (z=0;z<steps;z++) {

			if ((double)z/(steps-1) < 0.1)
				tf = 0;
			else if ((double)z/(steps-1) > 0.9)
				tf = 1.0;
			else
				tf = (((double)z/(steps-1))-0.1)/0.8;

			svg->SetFillColor( ColorInterpolate( col1, col2, tf, hsv ) );

			svgp = svg->AddPath();

			svgp->AddMoveTo( tcen1x + trad1 * cos(ta1-(dda-1.0)*ainc1), tcen1y + trad1 * sin(ta1-(dda-1.0)*ainc1) );

			svgp->AddArcTo( fabs(trad1), fabs(trad1), 0, 0, (ainc1>0), tcen1x + trad1 * cos(ta1+dda*ainc1), tcen1y + trad1 * sin(ta1+dda*ainc1) );

			if (z == steps-1)
				svgp->AddArcTo( mainrad, mainrad, 0, 0, 0, tcen2x + trad2 * cos(ta2+dda*ainc2), tcen2y + trad2 * sin(ta2+dda*ainc2)  );
			else
				svgp->AddLineTo( tcen2x + trad2 * cos(ta2+dda*ainc2), tcen2y + trad2 * sin(ta2+dda*ainc2)  );

			svgp->AddArcTo( fabs(trad2), fabs(trad2), 0, 0, (ainc2<0), tcen2x + trad2 * cos(ta2-(dda-1.0)*ainc2), tcen2y + trad2 * sin(ta2-(dda-1.0)*ainc2) );

			if (z == 0)
				svgp->AddArcTo( mainrad, mainrad, 0, 0, 0, tcen1x + trad1 * cos(ta1-(dda-1.0)*ainc1), tcen1y + trad1 * sin(ta1-(dda-1.0)*ainc1) );
			else
				svgp->AddLineTo( tcen1x + trad1 * cos(ta1-(dda-1.0)*ainc1), tcen1y + trad1 * sin(ta1-(dda-1.0)*ainc1) );

			svgp->AddClose();

			ta1 += ainc1;
			ta2 += ainc2;
		}

	} else if (fill && !smooth) {

		svg->SetFill( true );
		svg->SetStroke( false );
		svg->SetFillColor( ColorMix( col1, col2 ) );
		svg->SetFillOpacity(1.0);

		svgp = svg->AddPath();

		svgp->AddMoveTo( 500.0 + mainrad * cos(a1from), 500.0 + mainrad * sin(a1from) );

		ArcPart( svgp, mainrad, a1from, a2to );

		svgp->AddArcTo( mainrad, mainrad, 0, 0, 0, 500.0 + mainrad * cos(a2from), 500.0 + mainrad * sin(a2from) );

		ArcPart( svgp, mainrad, a2from, a1to );

		svgp->AddArcTo( mainrad, mainrad, 0, 0, 0, 500.0 + mainrad * cos(a1from), 500.0 + mainrad * sin(a1from) );

		svgp->AddClose();

	} else {

		svg->SetFill( false );
		svg->SetStroke( true );
		svg->SetStrokeColor( 0x000000 );
		svg->SetStrokeWidth( 6.0 );

		svgp = svg->AddPath();

		svgp->AddMoveTo( 500.0 + mainrad * cos(a1from), 500.0 + mainrad * sin(a1from) );

		ArcPart( svgp, mainrad, a1from, a2to );

		svgp = svg->AddPath();

		svgp->AddMoveTo( 500.0 + mainrad * cos(a2from), 500.0 + mainrad * sin(a2from) );

		ArcPart( svgp, mainrad, a2from, a1to );
	}
}



void CSankeyDiagramEngine::WriteSVG( const char *s ) {

	if (m_iType == 0)
		WriteSVG_Block( s );
	else if (m_iType == 1)
		WriteSVG_Circular( s );
}



void CSankeyDiagramEngine::WriteSVG_Circular( const char *s ) {

	CSVGWriter svg;
	int z, zi, z2, z2i, z2r, ti;
	double mainrad;
	double padabs;
	double unit;
	double tf;
	double drad1, drad2, dovl, txtrad;
	CSVGObject_Path *svgp;
	CSankeyNode *sn, *sn2;
	char buf[256];


	mainrad = 400.0;

	drad1 = 5.0;
	drad2 = 20.0;
	dovl = 0.005;
	txtrad = 70.0;


	m_fTotalSum = 0;
	for (z=0;z<m_iNodeCount;z++) {
		zi = m_iaGroupPermutation[z];
		m_oaNodes.push_back( CSankeyNode() );
		m_oaNodes[z].m_faValues.resize( m_iNodeCount );
		m_oaNodes[z].m_iaConnectIndex.resize( m_iNodeCount-1 );
		m_oaNodes[z].m_fSum = 0;
		for (z2=0;z2<m_iNodeCount;z2++) {
			z2i = m_iaGroupPermutation[z2];
			if ((m_faMatrix[ z2i*m_iNodeCount + zi ]+m_faMatrix[ zi*m_iNodeCount + z2i ])/2.0 < m_fZeroThreshold)
				continue;
			m_oaNodes[z].m_faValues[z2] = m_faMatrix[ z2i*m_iNodeCount + zi ];
			m_oaNodes[z].m_fSum += m_faMatrix[ z2i*m_iNodeCount + zi ];
		}
		ti = z;
		for (z2=0;z2<m_iNodeCount-1;z2++) {
			ti--;
			if (ti < 0)
				ti += m_iNodeCount;
			m_oaNodes[z].m_iaConnectIndex[z2] = ti;
		}
		m_fTotalSum += m_oaNodes[z].m_fSum;
		m_oaNodes[z].m_iColor = m_iaColors[zi];
	}

	unit = 2.0 * Pi / ((1.0 + m_fPadding/100.0) * m_fTotalSum);
	padabs = (m_fPadding/100.0 * m_fTotalSum) / m_iNodeCount;

	tf = 0;
	for (z=0;z<m_iNodeCount;z++) {
		m_oaNodes[z].m_fPosStart = tf;
		m_oaNodes[z].m_fPosEnd = tf + m_oaNodes[z].m_fSum * unit;
		tf += (m_oaNodes[z].m_fSum + padabs) * unit;
	}

	for (z=0;z<m_iNodeCount;z++) {
		zi = m_iaGroupPermutation[z];
		sn = &m_oaNodes[z];
		tf = 0;
		sn->m_faSubPos.resize( m_iNodeCount );
		for (z2=0;z2<m_iNodeCount;z2++) {
			sn->m_faSubPos[z2] = sn->m_fPosStart + tf * unit;
			if (z2 < m_iNodeCount-1) {
				z2i = m_iaGroupPermutation[sn->m_iaConnectIndex[z2]];
				//z2i = m_iaGroupPermutation[z2];
				if ((m_faMatrix[ z2i*m_iNodeCount + zi ]+m_faMatrix[ zi*m_iNodeCount + z2i ])/2.0 >= m_fZeroThreshold)
					tf += sn->m_faValues[ sn->m_iaConnectIndex[z2] ];
			}
		}
	}


	// Determine output size
	svg.SetUserSize( 1000.0, 1000.0 );

	svg.SetPixelSize(800);

	// Background Fill
	svg.SetFill(true);
	svg.SetStroke(false);
	svg.SetFillColor(0xE8E8FF);
	svg.SetFillColor(0xFFFFFF);
	//svg.AddRectangle( 10.0, 10.0, 990.0, 990.0 );

	svg.SetFill(false);
	svg.SetStroke(true);
	svg.SetStrokeColor(0x000000);
	svg.SetStrokeWidth(5.0);

	svg.AddCircle( 500.0, 500.0, mainrad );

	svg.SetStroke(true);
	svg.SetStrokeWidth(6.0);
	svg.SetFill(false);

	for (z=0;z<m_iNodeCount;z++) {

		sn = &m_oaNodes[z];

		svgp = svg.AddPath();
		svgp->AddMoveTo( 500.0 + (mainrad-drad1) * cos(sn->m_fPosStart), 500.0 + (mainrad-drad1) * sin(sn->m_fPosStart) );
		svgp->AddArcTo( mainrad-drad1, mainrad-drad1, 0, 0, 1, 500.0 + (mainrad-drad1) * cos(sn->m_fPosEnd), 500.0 + (mainrad-drad1) * sin(sn->m_fPosEnd) );
		svgp->AddLineTo( 500.0 + (mainrad+drad2) * cos(sn->m_fPosEnd), 500.0 + (mainrad+drad2) * sin(sn->m_fPosEnd) );
		svgp->AddArcTo( mainrad+drad2, mainrad+drad2, 0, 0, 0, 500.0 + (mainrad+drad2) * cos(sn->m_fPosStart), 500.0 + (mainrad+drad2) * sin(sn->m_fPosStart) );
		svgp->AddClose();
	}

	svg.SetStroke(false);
	svg.SetFill(true);

	for (z=0;z<m_iNodeCount;z++) {

		sn = &m_oaNodes[z];

		svg.SetFillColor( sn->m_iColor );

		svgp = svg.AddPath();
		svgp->AddMoveTo( 500.0 + (mainrad-drad1) * cos(sn->m_fPosStart), 500.0 + (mainrad-drad1) * sin(sn->m_fPosStart) );
		svgp->AddArcTo( mainrad-drad1, mainrad-drad1, 0, 0, 1, 500.0 + (mainrad-drad1) * cos(sn->m_fPosEnd), 500.0 + (mainrad-drad1) * sin(sn->m_fPosEnd) );
		svgp->AddLineTo( 500.0 + (mainrad+drad2) * cos(sn->m_fPosEnd), 500.0 + (mainrad+drad2) * sin(sn->m_fPosEnd) );
		svgp->AddArcTo( mainrad+drad2, mainrad+drad2, 0, 0, 0, 500.0 + (mainrad+drad2) * cos(sn->m_fPosStart), 500.0 + (mainrad+drad2) * sin(sn->m_fPosStart) );
		svgp->AddClose();
	}

	svg.SetFill(true);
	svg.SetFillColor(0x000000);
	svg.SetStroke(false);
	svg.SetFontSize(50);
	svg.SetFontWeight(SVG_FONTWEIGHT_BOLD);
	svg.SetFontAlign(SVG_FONTALIGN_CENTER);
	svg.SetFontBaseline( SVG_FONTBASELINE_CENTER );

	for (z=0;z<m_iNodeCount;z++) {
		sn = &m_oaNodes[z];
		sprintf(buf,"%d",m_iaGroupPermutation[z]+1);
		svg.AddText( 500.0 + (mainrad+txtrad) * cos((sn->m_fPosStart+sn->m_fPosEnd)/2.0), 500.0 + (mainrad+txtrad) * sin((sn->m_fPosStart+sn->m_fPosEnd)/2.0) , buf );
	}

	for (z=0;z<m_iNodeCount;z++) {

		zi = m_iaGroupPermutation[z];

		sn = &m_oaNodes[z];

		for (z2=0;z2<m_iNodeCount-1;z2++) {

			//z2i = m_iaGroupPermutation[sn->m_iaConnectIndex[z2]];
			z2i = sn->m_iaConnectIndex[z2];

			//if ((m_faMatrix[ z2i*m_iNodeCount + zi ]+m_faMatrix[ zi*m_iNodeCount + z2i ])/2.0 < m_fZeroThreshold)
			if ((m_faMatrix[ m_iaGroupPermutation[z2i]*m_iNodeCount + zi ]+m_faMatrix[ zi*m_iNodeCount + m_iaGroupPermutation[z2i] ])/2.0 < m_fZeroThreshold)
				continue;

			if (z2i < z)
				continue;

			sn2 = &m_oaNodes[z2i];

			for (z2r=0;z2r<m_iNodeCount-1;z2r++)
				if (sn2->m_iaConnectIndex[z2r] == z)
					break;

			DrawSankeyConnection(
				&svg,
				mainrad,
				sn->m_faSubPos[z2]-((z2!=0)?((sn->m_fPosEnd-sn->m_fPosStart)*dovl):0),
				sn->m_faSubPos[z2+1]+((z2!=m_iNodeCount-2)?((sn->m_fPosEnd-sn->m_fPosStart)*dovl):0),
				sn2->m_faSubPos[z2r]-((z2r!=0)?((sn2->m_fPosEnd-sn2->m_fPosStart)*dovl):0),
				sn2->m_faSubPos[z2r+1]+((z2r!=m_iNodeCount-2)?((sn2->m_fPosEnd-sn2->m_fPosStart)*dovl):0),
				sn->m_iColor,
				sn2->m_iColor,
				false,
				false,
				false
			);
		}
	}

	for (z=0;z<m_iNodeCount;z++) {

		zi = m_iaGroupPermutation[z];

		sn = &m_oaNodes[z];

		for (z2=0;z2<m_iNodeCount-1;z2++) {

			z2i = sn->m_iaConnectIndex[z2];
			//z2i = m_iaGroupPermutation[sn->m_iaConnectIndex[z2]];

			if ((m_faMatrix[ m_iaGroupPermutation[z2i]*m_iNodeCount + zi ]+m_faMatrix[ zi*m_iNodeCount + m_iaGroupPermutation[z2i] ])/2.0 < m_fZeroThreshold)
			//if ((m_faMatrix[ z2i*m_iNodeCount + zi ]+m_faMatrix[ zi*m_iNodeCount + z2i ])/2.0 < m_fZeroThreshold)
				continue;

			if (z2i < z)
				continue;

			sn2 = &m_oaNodes[z2i];

			for (z2r=0;z2r<m_iNodeCount-1;z2r++)
				if (sn2->m_iaConnectIndex[z2r] == z)
					break;

			DrawSankeyConnection(
				&svg,
				mainrad,
				sn->m_faSubPos[z2]-((z2!=0)?((sn->m_fPosEnd-sn->m_fPosStart)*dovl):0),
				sn->m_faSubPos[z2+1]+((z2!=m_iNodeCount-2)?((sn->m_fPosEnd-sn->m_fPosStart)*dovl):0),
				sn2->m_faSubPos[z2r]-((z2r!=0)?((sn2->m_fPosEnd-sn2->m_fPosStart)*dovl):0),
				sn2->m_faSubPos[z2r+1]+((z2r!=m_iNodeCount-2)?((sn2->m_fPosEnd-sn2->m_fPosStart)*dovl):0),
				sn->m_iColor,
				sn2->m_iColor,
				true,
				m_bSmooth,
				m_bSmoothHSV
			);
		}
	}

	svg.WriteSVG( s );
}



void DrawBezierConnection( CSVGWriter *svg, double x1, double y1a, double y1b, double x2, double y2a, double y2b, double bezier, unsigned long col1, unsigned long col2, bool correctwidth, bool fill, bool hsv ) {

	int z, nsteps;
	double t1, t2, tc, px1, py1a, py1b, py1c, px2, py2a, py2b, py2c, nr1, nr2, an;
	double ovl;
	CSVGObject_Path *svgp;


	nsteps = 400;
	ovl = 0.5 / nsteps;

	for (z=0;z<nsteps;z++) {

		t1 = (double)z/nsteps - ovl;
		t2 = (double)(z+1)/nsteps + ovl;

		px1  = pow3(1.0-t1)*x1  + 3.0*pow2(1.0-t1)*t1*(x1+bezier) + 3.0*(1.0-t1)*pow2(t1)*(x2-bezier) + pow3(t1)*x2;
		px2  = pow3(1.0-t2)*x1  + 3.0*pow2(1.0-t2)*t2*(x1+bezier) + 3.0*(1.0-t2)*pow2(t2)*(x2-bezier) + pow3(t2)*x2;

		py1a = pow3(1.0-t1)*y1a + 3.0*pow2(1.0-t1)*t1*y1a         + 3.0*(1.0-t1)*pow2(t1)*y2a         + pow3(t1)*y2a;
		py1b = pow3(1.0-t1)*y1b + 3.0*pow2(1.0-t1)*t1*y1b         + 3.0*(1.0-t1)*pow2(t1)*y2b         + pow3(t1)*y2b;

		py2a = pow3(1.0-t2)*y1a + 3.0*pow2(1.0-t2)*t2*y1a         + 3.0*(1.0-t2)*pow2(t2)*y2a         + pow3(t2)*y2a;
		py2b = pow3(1.0-t2)*y1b + 3.0*pow2(1.0-t2)*t2*y1b         + 3.0*(1.0-t2)*pow2(t2)*y2b         + pow3(t2)*y2b;

		if (fill) {
			tc = (x2-((px1+px2)/2.0))/(x2-x1);
			if (tc < 0.1)
				tc = 0;
			else if (tc > 0.9)
				tc = 1.0;
			else
				tc = (tc-0.1)/0.8;
			tc = 1.0 - tc;

			svg->SetFillColor( ColorInterpolate( col1, col2, tc, hsv ) );
		}

		py1c = pow3(1.0-t1)*(y1a+y1b)/2.0 + 3.0*pow2(1.0-t1)*t1*(y1a+y1b)/2.0 + 3.0*(1.0-t1)*pow2(t1)*(y2a+y2b)/2.0 + pow3(t1)*(y2a+y2b)/2.0;
		py2c = pow3(1.0-t2)*(y1a+y1b)/2.0 + 3.0*pow2(1.0-t2)*t2*(y1a+y1b)/2.0 + 3.0*(1.0-t2)*pow2(t2)*(y2a+y2b)/2.0 + pow3(t2)*(y2a+y2b)/2.0;

		nr1 = (py1b-py1a)/2.0;
		nr2 = (py2b-py2a)/2.0;

		if (correctwidth) {
			an = 1.0/cos(atan(fabs(py2c-py1c)/(px2-px1)));

			nr1 *= an;
			nr2 *= an;
		}

		if (fill) {
			svgp = svg->AddPath();
			svgp->AddMoveTo( px1, py1c-nr1 );
			svgp->AddLineTo( px2, py2c-nr2 );
			svgp->AddLineTo( px2, py2c+nr2 );
			svgp->AddLineTo( px1, py1c+nr1 );
			svgp->AddClose();
		} else {
			svgp = svg->AddPath();
			svgp->AddMoveTo( px1, py1c-nr1 );
			svgp->AddLineTo( px2, py2c-nr2 );
			svgp = svg->AddPath();
			svgp->AddMoveTo( px1, py1c+nr1 );
			svgp->AddLineTo( px2, py2c+nr2 );
		}
	}
}



void CSankeyDiagramEngine::WriteSVG_Block( const char *s ) {

	int z, zi, z2, z2i;
	CSVGWriter svg;
	double tf;
	double unitsource, padabssource, unitsink, padabssink;
	double blockheight, blockwidth, blockdist, blockpad, bezier;
	CSankeyNode *sn, *sn2;
	char buf[256];


	blockheight = 900.0;
	blockwidth = 100.0;
	blockdist = 700.0;
	blockpad = 50.0;
	bezier = blockdist / 2.0;


	m_fTotalSumSource = 0;
	for (z=0;z<m_iSourceCount;z++) {
		zi = m_iaGroupPermutationSource[z];
		m_oaNodesSource.push_back( CSankeyNode() );
		m_oaNodesSource[z].m_faValues.resize( m_iSinkCount );
		m_oaNodesSource[z].m_iaConnectIndex.resize( m_iSinkCount );
		m_oaNodesSource[z].m_fSum = 0;
		for (z2=0;z2<m_iSinkCount;z2++) {
			z2i = m_iaGroupPermutationSink[z2];
			if ((m_faMatrixSource[ zi*m_iSinkCount + z2i ]+m_faMatrixSink[ zi*m_iSinkCount + z2i ])/2.0 < m_fZeroThreshold)
				continue;
			m_oaNodesSource[z].m_faValues[z2] = m_faMatrixSource[ zi*m_iSinkCount + z2i ];
			m_oaNodesSource[z].m_fSum += m_faMatrixSource[ zi*m_iSinkCount + z2i ];
		}
		for (z2=0;z2<m_iSinkCount;z2++)
			m_oaNodesSource[z].m_iaConnectIndex[z2] = z2;
		m_fTotalSumSource += m_oaNodesSource[z].m_fSum;
		m_oaNodesSource[z].m_iColor = m_iaColorsSource[zi];
	}

	m_fTotalSumSink = 0;
	for (z=0;z<m_iSinkCount;z++) {
		zi = m_iaGroupPermutationSink[z];
		m_oaNodesSink.push_back( CSankeyNode() );
		m_oaNodesSink[z].m_faValues.resize( m_iSourceCount );
		m_oaNodesSink[z].m_iaConnectIndex.resize( m_iSourceCount );
		m_oaNodesSink[z].m_fSum = 0;
		for (z2=0;z2<m_iSourceCount;z2++) {
			z2i = m_iaGroupPermutationSource[z2];
			if ((m_faMatrixSource[ z2i*m_iSinkCount + zi ]+m_faMatrixSink[ z2i*m_iSinkCount + zi ])/2.0 < m_fZeroThreshold)
				continue;
			m_oaNodesSink[z].m_faValues[z2] = m_faMatrixSink[ z2i*m_iSinkCount + zi ];
			m_oaNodesSink[z].m_fSum += m_faMatrixSink[ z2i*m_iSinkCount + zi ];
		}
		for (z2=0;z2<m_iSourceCount;z2++)
			m_oaNodesSink[z].m_iaConnectIndex[z2] = z2;
		m_fTotalSumSink += m_oaNodesSink[z].m_fSum;
		m_oaNodesSink[z].m_iColor = m_iaColorsSink[zi];
	}


	unitsource = blockheight / ((1.0 + m_fPadding/100.0) * m_fTotalSumSource);
	padabssource = (m_fPadding/100.0 * m_fTotalSumSource) / (m_iSourceCount-1);

	unitsink = blockheight / ((1.0 + m_fPadding/100.0) * m_fTotalSumSink);
	padabssink = (m_fPadding/100.0 * m_fTotalSumSink) / (m_iSinkCount-1);

	if (m_bCorrectRatio) {

		if (unitsource < unitsink) {

			padabssink = (blockheight/unitsource - m_fTotalSumSink) / (m_iSinkCount-1);
			unitsink = unitsource;

		} else { // unitsource > unitsink

			padabssource = (blockheight/unitsink - m_fTotalSumSource) / (m_iSourceCount-1);
			unitsource = unitsink;
		}
	}

	tf = 0;
	for (z=0;z<m_iSourceCount;z++) {
		m_oaNodesSource[z].m_fPosStart = tf;
		m_oaNodesSource[z].m_fPosEnd = tf + m_oaNodesSource[z].m_fSum * unitsource;
		tf += (m_oaNodesSource[z].m_fSum + padabssource) * unitsource;
	}

	for (z=0;z<m_iSourceCount;z++) {
		zi = m_iaGroupPermutationSource[z];
		sn = &m_oaNodesSource[z];
		tf = 0;
		sn->m_faSubPos.resize( m_iSinkCount+1 );
		for (z2=0;z2<m_iSinkCount+1;z2++) {
			sn->m_faSubPos[z2] = sn->m_fPosStart + tf * unitsource;
			if (z2 < m_iSinkCount) {
				z2i = m_iaGroupPermutationSink[sn->m_iaConnectIndex[z2]];
				if ((m_faMatrixSource[ zi*m_iSinkCount + z2i ]+m_faMatrixSink[ zi*m_iSinkCount + z2i ])/2.0 >= m_fZeroThreshold)
					tf += sn->m_faValues[ sn->m_iaConnectIndex[z2] ];
			}
		}
	}

	tf = 0;
	for (z=0;z<m_iSinkCount;z++) {
		m_oaNodesSink[z].m_fPosStart = tf;
		m_oaNodesSink[z].m_fPosEnd = tf + m_oaNodesSink[z].m_fSum * unitsink;
		tf += (m_oaNodesSink[z].m_fSum + padabssink) * unitsink;
	}

	for (z=0;z<m_iSinkCount;z++) {
		zi = m_iaGroupPermutationSink[z];
		sn = &m_oaNodesSink[z];
		tf = 0;
		sn->m_faSubPos.resize( m_iSourceCount+1 );
		for (z2=0;z2<m_iSourceCount+1;z2++) {
			sn->m_faSubPos[z2] = sn->m_fPosStart + tf * unitsink;
			if (z2 < m_iSourceCount) {
				z2i = m_iaGroupPermutationSource[sn->m_iaConnectIndex[z2]];
				if ((m_faMatrixSource[ z2i*m_iSinkCount + zi ]+m_faMatrixSink[ z2i*m_iSinkCount + zi ])/2.0 >= m_fZeroThreshold)
					tf += sn->m_faValues[ sn->m_iaConnectIndex[z2] ];
			}
		}
	}


	// Determine output size
	svg.SetUserSize( 20.0+2.0*blockwidth+blockdist+2.0*blockpad, 20.0+blockheight+2.0*blockpad );

	svg.SetPixelSize(800);

	// Background Fill
/*	svg.SetFill(true);
	svg.SetStroke(false);
	svg.SetFillColor(0xE8E8FF);
	//svg.SetFillColor(0xFFFFFF);
	svg.AddRectangle( 10.0, 10.0, 10.0+2.0*blockwidth+blockdist+2.0*blockpad, 10.0+blockheight+2.0*blockpad );
*/

	svg.SetFill(true);
	svg.SetFillColor(0x000000);
	svg.SetStroke(false);
	svg.SetFontSize(50);
	svg.SetFontWeight(SVG_FONTWEIGHT_BOLD);
	svg.SetFontAlign(SVG_FONTALIGN_CENTER);
	svg.SetFontBaseline( SVG_FONTBASELINE_CENTER );

	for (z=0;z<m_iSourceCount;z++) {
		sn = &m_oaNodesSource[z];
		sprintf(buf,"%d",m_iaGroupPermutationSource[z]+1);
		svg.AddText( 10.0+blockpad-30.0, 10.0+blockpad+(sn->m_fPosStart+sn->m_fPosEnd)/2.0 , buf );
	}

	for (z=0;z<m_iSinkCount;z++) {
		sn = &m_oaNodesSink[z];
		sprintf(buf,"%d",m_iaGroupPermutationSink[z]+1);
		svg.AddText( 10.0+blockpad+2.0*blockwidth+blockdist+30.0, 10.0+blockpad+(sn->m_fPosStart+sn->m_fPosEnd)/2.0 , buf );
	}

	svg.SetStroke(true);
	svg.SetStrokeColor(0x000000);
	svg.SetStrokeWidth(5.0);
	svg.SetFill(false);

	for (z=0;z<m_iSourceCount;z++) {
		sn = &m_oaNodesSource[z];
		svg.AddRectangle( 10.0+blockpad, 10.0+blockpad+sn->m_fPosStart, 10.0+blockpad+blockwidth+5.0, 10.0+blockpad+sn->m_fPosEnd );
	}

	for (z=0;z<m_iSinkCount;z++) {
		sn = &m_oaNodesSink[z];
		svg.AddRectangle( 10.0+blockpad+blockwidth+blockdist-5.0, 10.0+blockpad+sn->m_fPosStart, 10.0+blockpad+2.0*blockwidth+blockdist, 10.0+blockpad+sn->m_fPosEnd );
	}

	for (z=0;z<m_iSourceCount;z++) {

		zi = m_iaGroupPermutationSource[z];

		sn = &m_oaNodesSource[z];

		for (z2=0;z2<m_iSinkCount;z2++) {

			z2i = m_iaGroupPermutationSink[z2];

			sn2 = &m_oaNodesSink[z2];

			if ((m_faMatrixSource[ zi*m_iSinkCount + z2i ]+m_faMatrixSink[ zi*m_iSinkCount + z2i ])/2.0 < m_fZeroThreshold)
				continue;

			DrawBezierConnection(
				&svg,
				10.0+blockpad+blockwidth,
				10.0+blockpad+sn->m_faSubPos[z2],
				10.0+blockpad+sn->m_faSubPos[z2+1],
				10.0+blockpad+blockwidth+blockdist,
				10.0+blockpad+sn2->m_faSubPos[z],
				10.0+blockpad+sn2->m_faSubPos[z+1],
				bezier,
				sn->m_iColor,
				sn2->m_iColor,
				m_bCorrectWidth,
				false,
				true
			);
		}
	}


	svg.SetStroke(false);
	svg.SetFill(true);

	for (z=0;z<m_iSourceCount;z++) {
		sn = &m_oaNodesSource[z];
		svg.SetFillColor( sn->m_iColor );
		svg.AddRectangle( 10.0+blockpad, 10.0+blockpad+sn->m_fPosStart, 10.0+blockpad+blockwidth+5.0, 10.0+blockpad+sn->m_fPosEnd );
	}

	for (z=0;z<m_iSinkCount;z++) {
		sn = &m_oaNodesSink[z];
		svg.SetFillColor( sn->m_iColor );
		svg.AddRectangle( 10.0+blockpad+blockwidth+blockdist-5.0, 10.0+blockpad+sn->m_fPosStart, 10.0+blockpad+2.0*blockwidth+blockdist, 10.0+blockpad+sn->m_fPosEnd );
	}

	for (z=0;z<m_iSourceCount;z++) {

		zi = m_iaGroupPermutationSource[z];

		sn = &m_oaNodesSource[z];

		for (z2=0;z2<m_iSinkCount;z2++) {

			z2i = m_iaGroupPermutationSink[z2];

			sn2 = &m_oaNodesSink[z2];

			if ((m_faMatrixSource[ zi*m_iSinkCount + z2i ]+m_faMatrixSink[ zi*m_iSinkCount + z2i ])/2.0 < m_fZeroThreshold)
				continue;

			DrawBezierConnection(
				&svg,
				10.0+blockpad+blockwidth,
				10.0+blockpad+sn->m_faSubPos[z2],
				10.0+blockpad+sn->m_faSubPos[z2+1],
				10.0+blockpad+blockwidth+blockdist,
				10.0+blockpad+sn2->m_faSubPos[z],
				10.0+blockpad+sn2->m_faSubPos[z+1],
				bezier,
				sn->m_iColor,
				sn2->m_iColor,
				m_bCorrectWidth,
				true,
				true
			);
		}
	}

	svg.WriteSVG( s );
}




