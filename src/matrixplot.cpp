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

#include "matrixplot.h"
#include "svgwriter.h"
#include "xstring.h"


const char *GetRevisionInfo_matrixplot(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_matrixplot() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



CMatrixPlot::CMatrixPlot() {

	m_iRows = 0;
	m_iCols = 0;

	m_bShowLabelTop = false;
	m_bShowLabelBottom = true;
	m_bShowLabelLeft = true;
	m_bShowLabelRight = false;
	m_bRangeSet = false;
	m_fRangeMin = 0;
	m_fRangeMax = 0;

	Set2DMode(false);
	SetValueLabel("Intensity");
	SetValueLabel2("Intensity 2");
	SetInvertColorScale(false);
	SetInvertColorScale2(false);
	SetPlotExponent(1.0);
	SetPlotExponent2(1.0);
	SetRangeZero(false);
	SetRangeZero2(false);
	ResetRangeCut();
	ResetRangeCut2();
}


CMatrixPlot::~CMatrixPlot() {
}


void CMatrixPlot::Init(int rows, int cols) {

	int z;

	m_iRows = rows;
	m_iCols = cols;

	m_saRowLabels.resize(m_iRows);
	m_saColumnLabels.resize(m_iCols);

	m_faBin.resize(m_iRows*m_iCols);
	m_faBin2.resize(m_iRows*m_iCols);
	m_iaActive.resize(m_iRows*m_iCols);
	for (z=0;z<m_iRows*m_iCols;z++) {
		m_faBin[z] = 0;
		m_faBin2[z] = 0;
		m_iaActive[z] = 1;
	}
}


void CMatrixPlot::SetRowLabel(int row, const char *label) {

	m_saRowLabels[row] = label;
}


void CMatrixPlot::SetColumnLabel(int column, const char *label) {

	m_saColumnLabels[column] = label;
}


void CMatrixPlot::SetLabelShow(bool top, bool bottom, bool left, bool right) {

	m_bShowLabelTop = top;
	m_bShowLabelBottom = bottom;
	m_bShowLabelLeft = left;
	m_bShowLabelRight = right;
}


void CMatrixPlot::AddRowCategory(int i) {

	m_iaRowCategories.push_back(i);
}


void CMatrixPlot::AddColumnCategory(int i) {

	m_iaColumnCategories.push_back(i);
}


void CMatrixPlot::AddRowSubCategory(int i) {

	m_iaRowSubCategories.push_back(i);
}


void CMatrixPlot::AddColumnSubCategory(int i) {

	m_iaColumnSubCategories.push_back(i);
}


void CMatrixPlot::SetPlotExponent(double e) {

	m_fPlotExponent = e;
}


void CMatrixPlot::SetPlotExponent2(double e) {

	m_fPlotExponent2 = e;
}


void CMatrixPlot::SetRangeZero(bool b) {

	m_bRangeZero = b;
}


void CMatrixPlot::SetRangeZero2(bool b) {

	m_bRangeZero2 = b;
}


void CMatrixPlot::SetValueLabel(const char *p) {

	m_sValueLabel = p;
}


void CMatrixPlot::SetValueLabel2(const char *p) {

	m_sValueLabel2 = p;
}


void CMatrixPlot::SetInvertColorScale(bool b) {

	m_bInvertColorScale = b;
}


void CMatrixPlot::SetInvertColorScale2(bool b) {

	m_bInvertColorScale2 = b;
}


void CMatrixPlot::Set2DMode(bool b) {

	m_b2DMode = b;
}


void CMatrixPlot::SetRangeCut(double d) {

	m_fRangeCut = d;
}


void CMatrixPlot::SetRangeCut2(double d) {

	m_fRangeCut2 = d;
}


void CMatrixPlot::SetRange(double mi, double ma) {

	m_fRangeMin = mi;
	m_fRangeMax = ma;
	m_bRangeSet = true;
}


void CMatrixPlot::ResetRangeCut() {

	m_fRangeCut = 1.0e80;
}


void CMatrixPlot::ResetRangeCut2() {

	m_fRangeCut2 = 1.0e80;
}


unsigned long CMatrixPlot::ColorFunction(double v) {

	double r, g, b;
	int ir, ig, ib;


	if (v < 0)
		v = 0;
	if (v > 1.0)
		v = 1.0;

	if (m_bInvertColorScale)
		v = 1.0 - v;

	if (v <= 0.0) {
		r = 1.0;
		g = 1.0;
		b = 1.0;
	} else if (v < 0.2) {
		r = 1.0-(mypow(5.0,1.5)*mypow(v,1.5))*0.8;
		g = 1.0-(mypow(5.0,1.5)*mypow(v,1.5))*0.8;
		b = 1.0;
	} else if (v < 0.4) {
		r = 0.2-(v-0.2);
		g = 0.2+0.8*(mypow(5.0,0.75)*mypow(v-0.2,0.75));
		b = 1.0-mypow(5.0,1.33)*mypow(v-0.2,1.33);
	} else if (v < 0.6) {
		r = mypow(5.0,0.5)*mypow(v-0.4,0.5);
		g = 1.0;
		b = 0.0;
	} else if (v < 0.8) {
		r = 1.0;
		g = 1.0-5.0*(v-0.6);
		b = 0.0;
	} else if (v < 1.0) {
		r = 1.0;
		g = 0.0;
		b = 5.0*(v-0.8);
	} else {
		r = 1.0;
		g = 0.0;
		b = 1.0;
	}

	ir = (int)(255.0*r);
	ig = (int)(255.0*g);
	ib = (int)(255.0*b);

	if (ir < 0)
		ir = 0;
	if (ir > 255)
		ir = 255;
	if (ig < 0)
		ig = 0;
	if (ig > 255)
		ig = 255;
	if (ib < 0)
		ib = 0;
	if (ib > 255)
		ib = 255;

	return ir*65536 + ig*256 + ib;
}


unsigned long CMatrixPlot::ColorFunction2D(double x, double y) {

	double r, g, b;
	int ir, ig, ib;


	if (x < 0)
		x = 0;
	if (x > 1.0)
		x = 1.0;
	if (y < 0)
		y = 0;
	if (y > 1.0)
		y = 1.0;

	if (m_bInvertColorScale)
		x = 1.0 - x;
	if (m_bInvertColorScale2)
		y = 1.0 - y;

	b = 1.0-y;
	g = 1.0-x;
	r = 1.0-x + 1.6*mypow(x,1.0)*mypow(y,1.3);

	ir = (int)(255.0*r);
	ig = (int)(255.0*g);
	ib = (int)(255.0*b);

	if (ir < 0)
		ir = 0;
	if (ir > 255)
		ir = 255;
	if (ig < 0)
		ig = 0;
	if (ig > 255)
		ig = 255;
	if (ib < 0)
		ib = 0;
	if (ib > 255)
		ib = 255;

	return ir*65536 + ig*256 + ib;
}


void CMatrixPlot::WriteSVGPlot(const char *s) {

	CSVGWriter sw;
	int z, ix, iy, major, minor, major2, minor2, i;
	CxString buf;
	double rowlabellength, collabellength;
	double fontsize;
	double margin;
	double binsize;
	double ticklength, tickdist;
	double xll, xlr, xplot;
	double ylt, ylb, yplot, xleg, yleg;
	double minval, maxval, minval2, maxval2;


	// Design settings
	if ((m_iRows < 8) || (m_iCols < 8))
		fontsize = 48.0;
	else
		fontsize = 64.0;
	margin = 50.0;
	if ((m_iRows < 5) || (m_iCols < 5))
		binsize = 160.0;
	else if ((m_iRows < 10) || (m_iCols < 10))
		binsize = 120.0;
	else
		binsize = 100.0;
	ticklength = 20.0;
	tickdist = 10.0;


	// Computations
	rowlabellength = 0;
	for (z=0;z<m_iRows;z++)
		if (m_saRowLabels[z].length() > rowlabellength)
			rowlabellength = (double)m_saRowLabels[z].length();
	collabellength = 0;
	for (z=0;z<m_iCols;z++)
		if (m_saColumnLabels[z].length() > collabellength)
			collabellength = (double)m_saColumnLabels[z].length();
	rowlabellength *= 0.7 * fontsize;
	collabellength *= 0.7 * fontsize;
	xplot = binsize * m_iCols;
	yplot = binsize * m_iRows;
	if (m_bShowLabelTop)
		ylt = collabellength + ticklength + tickdist;
	else
		ylt = 0;
	if (m_bShowLabelBottom)
		ylb = collabellength + ticklength + tickdist;
	else
		ylb = 0;
	if (m_bShowLabelLeft)
		xll = rowlabellength + ticklength + tickdist;
	else
		xll = 0;
	if (m_bShowLabelRight)
		xlr = rowlabellength + ticklength + tickdist;
	else
		xlr = 0;
	if (m_b2DMode) {
		xleg = yplot*2.0/3.0 + margin + 20.0 + fontsize*0.7*11 + 20.0;
		yleg = 0;
	} else {
		xleg = 3.0 * fontsize;
		yleg = 120.0 + 2.5*fontsize + margin;
	}


	// Bin value processing
	if (m_bRangeSet) {
		minval = m_fRangeMin;
		maxval = m_fRangeMax;
	} else {
		minval = 1.0e30;
		maxval = -1.0e30;
		for (z=0;z<m_iRows*m_iCols;z++) {
			if (!m_iaActive[z])
				continue;
			if (m_faBin[z] < minval)
				minval = m_faBin[z];
			if (m_faBin[z] > maxval)
				maxval = m_faBin[z];
		}
		if (m_bRangeZero && (minval > 0))
			minval = 0;
		if (maxval > m_fRangeCut)
			maxval = m_fRangeCut;
		if (maxval-minval < 1.0e-8) {
			maxval += 1.0e-8;
			minval -= 1.0e-8;
		}
	}
	minval2 = 1.0e30;
	maxval2 = -1.0e30;
	if (m_b2DMode) {
		for (z=0;z<m_iRows*m_iCols;z++) {
			if (!m_iaActive[z])
				continue;
			if (m_faBin2[z] < minval2)
				minval2 = m_faBin2[z];
			if (m_faBin2[z] > maxval2)
				maxval2 = m_faBin2[z];
		}
		if (m_bRangeZero2 && (minval2 > 0))
			minval2 = 0;
		if (maxval2 > m_fRangeCut2)
			maxval2 = m_fRangeCut2;
		if (maxval2-minval2 < 1.0e-8) {
			maxval2 += 1.0e-8;
			minval2 -= 1.0e-8;
		}
	}


	// Determine output size
	sw.SetUserSize(
		2.0*margin+xll+xlr+xplot+xleg,
		2.0*margin+ylt+ylb+yplot+yleg
	);
	// 1680 x 768
	if ((2.0*margin+xll+xlr+xplot+xleg)/(2.0*margin+ylt+ylb+yplot+yleg) > 2.1875)
		sw.SetPixelSize(1680);
	else
		sw.SetPixelSize((int)(768*(2.0*margin+xll+xlr+xplot+xleg)/(2.0*margin+ylt+ylb+yplot+yleg)));

	// Background Fill
	sw.SetFill(true);
	sw.SetStroke(false);
	//sw.SetFillColor(0xE8E8FF);
	sw.SetFillColor(0xFFFFFF);
	sw.AddRectangle(
		0,
		0,
		2.0*margin+xll+xlr+xplot+xleg,
		2.0*margin+ylt+ylb+yplot+yleg
	);


	// Boxes for the labels (debug)
/*	sw.SetFill(false);
	sw.SetStroke(true);
	sw.SetStrokeColor(0xFF0000);
	sw.SetStrokeWidth(1.0);
	if (m_bShowLabelTop)
		sw.AddRectangle(margin+xll,margin,margin+xll+xplot,margin+collabellength);
	if (m_bShowLabelBottom)
		sw.AddRectangle(margin+xll,margin+ylt+yplot+ticklength+tickdist,margin+xll+xplot,margin+ylt+yplot+ylb);
	if (m_bShowLabelLeft)
		sw.AddRectangle(margin,margin+ylt,margin+rowlabellength,margin+ylt+yplot);
	if (m_bShowLabelRight)
		sw.AddRectangle(margin+xll+xplot+ticklength+tickdist,margin+ylt,margin+xll+xplot+xlr,margin+ylt+yplot);
	sw.AddRectangle(margin+xll,2.0*margin+ylt+yplot+ylb,margin+xll+xplot,2.0*margin+ylt+yplot+ylb+yleg);*/


	// The bin fills
	sw.SetFill(true);
	sw.SetStroke(false);
	for (iy=0;iy<m_iRows;iy++) {
		for (ix=0;ix<m_iCols;ix++) {
			if (!m_iaActive[iy*m_iCols+ix])
				continue;
			if (m_b2DMode) {
				sw.SetFillColor(
					ColorFunction2D(
						mypow( (m_faBin[iy*m_iCols+ix]-minval)/(maxval-minval), m_fPlotExponent ),
						mypow( (m_faBin2[iy*m_iCols+ix]-minval2)/(maxval2-minval2), m_fPlotExponent2 )
					) 
				);
				buf.Format("%.3f | %.3f",m_faBin[iy*m_iCols+ix],m_faBin2[iy*m_iCols+ix]);
			} else {
				sw.SetFillColor( ColorFunction( mypow( (m_faBin[iy*m_iCols+ix]-minval)/(maxval-minval), m_fPlotExponent ) ) );
				buf.Format("%.3f",m_faBin[iy*m_iCols+ix]);
			}
			sw.AddRectangle(
				margin+xll+binsize*ix,
				margin+ylt+binsize*iy,
				margin+xll+binsize*(ix+1),
				margin+ylt+binsize*(iy+1),
				(const char*)buf
			);
		}
	}


	// Inactive bins
	sw.SetFill(false);
	sw.SetStroke(true);
	sw.SetStrokeColor(0x000000);
	sw.SetStrokeWidth(1.5);
	for (iy=0;iy<m_iRows;iy++) {
		for (ix=0;ix<m_iCols;ix++) {
			if (m_iaActive[iy*m_iCols+ix])
				continue;
			sw.AddLine(
				margin+xll+binsize*ix,
				margin+ylt+binsize*iy,
				margin+xll+binsize*(ix+1),
				margin+ylt+binsize*(iy+1)
			);
			sw.AddLine(
				margin+xll+binsize*(ix+1),
				margin+ylt+binsize*iy,
				margin+xll+binsize*ix,
				margin+ylt+binsize*(iy+1)
			);
		}
	}


	// Plot Outline Stroke
	sw.SetFill(false);
	sw.SetStroke(true);
	sw.SetStrokeWidth(10.0);
	sw.SetStrokeColor(0x000000);
	sw.AddRectangle(
		margin+xll,
		margin+ylt,
		margin+xll+xplot,
		margin+ylt+yplot
	);


	// Grid Lines
	sw.SetStrokeWidth(1.5);
	for (ix=1;ix<m_iCols;ix++)
		sw.AddLine(
			margin+xll+binsize*ix,
			margin+ylt,
			margin+xll+binsize*ix,
			margin+ylt+yplot
		);
	for (ix=1;ix<m_iRows;ix++)
		sw.AddLine(
			margin+xll,
			margin+ylt+binsize*ix,
			margin+xll+xplot,
			margin+ylt+binsize*ix
		);


	// Category Lines
	sw.SetStrokeWidth(14.0);
	for (ix=0;ix<(int)m_iaColumnCategories.size();ix++)
		sw.AddLine(
			margin+xll+binsize*(m_iaColumnCategories[ix]+1),
			margin+ylt,
			margin+xll+binsize*(m_iaColumnCategories[ix]+1),
			margin+ylt+yplot
		);
	for (ix=0;ix<(int)m_iaRowCategories.size();ix++)
		sw.AddLine(
			margin+xll,
			margin+ylt+binsize*(m_iaRowCategories[ix]+1),
			margin+xll+xplot,
			margin+ylt+binsize*(m_iaRowCategories[ix]+1)
		);


	// Subcategory Lines
	sw.SetStrokeWidth(6.0);
	for (ix=0;ix<(int)m_iaColumnSubCategories.size();ix++)
		sw.AddLine(
			margin+xll+binsize*(m_iaColumnSubCategories[ix]+1),
			margin+ylt,
			margin+xll+binsize*(m_iaColumnSubCategories[ix]+1),
			margin+ylt+yplot
		);
	for (ix=0;ix<(int)m_iaRowSubCategories.size();ix++)
		sw.AddLine(
			margin+xll,
			margin+ylt+binsize*(m_iaRowSubCategories[ix]+1),
			margin+xll+xplot,
			margin+ylt+binsize*(m_iaRowSubCategories[ix]+1)
		);


	// Ticks
	sw.SetStrokeWidth(6.0);
	for (ix=0;ix<m_iCols;ix++) {
		if (m_saColumnLabels[ix].length() == 0)
			continue;
		if (m_bShowLabelTop)
			sw.AddLine(
				margin+xll+binsize*(ix+0.5),
				margin+ylt-ticklength,
				margin+xll+binsize*(ix+0.5),
				margin+ylt
			);
		if (m_bShowLabelBottom)
			sw.AddLine(
				margin+xll+binsize*(ix+0.5),
				margin+ylt+yplot,
				margin+xll+binsize*(ix+0.5),
				margin+ylt+yplot+ticklength
			);
	}
	for (ix=0;ix<m_iRows;ix++) {
		if (m_saRowLabels[ix].length() == 0)
			continue;
		if (m_bShowLabelLeft)
			sw.AddLine(
				margin+xll-ticklength,
				margin+ylt+binsize*(ix+0.5),
				margin+xll,
				margin+ylt+binsize*(ix+0.5)
			);
		if (m_bShowLabelRight)
			sw.AddLine(
			margin+xll+xplot,
			margin+ylt+binsize*(ix+0.5),
			margin+xll+xplot+ticklength,
			margin+ylt+binsize*(ix+0.5)
		);
	}


	// Tick Labels
	sw.SetFill(true);
	sw.SetFillColor(0x000000);
	sw.SetStroke(false);
	sw.SetFontSize(fontsize);
	sw.SetFontWeight(SVG_FONTWEIGHT_BOLD);
	for (ix=0;ix<m_iCols;ix++) {
		if (m_saColumnLabels[ix].length() == 0)
			continue;
		if (m_bShowLabelTop) {
			sw.SetFontAlign(SVG_FONTALIGN_LEFT);
			sw.AddText(
				margin+xll+binsize*(ix+0.5)+fontsize*0.325,
				margin+collabellength,
				m_saColumnLabels[ix].c_str(),
				-90.0
			);
		}
		if (m_bShowLabelBottom) {
			sw.SetFontAlign(SVG_FONTALIGN_RIGHT);
			sw.AddText(
				margin+xll+binsize*(ix+0.5)+fontsize*0.325,
				margin+ylt+yplot+ticklength+tickdist,
				m_saColumnLabels[ix].c_str(),
				-90.0
			);
		}
	}
	for (ix=0;ix<m_iRows;ix++) {
		if (m_saRowLabels[ix].length() == 0)
			continue;
		if (m_bShowLabelLeft) {
			sw.SetFontAlign(SVG_FONTALIGN_RIGHT);
			sw.AddText(
				margin+rowlabellength,
				margin+ylt+binsize*(ix+0.5)+fontsize*0.325,
				m_saRowLabels[ix].c_str()
			);
		}
		if (m_bShowLabelRight) {
			sw.SetFontAlign(SVG_FONTALIGN_LEFT);
			sw.AddText(
				margin+xll+xplot+ticklength+tickdist,
				margin+ylt+binsize*(ix+0.5)+fontsize*0.325,
				m_saRowLabels[ix].c_str()
			);
		}
	}


	if (m_b2DMode) {

		// Legend Fill
		CreateTicks(minval,maxval,major,minor,false);
		CreateTicks(minval2,maxval2,major2,minor2,false);
		sw.SetFill(true);
		sw.SetStroke(false);

		if (m_iRows < 5) {
			if (major > 4)
				major = 4;
			if (major2 > 4)
				major2 = 4;
		}

		for (iy=0;iy<100;iy++) {
			for (ix=0;ix<100;ix++) {
				if ((ix == 99) || (iy == 99))
					i = 1;
				else
					i = 2;
				sw.SetFillColor( ColorFunction2D( mypow(ix/99.0,m_fPlotExponent), mypow(1.0-iy/99.0,m_fPlotExponent2) ) );
				buf.Format("%.3f | %.3f",minval+ix/99.0*(maxval-minval),minval2+iy/99.0*(maxval2-minval2));
				sw.AddRectangle(
					2.0*margin+xll+xplot+xlr+40.0+fontsize*0.7*9+yplot/3.0*2.0*(ix/100.0),
					margin+ylt+yplot/6.0+yplot/3.0*2.0*(iy/100.0),
					2.0*margin+xll+xplot+xlr+40.0+fontsize*0.7*9+yplot/3.0*2.0*((ix+i)/100.0),
					margin+ylt+yplot/6.0+yplot/3.0*2.0*((iy+i)/100.0),
					(const char*)buf
				);
			}
		}


		// Legend Frame
		sw.SetFill(false);
		sw.SetStroke(true);
		sw.SetStrokeWidth(4.0);
		sw.SetStrokeColor(0x000000);
		sw.AddRectangle(
			2.0*margin+xll+xplot+xlr+40.0+fontsize*0.7*9,
			margin+ylt+yplot/6.0,
			2.0*margin+xll+xplot+xlr+40.0+fontsize*0.7*9+yplot/3.0*2.0,
			margin+ylt+5.0*yplot/6.0
		);


		// Legend Ticks
		sw.SetStrokeWidth(4.0);
		for (ix=0;ix<major;ix++)
			sw.AddLine(
				2.0*margin+xll+xplot+xlr+40.0+fontsize*0.7*9+yplot/3.0*2.0*ix/(major-1.0),
				margin+ylt+5.0*yplot/6.0,
				2.0*margin+xll+xplot+xlr+40.0+fontsize*0.7*9+yplot/3.0*2.0*ix/(major-1.0),
				margin+ylt+5.0*yplot/6.0+ticklength
			);
		for (ix=0;ix<major2;ix++)
			sw.AddLine(
				2.0*margin+xll+xplot+xlr+40.0+fontsize*0.7*9-ticklength,
				margin+ylt+yplot/6.0+yplot/3.0*2.0*ix/(major2-1.0),
				2.0*margin+xll+xplot+xlr+40.0+fontsize*0.7*9,
				margin+ylt+yplot/6.0+yplot/3.0*2.0*ix/(major2-1.0)
			);


		// Legend Labels
		sw.SetFill(true);
		sw.SetFillColor(0x000000);
		sw.SetStroke(false);
		sw.SetFontSize(fontsize);
		sw.SetFontWeight(SVG_FONTWEIGHT_BOLD);
		sw.SetFontAlign(SVG_FONTALIGN_CENTER);
		for (ix=0;ix<major;ix++) {
			if (maxval >= 100.0)
				buf.Format("%.0f",minval+((double)ix/(major-1.0))*(maxval-minval));
			else if (maxval >= 10.0)
				buf.Format("%.1f",minval+((double)ix/(major-1.0))*(maxval-minval));
			else
				buf.Format("%.2f",minval+((double)ix/(major-1.0))*(maxval-minval));
			sw.AddText(
				2.0*margin+xll+xplot+xlr+40.0+fontsize*0.7*9+yplot/3.0*2.0*ix/(major-1.0),
				margin+ylt+5.0*yplot/6.0+ticklength+tickdist+fontsize,
				(const char*)buf
			);
		}
		sw.AddText(
			2.0*margin+xll+xplot+xlr+40.0+fontsize*0.7*9+yplot/3.0,
			margin+ylt+5.0*yplot/6.0+ticklength+tickdist+fontsize+1.5*fontsize,
			m_sValueLabel.c_str()
		);

		sw.SetFontAlign(SVG_FONTALIGN_RIGHT);
		for (ix=0;ix<major2;ix++) {
			if (maxval2 >= 100.0)
				buf.Format("%.0f",minval2+(1.0-(double)ix/(major2-1.0))*(maxval2-minval2));
			else if (maxval2 >= 10.0)
				buf.Format("%.1f",minval2+(1.0-(double)ix/(major2-1.0))*(maxval2-minval2));
			else
				buf.Format("%.2f",minval2+(1.0-(double)ix/(major2-1.0))*(maxval2-minval2));
			sw.AddText(
				2.0*margin+xll+xplot+xlr+40.0+fontsize*0.7*9-ticklength-tickdist,
				margin+ylt+yplot/6.0+yplot/3.0*2.0*ix/(major2-1.0)+fontsize*0.325,
				(const char*)buf
			);
		}
		sw.SetFontAlign(SVG_FONTALIGN_CENTER);
		sw.AddText(
			2.0*margin+xll+xplot+xlr+40.0+fontsize*0.7*9-ticklength-tickdist-5.0*0.7*fontsize,
			margin+ylt+yplot/2.0,
			m_sValueLabel2.c_str(),
			-90.0
		);

	} else {

		// Legend Fill
		CreateTicks(minval,maxval,major,minor,false);
		if (m_iCols < 5)
			if (major > 4)
				major = 4;
		sw.SetFill(true);
		sw.SetStroke(false);
		for (z=0;z<400;z++) {
			if (z == 399)
				i = 1;
			else
				i = 2;
			sw.SetFillColor( ColorFunction( mypow(z/399.0,m_fPlotExponent) ) );
			buf.Format("%.3f",minval+z/399.0*(maxval-minval));
			sw.AddRectangle(
				margin+xll+xplot*(z/400.0),
				2.0*margin+10.0+ylt+yplot+ylb,
				margin+xll+xplot*((z+i)/400.0),
				2.0*margin+10.0+ylt+yplot+ylb+50.0,
				(const char*)buf
			);
		}


		// Legend Frame
		sw.SetFill(false);
		sw.SetStroke(true);
		sw.SetStrokeWidth(4.0);
		sw.SetStrokeColor(0x000000);
		sw.AddRectangle(
			margin+xll,
			2.0*margin+10.0+ylt+yplot+ylb,
			margin+xll+xplot,
			2.0*margin+10.0+ylt+yplot+ylb+50.0
		);


		// Legend Ticks
		sw.SetStrokeWidth(4.0);
		for (ix=0;ix<major;ix++)
			sw.AddLine(
				margin+xll+(xplot/(major-1.0))*ix,
				2.0*margin+10.0+ylt+yplot+ylb+50.0,
				margin+xll+(xplot/(major-1.0))*ix,
				2.0*margin+10.0+ylt+yplot+ylb+50.0+ticklength
			);


		// Legend Labels
		sw.SetFill(true);
		sw.SetFillColor(0x000000);
		sw.SetStroke(false);
		sw.SetFontSize(fontsize);
		sw.SetFontWeight(SVG_FONTWEIGHT_BOLD);
		sw.SetFontAlign(SVG_FONTALIGN_CENTER);
		for (ix=0;ix<major;ix++) {
			//buf.Format("%.2f",minval+(mypow((double)ix/major,m_fPlotExponent))*(maxval-minval));
			if (maxval >= 100.0)
				buf.Format("%.0f",minval+((double)ix/(major-1.0))*(maxval-minval));
			else if (maxval >= 10.0)
				buf.Format("%.1f",minval+((double)ix/(major-1.0))*(maxval-minval));
			else
				buf.Format("%.2f",minval+((double)ix/(major-1.0))*(maxval-minval));
			sw.AddText(
				margin+xll+(xplot/(major-1.0))*ix,
				2.0*margin+10.0+ylt+yplot+ylb+50.0+ticklength+tickdist+fontsize,
				(const char*)buf
			);
		}
		sw.AddText(
			margin+xll+xplot/2.0,
			2.0*margin+10.0+ylt+yplot+ylb+50.0+ticklength+tickdist+fontsize+1.5*fontsize,
			m_sValueLabel.c_str()
		);
	}


	// Final Output
	sw.WriteSVG(s);
}





