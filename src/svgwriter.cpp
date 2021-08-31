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

#include "svgwriter.h"


const char *GetRevisionInfo_svgwriter(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_svgwriter() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



CSVGWriter::CSVGWriter() {

	m_fUserWidth = 1.0;
	m_fUserHeight = 1.0;
	SetPixelSize(768);

	SetFill(true);
	SetFillColor(0x000000);
	SetFillOpacity(1.0);
	SetStroke(false);
	SetStrokeWidth(1.0);
	SetStrokeColor(0x000000);
	SetStrokeOpacity(1.0);
	SetStrokeDash(false);
	SetStrokeDashArray(1,0);
	SetStrokeDashOffset(0);

	SetLineCap(SVG_LINECAP_BUTT);
	SetLineJoin(SVG_LINEJOIN_MITER);

	SetFontFamily("sans-serif");
	SetFontSize(10.0);
	SetFontWeight(SVG_FONTWEIGHT_NORMAL);
	SetFontStyle(SVG_FONTSTYLE_NORMAL);
	SetFontStretch(SVG_FONTSTRETCH_NORMAL);
	SetFontAlign(SVG_FONTALIGN_LEFT);
	SetFontBaseline(SVG_FONTBASELINE_BASELINE);

	ResetClipPolygon();
}


CSVGWriter::~CSVGWriter() {
}


void CSVGWriter::SetUserSize(double width, double height) {

	m_fUserWidth = width;
	m_fUserHeight = height;
	m_iPixelHeight = (int)floor((double)m_iPixelWidth * m_fUserHeight / m_fUserWidth + 0.5);
}


void CSVGWriter::SetPixelSize(int width) {

	m_iPixelWidth = width;
	m_iPixelHeight = (int)floor((double)width * m_fUserHeight / m_fUserWidth + 0.5);
}


bool CSVGWriter::WriteSVG(const char *s) {

	FILE *a;
	int z, z2;


	a = OpenFileWrite(s,true);

	mfprintf(a,"<?xml version=\"1.0\" standalone=\"no\" ?>\n");
	mfprintf(a,"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20010904//EN\"\n");
	mfprintf(a,"  \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n");
	mfprintf(a,"<svg width=\"%d\" height=\"%d\" viewBox=\"0 0 %.6f %.6f\" xmlns=\"http://www.w3.org/2000/svg\"\n",
		m_iPixelWidth,m_iPixelHeight,m_fUserWidth,m_fUserHeight);
	mfprintf(a,"  xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n");

	if (m_oaClipPolygons.size() != 0) {
		mfprintf(a,"<defs>\n");
		for (z=0;z<(int)m_oaClipPolygons.size();z++) {
			mfprintf(a,"  <clipPath id=\"myclip%d\">\n",z+1);
			mfprintf(a,"    <polygon points=\"");
			for (z2=0;z2<(int)m_oaClipPolygons[z]->m_faPoints.size();z2++) {
				mfprintf(a,"%.6f",m_oaClipPolygons[z]->m_faPoints[z2]);
				if ((z2+1 < (int)m_oaClipPolygons[z]->m_faPoints.size()) && ((z2 % 2) == 0))
					mfprintf(a,", ");
				else if (z2+1 < (int)m_oaClipPolygons[z]->m_faPoints.size())
					mfprintf(a," ");
			}
			mfprintf(a,"\" />\n");
			mfprintf(a,"  </clipPath>\n");
		}
		mfprintf(a,"</defs>\n");
	}

	for (z=0;z<(int)m_oaObjects.size();z++)
		m_oaObjects[z]->Print(a);

	mfprintf(a,"</svg>\n");

	fclose(a);

	return true;
}


void CSVGWriter::AddRectangle(double x1, double y1, double x2, double y2) {

	AddRectangle(x1,y1,x2,y2,"");
}

	
void CSVGWriter::AddRectangle(double x1, double y1, double x2, double y2, const char *title) {

	CSVGObject_Rectangle *o;

	o = new CSVGObject_Rectangle();
	o->m_fX1 = x1;
	o->m_fY1 = y1;
	o->m_fX2 = x2;
	o->m_fY2 = y2;
	o->m_fRX = 0;
	o->m_fRY = 0;
	o->m_bStroke = (m_iaStroke.back()!=0);
	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_bFill = (m_iaFill.back()!=0);
	o->m_iFillColor = m_iaFillColor.back();
	o->m_fFillOpacity = m_faFillOpacity.back();
	o->m_iLineJoin = m_iaLineJoin.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	o->m_sTitle = title;
	m_oaObjects.push_back(o);
}


void CSVGWriter::AddRectangle(double x1, double y1, double x2, double y2, double radius) {

	AddRectangle(x1,y1,x2,y2,radius,"");
}


void CSVGWriter::AddRectangle(double x1, double y1, double x2, double y2, double radius, const char *title) {

	CSVGObject_Rectangle *o;

	o = new CSVGObject_Rectangle();
	o->m_fX1 = x1;
	o->m_fY1 = y1;
	o->m_fX2 = x2;
	o->m_fY2 = y2;
	o->m_fRX = radius;
	o->m_fRY = radius;
	o->m_bStroke = (m_iaStroke.back()!=0);
	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_bFill = (m_iaFill.back()!=0);
	o->m_iFillColor = m_iaFillColor.back();
	o->m_fFillOpacity = m_faFillOpacity.back();
	o->m_iLineJoin = m_iaLineJoin.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	o->m_sTitle = title;
	m_oaObjects.push_back(o);
}


void CSVGWriter::AddRectangle(double x1, double y1, double x2, double y2, double radx, double rady) {

	AddRectangle(x1,y1,x2,y2,radx,rady,"");
}


void CSVGWriter::AddRectangle(double x1, double y1, double x2, double y2, double radx, double rady, const char *title) {

	CSVGObject_Rectangle *o;

	o = new CSVGObject_Rectangle();
	o->m_fX1 = x1;
	o->m_fY1 = y1;
	o->m_fX2 = x2;
	o->m_fY2 = y2;
	o->m_fRX = radx;
	o->m_fRY = rady;
	o->m_bStroke = (m_iaStroke.back()!=0);
	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_bFill = (m_iaFill.back()!=0);
	o->m_iFillColor = m_iaFillColor.back();
	o->m_fFillOpacity = m_faFillOpacity.back();
	o->m_iLineJoin = m_iaLineJoin.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	o->m_sTitle = title;
	m_oaObjects.push_back(o);
}


void CSVGWriter::AddLine(double x1, double y1, double x2, double y2) {

	AddLine(x1,y1,x2,y2,"");
}


void CSVGWriter::AddLine(double x1, double y1, double x2, double y2, const char *title) {

	CSVGObject_Line *o;

	o = new CSVGObject_Line();
	o->m_fX1 = x1;
	o->m_fY1 = y1;
	o->m_fX2 = x2;
	o->m_fY2 = y2;
	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_iLineCap = m_iaLineCap.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	o->m_sTitle = title;
	m_oaObjects.push_back(o);
}


void CSVGWriter::AddCircle(double cx, double cy, double r) {

	AddCircle(cx,cy,r,"");
}


void CSVGWriter::AddCircle(double cx, double cy, double r, const char *title) {

	AddEllipse(cx,cy,r,r,title);
}


void CSVGWriter::AddEllipse(double cx, double cy, double rx, double ry) {

	AddEllipse(cx,cy,rx,ry,"");
}


void CSVGWriter::AddEllipse(double cx, double cy, double rx, double ry, const char *title) {

	CSVGObject_Ellipse *o;

	o = new CSVGObject_Ellipse();
	o->m_fCX = cx;
	o->m_fCY = cy;
	o->m_fRX = rx;
	o->m_fRY = ry;
	o->m_bStroke = (m_iaStroke.back()!=0);
	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_bFill = (m_iaFill.back()!=0);
	o->m_iFillColor = m_iaFillColor.back();
	o->m_fFillOpacity = m_faFillOpacity.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	o->m_sTitle = title;
	m_oaObjects.push_back(o);
}


void CSVGWriter::AddPolyLine_Variadic(int n, ...) {

	CSVGObject_PolyLine *o;
	int z;
	va_list vl;

	o = new CSVGObject_PolyLine();

	va_start(vl,n);
	for (z=0;z<n;z++)
		o->m_faPoints.push_back( va_arg(vl,double) );
	va_end(vl);

	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_iLineJoin = m_iaLineJoin.back();
	o->m_iLineCap = m_iaLineCap.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	o->m_sTitle = "";
	m_oaObjects.push_back(o);
}


void CSVGWriter::AddPolyLine_Array(int n, double *pts) {

	CSVGObject_PolyLine *o;
	int z;

	o = new CSVGObject_PolyLine();

	for (z=0;z<n;z++)
		o->m_faPoints.push_back( pts[z] );

	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_iLineJoin = m_iaLineJoin.back();
	o->m_iLineCap = m_iaLineCap.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	o->m_sTitle = "";
	m_oaObjects.push_back(o);
}


void CSVGWriter::AddPolyLine_Variadic(const char *title, int n, ...) {

	CSVGObject_PolyLine *o;
	int z;
	va_list vl;

	o = new CSVGObject_PolyLine();

	va_start(vl,n);
	for (z=0;z<n;z++)
		o->m_faPoints.push_back( va_arg(vl,double) );
	va_end(vl);

	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_iLineJoin = m_iaLineJoin.back();
	o->m_iLineCap = m_iaLineCap.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	o->m_sTitle = title;
	m_oaObjects.push_back(o);
}


void CSVGWriter::AddPolygon_Variadic(int n, ...) {

	CSVGObject_Polygon *o;
	int z;
	va_list vl;

	o = new CSVGObject_Polygon();

	va_start(vl,n);
	for (z=0;z<n;z++)
		o->m_faPoints.push_back( va_arg(vl,double) );
	va_end(vl);

	o->m_bStroke = (m_iaStroke.back()!=0);
	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_bFill = (m_iaFill.back()!=0);
	o->m_iFillColor = m_iaFillColor.back();
	o->m_fFillOpacity = m_faFillOpacity.back();
	o->m_iLineJoin = m_iaLineJoin.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	o->m_sTitle = "";
	m_oaObjects.push_back(o);
}


void CSVGWriter::AddPolygon_Array(int n, double *pts) {

	CSVGObject_Polygon *o;
	int z;

	o = new CSVGObject_Polygon();

	for (z=0;z<n;z++)
		o->m_faPoints.push_back( pts[z] );

	o->m_bStroke = (m_iaStroke.back()!=0);
	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_bFill = (m_iaFill.back()!=0);
	o->m_iFillColor = m_iaFillColor.back();
	o->m_fFillOpacity = m_faFillOpacity.back();
	o->m_iLineJoin = m_iaLineJoin.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	o->m_sTitle = "";
	m_oaObjects.push_back(o);
}


void CSVGWriter::AddPolygon_Array(const std::vector<double> *vec) {

	CSVGObject_Polygon *o;
	int z;

	o = new CSVGObject_Polygon();

	for (z=0;z<(int)vec->size();z++)
		o->m_faPoints.push_back( (*vec)[z] );

	o->m_bStroke = (m_iaStroke.back()!=0);
	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_bFill = (m_iaFill.back()!=0);
	o->m_iFillColor = m_iaFillColor.back();
	o->m_fFillOpacity = m_faFillOpacity.back();
	o->m_iLineJoin = m_iaLineJoin.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	o->m_sTitle = "";
	m_oaObjects.push_back(o);
}


void CSVGWriter::AddPolygon_Variadic(const char *title, int n, ...) {

	CSVGObject_Polygon *o;
	int z;
	va_list vl;

	o = new CSVGObject_Polygon();

	va_start(vl,n);
	for (z=0;z<n;z++)
		o->m_faPoints.push_back( va_arg(vl,double) );
	va_end(vl);

	o->m_bStroke = (m_iaStroke.back()!=0);
	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_bFill = (m_iaFill.back()!=0);
	o->m_iFillColor = m_iaFillColor.back();
	o->m_fFillOpacity = m_faFillOpacity.back();
	o->m_iLineJoin = m_iaLineJoin.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	o->m_sTitle = title;
	m_oaObjects.push_back(o);
}


void CSVGWriter::AddText(double x, double y, const char *text) {

	AddText(x,y,text,"");
}


void CSVGWriter::AddText(double x, double y, const char *text, const char *title) {

	CSVGObject_Text *o;

	o = new CSVGObject_Text();
	o->m_fPosX = x;
	o->m_fPosY = y;
	o->m_sText = text;
	o->m_fAngle = 0;
	o->m_bStroke = (m_iaStroke.back()!=0);
	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_bFill = (m_iaFill.back()!=0);
	o->m_iFillColor = m_iaFillColor.back();
	o->m_fFillOpacity = m_faFillOpacity.back();
	o->m_iLineJoin = m_iaLineJoin.back();
	o->m_sFontFamily = m_saFontFamily.back();
	o->m_fFontSize = m_faFontSize.back();
	o->m_iFontWeight = m_iaFontWeight.back();
	o->m_iFontStyle = m_iaFontStyle.back();
	o->m_iFontStretch = m_iaFontStretch.back();
	o->m_iFontAlign = m_iaFontAlign.back();
	o->m_iFontBaseline = m_iaFontBaseline.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	o->m_sTitle = title;
	m_oaObjects.push_back(o);
}


void CSVGWriter::AddText(double x, double y, const char *text, double angle) {

	AddText(x,y,text,angle,"");
}


void CSVGWriter::AddText(double x, double y, const char *text, double angle, const char *title) {

	CSVGObject_Text *o;

	o = new CSVGObject_Text();
	o->m_fPosX = x;
	o->m_fPosY = y;
	o->m_sText = text;
	o->m_fAngle = angle;
	o->m_bStroke = (m_iaStroke.back()!=0);
	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_bFill = (m_iaFill.back()!=0);
	o->m_iFillColor = m_iaFillColor.back();
	o->m_fFillOpacity = m_faFillOpacity.back();
	o->m_iLineJoin = m_iaLineJoin.back();
	o->m_sFontFamily = m_saFontFamily.back();
	o->m_fFontSize = m_faFontSize.back();
	o->m_iFontWeight = m_iaFontWeight.back();
	o->m_iFontStyle = m_iaFontStyle.back();
	o->m_iFontStretch = m_iaFontStretch.back();
	o->m_iFontAlign = m_iaFontAlign.back();
	o->m_iFontBaseline = m_iaFontBaseline.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	o->m_sTitle = title;
	m_oaObjects.push_back(o);
}


CSVGObject_Path* CSVGWriter::AddPath() {

	CSVGObject_Path *o;

	o = new CSVGObject_Path();
	o->m_bStroke = (m_iaStroke.back()!=0);
	o->m_fStrokeWidth = m_faStrokeWidth.back();
	o->m_iStrokeColor = m_iaStrokeColor.back();
	o->m_fStrokeOpacity = m_faStrokeOpacity.back();
	o->m_bStrokeDash = (m_iaStrokeDash.back()!=0);
	o->m_faStrokeDashArray.assign( m_faaStrokeDashArray.back().begin(), m_faaStrokeDashArray.back().end() );
	o->m_fStrokeDashOffset = m_faStrokeDashOffset.back();
	o->m_bFill = (m_iaFill.back()!=0);
	o->m_iFillColor = m_iaFillColor.back();
	o->m_fFillOpacity = m_faFillOpacity.back();
	o->m_iLineJoin = m_iaLineJoin.back();
	o->m_iClipPath = m_iaClipPolygon.back();
	m_oaObjects.push_back(o);

	return o;
}


/*****************************************************************************************************************************/



void CSVGObject_Rectangle::Print(FILE *a) {

	int z;

	mfprintf(a,"  <rect x=\"%.6f\" y=\"%.6f\" width=\"%.6f\" height=\"%.6f\"",
		m_fX1, m_fY1, m_fX2-m_fX1, m_fY2-m_fY1 );

	if ((m_fRX != 0) || (m_fRY != 0))
		mfprintf(a," rx=\"%.6f\" ry=\"%.6f\"",m_fRX,m_fRY);

	if (m_bFill) {
		mfprintf(a," fill=\"#%06lX\"",m_iFillColor);
		if (m_fFillOpacity != 1)
			mfprintf(a," fill-opacity=\"%.6f\"",m_fFillOpacity);
	} else
		mfprintf(a," fill=\"none\"");

	if (m_bStroke) {
		mfprintf(a," stroke=\"#%06lX\" stroke-width=\"%.6f\"",m_iStrokeColor,m_fStrokeWidth);
		if (m_fStrokeOpacity != 1)
			mfprintf(a," stroke-opacity=\"%.6f\"",m_fStrokeOpacity);
		if (m_bStrokeDash) {
			mfprintf(a," stroke-dasharray=\"");
			for (z=0;z<(int)m_faStrokeDashArray.size();z++) {
				mfprintf(a,"%.6f",m_faStrokeDashArray[z]);
				if (z+1 < (int)m_faStrokeDashArray.size())
					mfprintf(a,",");
			}
			mfprintf(a,"\" stroke-dashoffset=\"%.6f\"",m_fStrokeDashOffset);
		}
		if (m_iLineJoin != SVG_LINEJOIN_MITER) {
			switch(m_iLineJoin) {
				case SVG_LINEJOIN_MITER:
					mfprintf(a," stroke-linejoin=\"miter\"");
					break;
				case SVG_LINEJOIN_ROUND:
					mfprintf(a," stroke-linejoin=\"round\"");
					break;
				case SVG_LINEJOIN_BEVEL:
					mfprintf(a," stroke-linejoin=\"bevel\"");
					break;
			}
		}
	} else
		mfprintf(a," stroke=\"none\"");

	if (m_iClipPath != -1)
		mfprintf(a," clip-path=\"url(#myclip%d)\"",m_iClipPath+1);

	if (m_sTitle.length() != 0)
		mfprintf(a,"><title>%s</title></rect>\n",m_sTitle.c_str());
	else
		mfprintf(a," />\n");
}


void CSVGObject_Line::Print(FILE *a) {

	int z;

	mfprintf(a,"  <line x1=\"%.6f\" y1=\"%.6f\" x2=\"%.6f\" y2=\"%.6f\"",
		m_fX1, m_fY1, m_fX2, m_fY2 );

	if (m_iLineCap != SVG_LINECAP_BUTT) {
		switch(m_iLineCap) {
			case SVG_LINECAP_BUTT:
				mfprintf(a," stroke-linecap=\"butt\"");
				break;
			case SVG_LINECAP_ROUND:
				mfprintf(a," stroke-linecap=\"round\"");
				break;
			case SVG_LINECAP_SQUARE:
				mfprintf(a," stroke-linecap=\"square\"");
				break;
		}
	}

	mfprintf(a," stroke=\"#%06lX\" stroke-width=\"%.6f\"",m_iStrokeColor,m_fStrokeWidth);
	if (m_fStrokeOpacity != 1)
		mfprintf(a," stroke-opacity=\"%.6f\"",m_fStrokeOpacity);
	if (m_bStrokeDash) {
		mfprintf(a," stroke-dasharray=\"");
		for (z=0;z<(int)m_faStrokeDashArray.size();z++) {
			mfprintf(a,"%.6f",m_faStrokeDashArray[z]);
			if (z+1 < (int)m_faStrokeDashArray.size())
				mfprintf(a,",");
		}
		mfprintf(a,"\" stroke-dashoffset=\"%.6f\"",m_fStrokeDashOffset);
	}

	if (m_iClipPath != -1)
		mfprintf(a," clip-path=\"url(#myclip%d)\"",m_iClipPath+1);

	if (m_sTitle.length() != 0)
		mfprintf(a,"><title>%s</title></line>\n",m_sTitle.c_str());
	else
		mfprintf(a," />\n");
}


void CSVGObject_Ellipse::Print(FILE *a) {

	int z;

	mfprintf(a,"  <ellipse cx=\"%.6f\" cy=\"%.6f\" rx=\"%.6f\" ry=\"%.6f\"",
		m_fCX, m_fCY, m_fRX, m_fRY );

	if (m_bFill) {
		mfprintf(a," fill=\"#%06lX\"",m_iFillColor);
		if (m_fFillOpacity != 1)
			mfprintf(a," fill-opacity=\"%.6f\"",m_fFillOpacity);
	} else
		mfprintf(a," fill=\"none\"");

	if (m_bStroke) {
		mfprintf(a," stroke=\"#%06lX\" stroke-width=\"%.6f\"",m_iStrokeColor,m_fStrokeWidth);
		if (m_fStrokeOpacity != 1)
			mfprintf(a," stroke-opacity=\"%.6f\"",m_fStrokeOpacity);
		if (m_bStrokeDash) {
			mfprintf(a," stroke-dasharray=\"");
			for (z=0;z<(int)m_faStrokeDashArray.size();z++) {
				mfprintf(a,"%.6f",m_faStrokeDashArray[z]);
				if (z+1 < (int)m_faStrokeDashArray.size())
					mfprintf(a,",");
			}
			mfprintf(a,"\" stroke-dashoffset=\"%.6f\"",m_fStrokeDashOffset);
		}
	} else
		mfprintf(a," stroke=\"none\"");

	if (m_iClipPath != -1)
		mfprintf(a," clip-path=\"url(#myclip%d)\"",m_iClipPath+1);

	if (m_sTitle.length() != 0)
		mfprintf(a,"><title>%s</title></ellipse>\n",m_sTitle.c_str());
	else
		mfprintf(a," />\n");
}


void CSVGObject_PolyLine::Print(FILE *a) {

	int z;

	mfprintf(a,"  <polyline points=\"");
	for (z=0;z<(int)m_faPoints.size();z++) {
		mfprintf(a,"%.6f",m_faPoints[z]);
		if ((z+1 < (int)m_faPoints.size()) && ((z % 2) == 0))
			mfprintf(a,", ");
		else if (z+1 < (int)m_faPoints.size())
			mfprintf(a," ");
	}
	mfprintf(a,"\"");

	if (m_iLineCap != SVG_LINECAP_BUTT) {
		switch(m_iLineCap) {
			case SVG_LINECAP_BUTT:
				mfprintf(a," stroke-linecap=\"butt\"");
				break;
			case SVG_LINECAP_ROUND:
				mfprintf(a," stroke-linecap=\"round\"");
				break;
			case SVG_LINECAP_SQUARE:
				mfprintf(a," stroke-linecap=\"square\"");
				break;
		}
	}

	if (m_iLineJoin != SVG_LINEJOIN_MITER) {
		switch(m_iLineJoin) {
			case SVG_LINEJOIN_MITER:
				mfprintf(a," stroke-linejoin=\"miter\"");
				break;
			case SVG_LINEJOIN_ROUND:
				mfprintf(a," stroke-linejoin=\"round\"");
				break;
			case SVG_LINEJOIN_BEVEL:
				mfprintf(a," stroke-linejoin=\"bevel\"");
				break;
		}
	}

	mfprintf(a," fill=\"none\" stroke=\"#%06lX\" stroke-width=\"%.6f\"",m_iStrokeColor,m_fStrokeWidth);
	if (m_fStrokeOpacity != 1)
		mfprintf(a," stroke-opacity=\"%.6f\"",m_fStrokeOpacity);
	if (m_bStrokeDash) {
		mfprintf(a," stroke-dasharray=\"");
		for (z=0;z<(int)m_faStrokeDashArray.size();z++) {
			mfprintf(a,"%.6f",m_faStrokeDashArray[z]);
			if (z+1 < (int)m_faStrokeDashArray.size())
				mfprintf(a,",");
		}
		mfprintf(a,"\" stroke-dashoffset=\"%.6f\"",m_fStrokeDashOffset);
	}

	if (m_iClipPath != -1)
		mfprintf(a," clip-path=\"url(#myclip%d)\"",m_iClipPath+1);

	if (m_sTitle.length() != 0)
		mfprintf(a,"><title>%s</title></polyline>\n",m_sTitle.c_str());
	else
		mfprintf(a," />\n");
}


void CSVGObject_Polygon::Print(FILE *a) {

	int z;

	mfprintf(a,"  <polygon points=\"");
	for (z=0;z<(int)m_faPoints.size();z++) {
		mfprintf(a,"%.6f",m_faPoints[z]);
		if ((z+1 < (int)m_faPoints.size()) && ((z % 2) == 0))
			mfprintf(a,",");
		else if (z+1 < (int)m_faPoints.size())
			mfprintf(a," ");
	}
	mfprintf(a,"\"");

	if (m_bFill) {
		mfprintf(a," fill=\"#%06lX\"",m_iFillColor);
		if (m_fFillOpacity != 1)
			mfprintf(a," fill-opacity=\"%.6f\"",m_fFillOpacity);
	} else
		mfprintf(a," fill=\"none\"");

	if (m_bStroke) {
		mfprintf(a," stroke=\"#%06lX\" stroke-width=\"%.6f\"",m_iStrokeColor,m_fStrokeWidth);
		if (m_fStrokeOpacity != 1)
			mfprintf(a," stroke-opacity=\"%.6f\"",m_fStrokeOpacity);
		if (m_bStrokeDash) {
			mfprintf(a," stroke-dasharray=\"");
			for (z=0;z<(int)m_faStrokeDashArray.size();z++) {
				mfprintf(a,"%.6f",m_faStrokeDashArray[z]);
				if (z+1 < (int)m_faStrokeDashArray.size())
					mfprintf(a,",");
			}
			mfprintf(a,"\" stroke-dashoffset=\"%.6f\"",m_fStrokeDashOffset);
		}
		if (m_iLineJoin != SVG_LINEJOIN_MITER) {
			switch(m_iLineJoin) {
				case SVG_LINEJOIN_MITER:
					mfprintf(a," stroke-linejoin=\"miter\"");
					break;
				case SVG_LINEJOIN_ROUND:
					mfprintf(a," stroke-linejoin=\"round\"");
					break;
				case SVG_LINEJOIN_BEVEL:
					mfprintf(a," stroke-linejoin=\"bevel\"");
					break;
			}
		}
	} else
		mfprintf(a," stroke=\"none\"");

	if (m_iClipPath != -1)
		mfprintf(a," clip-path=\"url(#myclip%d)\"",m_iClipPath+1);

	if (m_sTitle.length() != 0)
		mfprintf(a,"><title>%s</title></polygon>\n",m_sTitle.c_str());
	else
		mfprintf(a," />\n");
}


void CSVGObject_Text::Print(FILE *a) {

	int z;

	mfprintf(a,"  <text x=\"%.6f\" y=\"%.6f\" font-family=\"%s\" font-size=\"%.6f\"",
		m_fPosX, m_fPosY, m_sFontFamily.c_str(), m_fFontSize );

	if (m_iFontWeight != SVG_FONTWEIGHT_NORMAL) {
		switch(m_iFontWeight) {
			case SVG_FONTWEIGHT_NORMAL:
				mfprintf(a," font-weight=\"normal\"");
				break;
			case SVG_FONTWEIGHT_LIGHTER:
				mfprintf(a," font-weight=\"lighter\"");
				break;
			case SVG_FONTWEIGHT_BOLD:
				mfprintf(a," font-weight=\"bold\"");
				break;
			case SVG_FONTWEIGHT_BOLDER:
				mfprintf(a," font-weight=\"bolder\"");
				break;
			case SVG_FONTWEIGHT_100:
				mfprintf(a," font-weight=\"100\"");
				break;
			case SVG_FONTWEIGHT_200:
				mfprintf(a," font-weight=\"200\"");
				break;
			case SVG_FONTWEIGHT_300:
				mfprintf(a," font-weight=\"300\"");
				break;
			case SVG_FONTWEIGHT_400:
				mfprintf(a," font-weight=\"400\"");
				break;
			case SVG_FONTWEIGHT_500:
				mfprintf(a," font-weight=\"500\"");
				break;
			case SVG_FONTWEIGHT_600:
				mfprintf(a," font-weight=\"600\"");
				break;
			case SVG_FONTWEIGHT_700:
				mfprintf(a," font-weight=\"700\"");
				break;
			case SVG_FONTWEIGHT_800:
				mfprintf(a," font-weight=\"800\"");
				break;
			case SVG_FONTWEIGHT_900:
				mfprintf(a," font-weight=\"900\"");
				break;
		}
	}

	if (m_iFontStyle != SVG_FONTSTYLE_NORMAL) {
		switch(m_iFontStyle) {
			case SVG_FONTSTYLE_NORMAL:
				mfprintf(a," font-style=\"normal\"");
				break;
			case SVG_FONTSTYLE_ITALIC:
				mfprintf(a," font-style=\"italic\"");
				break;
		}
	}

	if (m_iFontStretch != SVG_FONTSTRETCH_NORMAL) {
		switch(m_iFontStretch) {
			case SVG_FONTSTRETCH_NORMAL:
				mfprintf(a," font-stretch=\"normal\"");
				break;
			case SVG_FONTSTRETCH_WIDER:
				mfprintf(a," font-stretch=\"wider\"");
				break;
			case SVG_FONTSTRETCH_NARROWER:
				mfprintf(a," font-stretch=\"narrower\"");
				break;
			case SVG_FONTSTRETCH_ULTRA_CONDENSED:
				mfprintf(a," font-stretch=\"ultra-condensed\"");
				break;
			case SVG_FONTSTRETCH_EXTRA_CONDENSED:
				mfprintf(a," font-stretch=\"extra-condensed\"");
				break;
			case SVG_FONTSTRETCH_CONDENSED:
				mfprintf(a," font-stretch=\"condensed\"");
				break;
			case SVG_FONTSTRETCH_SEMI_CONDENSED:
				mfprintf(a," font-stretch=\"semi-condensed\"");
				break;
			case SVG_FONTSTRETCH_SEMI_EXPANDED:
				mfprintf(a," font-stretch=\"semi-expanded\"");
				break;
			case SVG_FONTSTRETCH_EXPANDED:
				mfprintf(a," font-stretch=\"expanded\"");
				break;
			case SVG_FONTSTRETCH_EXTRA_EXPANDED:
				mfprintf(a," font-stretch=\"extra-expanded\"");
				break;
			case SVG_FONTSTRETCH_ULTRA_EXPANDED:
				mfprintf(a," font-stretch=\"ultra-expanded\"");
				break;
		}
	}

	if (m_iFontAlign != SVG_FONTALIGN_LEFT) {
		switch(m_iFontAlign) {
			case SVG_FONTALIGN_LEFT:
				mfprintf(a," text-anchor=\"begin\"");
				break;
			case SVG_FONTALIGN_CENTER:
				mfprintf(a," text-anchor=\"middle\"");
				break;
			case SVG_FONTALIGN_RIGHT:
				mfprintf(a," text-anchor=\"end\"");
				break;
		}
	}

	if (m_iFontBaseline != SVG_FONTBASELINE_BASELINE) {
		switch(m_iFontBaseline) {
			case SVG_FONTBASELINE_TOP:
				mfprintf(a," dominant-baseline=\"text-before-edge\"");
				break;
			case SVG_FONTBASELINE_CENTER:
				mfprintf(a," dominant-baseline=\"central\"");
				break;
			case SVG_FONTBASELINE_BOTTOM:
				mfprintf(a," dominant-baseline=\"text-after-edge\"");
				break;
		}
	}

	if (m_fAngle != 0)
		mfprintf(a," transform=\"rotate(%.6f %.6f %.6f)\"",m_fAngle,m_fPosX,m_fPosY);

	if (m_bFill) {
		mfprintf(a," fill=\"#%06lX\"",m_iFillColor);
		if (m_fFillOpacity != 1)
			mfprintf(a," fill-opacity=\"%.6f\"",m_fFillOpacity);
	} else
		mfprintf(a," fill=\"none\"");

	if (m_bStroke) {
		mfprintf(a," stroke=\"#%06lX\" stroke-width=\"%.6f\"",m_iStrokeColor,m_fStrokeWidth);
		if (m_fStrokeOpacity != 1)
			mfprintf(a," stroke-opacity=\"%.6f\"",m_fStrokeOpacity);
		if (m_bStrokeDash) {
			mfprintf(a," stroke-dasharray=\"");
			for (z=0;z<(int)m_faStrokeDashArray.size();z++) {
				mfprintf(a,"%.6f",m_faStrokeDashArray[z]);
				if (z+1 < (int)m_faStrokeDashArray.size())
					mfprintf(a,",");
			}
			mfprintf(a,"\" stroke-dashoffset=\"%.6f\"",m_fStrokeDashOffset);
		}
		if (m_iLineJoin != SVG_LINEJOIN_MITER) {
			switch(m_iLineJoin) {
				case SVG_LINEJOIN_MITER:
					mfprintf(a," stroke-linejoin=\"miter\"");
					break;
				case SVG_LINEJOIN_ROUND:
					mfprintf(a," stroke-linejoin=\"round\"");
					break;
				case SVG_LINEJOIN_BEVEL:
					mfprintf(a," stroke-linejoin=\"bevel\"");
					break;
			}
		}
	} else
		mfprintf(a," stroke=\"none\"");

	if (m_iClipPath != -1)
		mfprintf(a," clip-path=\"url(#myclip%d)\"",m_iClipPath+1);

	if (m_sTitle.length() != 0)
		mfprintf(a,"><title>%s</title>%s</text>\n",m_sTitle.c_str(),m_sText.c_str());
	else
		mfprintf(a,">%s</text>\n",m_sText.c_str());
}


void CSVGObject_Path::Print(FILE *a) {

	int z;

	mfprintf(a,"  <path d=\"");
	for (z=0;z<(int)m_oaPathElements.size();z++) {
		m_oaPathElements[z].Print(a);
		if (z+1 < (int)m_oaPathElements.size())
			mfprintf(a," ");
	}
	mfprintf(a,"\"");

	if (m_bFill) {
		mfprintf(a," fill=\"#%06lX\"",m_iFillColor);
		if (m_fFillOpacity != 1)
			mfprintf(a," fill-opacity=\"%.6f\"",m_fFillOpacity);
	} else
		mfprintf(a," fill=\"none\"");

	if (m_bStroke) {
		mfprintf(a," stroke=\"#%06lX\" stroke-width=\"%.6f\"",m_iStrokeColor,m_fStrokeWidth);
		if (m_fStrokeOpacity != 1)
			mfprintf(a," stroke-opacity=\"%.6f\"",m_fStrokeOpacity);
		if (m_bStrokeDash) {
			mfprintf(a," stroke-dasharray=\"");
			for (z=0;z<(int)m_faStrokeDashArray.size();z++) {
				mfprintf(a,"%.6f",m_faStrokeDashArray[z]);
				if (z+1 < (int)m_faStrokeDashArray.size())
					mfprintf(a,",");
			}
			mfprintf(a,"\" stroke-dashoffset=\"%.6f\"",m_fStrokeDashOffset);
		}
		if (m_iLineJoin != SVG_LINEJOIN_MITER) {
			switch(m_iLineJoin) {
				case SVG_LINEJOIN_MITER:
					mfprintf(a," stroke-linejoin=\"miter\"");
					break;
				case SVG_LINEJOIN_ROUND:
					mfprintf(a," stroke-linejoin=\"round\"");
					break;
				case SVG_LINEJOIN_BEVEL:
					mfprintf(a," stroke-linejoin=\"bevel\"");
					break;
			}
		}
	} else
		mfprintf(a," stroke=\"none\"");

	if (m_iClipPath != -1)
		mfprintf(a," clip-path=\"url(#myclip%d)\"",m_iClipPath+1);

	if (m_sTitle.length() != 0)
		mfprintf(a,"><title>%s</title></polygon>\n",m_sTitle.c_str());
	else
		mfprintf(a," />\n");
}



void CSVGPathElement::Print(FILE *a) {
	switch(m_iType) {
		case SVG_PATHELEM_MOVETO:
			mfprintf( a, "M %.6f %.6f", m_faValues[0], m_faValues[1] );
			break;
		case SVG_PATHELEM_MOVEREL:
			mfprintf( a, "m %.6f %.6f", m_faValues[0], m_faValues[1] );
			break;
		case SVG_PATHELEM_LINETO:
			mfprintf( a, "L %.6f %.6f", m_faValues[0], m_faValues[1] );
			break;
		case SVG_PATHELEM_LINEREL:
			mfprintf( a, "l %.6f %.6f", m_faValues[0], m_faValues[1] );
			break;
		case SVG_PATHELEM_HORILINETO:
			mfprintf( a, "H %.6f", m_faValues[0] );
			break;
		case SVG_PATHELEM_HORILINEREL:
			mfprintf( a, "h %.6f", m_faValues[0] );
			break;
		case SVG_PATHELEM_VERTLINETO:
			mfprintf( a, "V %.6f", m_faValues[0] );
			break;
		case SVG_PATHELEM_VERTLINEREL:
			mfprintf( a, "v %.6f", m_faValues[0] );
			break;
		case SVG_PATHELEM_CLOSE:
			mfprintf( a, "Z" );
			break;
		case SVG_PATHELEM_CUBICBEZIERTO:
			mfprintf( a, "C %.6f %.6f %.6f %.6f %.6f %.6f", m_faValues[0], m_faValues[1], m_faValues[2], m_faValues[3], m_faValues[4], m_faValues[5] );
			break;
		case SVG_PATHELEM_CUBICBEZIERREL:
			mfprintf( a, "c %.6f %.6f %.6f %.6f %.6f %.6f", m_faValues[0], m_faValues[1], m_faValues[2], m_faValues[3], m_faValues[4], m_faValues[5] );
			break;
		case SVG_PATHELEM_CUBICBEZIERCONTINUETO:
			mfprintf( a, "S %.6f %.6f %.6f %.6f", m_faValues[0], m_faValues[1], m_faValues[2], m_faValues[3] );
			break;
		case SVG_PATHELEM_CUBICBEZIERCONTINUEREL:
			mfprintf( a, "s %.6f %.6f %.6f %.6f", m_faValues[0], m_faValues[1], m_faValues[2], m_faValues[3] );
			break;
		case SVG_PATHELEM_QUADRATICBEZIERTO:
			mfprintf( a, "Q %.6f %.6f %.6f %.6f", m_faValues[0], m_faValues[1], m_faValues[2], m_faValues[3] );
			break;
		case SVG_PATHELEM_QUADRATICBEZIERREL:
			mfprintf( a, "q %.6f %.6f %.6f %.6f", m_faValues[0], m_faValues[1], m_faValues[2], m_faValues[3] );
			break;
		case SVG_PATHELEM_QUADRATICBEZIERCONTINUETO:
			mfprintf( a, "T %.6f %.6f", m_faValues[0], m_faValues[1] );
			break;
		case SVG_PATHELEM_QUADRATICBEZIERCONTINUEREL:
			mfprintf( a, "t %.6f %.6f", m_faValues[0], m_faValues[1] );
			break;
		case SVG_PATHELEM_ARCTO:
			mfprintf( a, "A %.6f %.6f %.6f %d %d %.6f %.6f", m_faValues[0], m_faValues[1], m_faValues[2], m_iaValues[0], m_iaValues[1], m_faValues[3], m_faValues[4] );
			break;
		case SVG_PATHELEM_ARCREL:
			mfprintf( a, "a %.6f %.6f %.6f %d %d %.6f %.6f", m_faValues[0], m_faValues[1], m_faValues[2], m_iaValues[0], m_iaValues[1], m_faValues[3], m_faValues[4] );
			break;
	}
}



void CSVGWriter::SetClipPolygon(int n, ...) {

	CSVGObject_ClipPolygon *o;
	int z;
	va_list vl;

	o = new CSVGObject_ClipPolygon();

	va_start(vl,n);
	for (z=0;z<n;z++)
		o->m_faPoints.push_back( va_arg(vl,double) );
	va_end(vl);

	m_oaClipPolygons.push_back(o);

	m_iaClipPolygon.resize(1);
	m_iaClipPolygon[0] = (int)m_oaClipPolygons.size() - 1;
}



void CSVGWriter::ResetClipPolygon() {

	m_iaClipPolygon.resize(1);
	m_iaClipPolygon[0] = -1;
}



void CSVGWriter::PushClipPolygon(int n, ...) {

	CSVGObject_ClipPolygon *o;
	int z;
	va_list vl;

	o = new CSVGObject_ClipPolygon();

	va_start(vl,n);
	for (z=0;z<n;z++)
		o->m_faPoints.push_back( va_arg(vl,double) );
	va_end(vl);

	m_oaClipPolygons.push_back(o);

	m_iaClipPolygon.push_back((int)m_oaClipPolygons.size() - 1);
}



void CSVGWriter::PopClipPolygon() {

	if (m_iaClipPolygon.size() > 1)
		m_iaClipPolygon.pop_back();
}







