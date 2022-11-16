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


#ifndef SVGWRITER_H
#define SVGWRITER_H


// This must always be the first include directive
#include "config.h"

#include <stdio.h>
#include <vector>
#include <string>
#include "tools.h"


#define SVG_LINECAP_BUTT    0
#define SVG_LINECAP_SQUARE  1
#define SVG_LINECAP_ROUND   2

#define SVG_LINEJOIN_MITER  0
#define SVG_LINEJOIN_ROUND  1
#define SVG_LINEJOIN_BEVEL  2

#define SVG_FONTWEIGHT_NORMAL    0
#define SVG_FONTWEIGHT_LIGHTER   1
#define SVG_FONTWEIGHT_BOLD      2
#define SVG_FONTWEIGHT_BOLDER    3
#define SVG_FONTWEIGHT_100       4
#define SVG_FONTWEIGHT_200       5
#define SVG_FONTWEIGHT_300       6
#define SVG_FONTWEIGHT_400       7
#define SVG_FONTWEIGHT_500       8
#define SVG_FONTWEIGHT_600       9
#define SVG_FONTWEIGHT_700      10
#define SVG_FONTWEIGHT_800      11
#define SVG_FONTWEIGHT_900      12

#define SVG_FONTSTYLE_NORMAL  0
#define SVG_FONTSTYLE_ITALIC  1

#define SVG_FONTSTRETCH_NORMAL            0
#define SVG_FONTSTRETCH_WIDER             1
#define SVG_FONTSTRETCH_NARROWER          2
#define SVG_FONTSTRETCH_ULTRA_CONDENSED   3
#define SVG_FONTSTRETCH_EXTRA_CONDENSED   4
#define SVG_FONTSTRETCH_CONDENSED         5
#define SVG_FONTSTRETCH_SEMI_CONDENSED    6
#define SVG_FONTSTRETCH_SEMI_EXPANDED     7
#define SVG_FONTSTRETCH_EXPANDED          8
#define SVG_FONTSTRETCH_EXTRA_EXPANDED    9
#define SVG_FONTSTRETCH_ULTRA_EXPANDED   10

#define SVG_FONTALIGN_LEFT    0
#define SVG_FONTALIGN_CENTER  1
#define SVG_FONTALIGN_RIGHT   2

#define SVG_FONTBASELINE_TOP       0
#define SVG_FONTBASELINE_CENTER    1
#define SVG_FONTBASELINE_BASELINE  2
#define SVG_FONTBASELINE_BOTTOM    3

#define SVG_PATHELEM_MOVETO                       0
#define SVG_PATHELEM_MOVEREL                      1
#define SVG_PATHELEM_LINETO                       2
#define SVG_PATHELEM_LINEREL                      3
#define SVG_PATHELEM_HORILINETO                   4
#define SVG_PATHELEM_HORILINEREL                  5
#define SVG_PATHELEM_VERTLINETO                   6
#define SVG_PATHELEM_VERTLINEREL                  7
#define SVG_PATHELEM_CLOSE                        8
#define SVG_PATHELEM_CUBICBEZIERTO                9
#define SVG_PATHELEM_CUBICBEZIERREL              10
#define SVG_PATHELEM_CUBICBEZIERCONTINUETO       11
#define SVG_PATHELEM_CUBICBEZIERCONTINUEREL      12
#define SVG_PATHELEM_QUADRATICBEZIERTO           13
#define SVG_PATHELEM_QUADRATICBEZIERREL          14
#define SVG_PATHELEM_QUADRATICBEZIERCONTINUETO   15
#define SVG_PATHELEM_QUADRATICBEZIERCONTINUEREL  16
#define SVG_PATHELEM_ARCTO                       17
#define SVG_PATHELEM_ARCREL                      18



class CSVGObject {
public:

	virtual void Print(FILE *a) = 0;
};



class CSVGObject_Rectangle : public CSVGObject {
public:

	void Print(FILE *a);

	double m_fX1;
	double m_fY1;
	double m_fX2;
	double m_fY2;
	double m_fRX;
	double m_fRY;

	bool m_bFill;
	unsigned long m_iFillColor;
	double m_fFillOpacity;

	bool m_bStroke;
	double m_fStrokeWidth;
	unsigned long m_iStrokeColor;
	double m_fStrokeOpacity;
	bool m_bStrokeDash;
	std::vector<double> m_faStrokeDashArray;
	double m_fStrokeDashOffset;
	char m_iLineJoin;
	int m_iClipPath;

	std::string m_sTitle;
};



class CSVGObject_Line : public CSVGObject {
public:

	void Print(FILE *a);

	double m_fX1;
	double m_fY1;
	double m_fX2;
	double m_fY2;

	double m_fStrokeWidth;
	unsigned long m_iStrokeColor;
	double m_fStrokeOpacity;
	bool m_bStrokeDash;
	std::vector<double> m_faStrokeDashArray;
	double m_fStrokeDashOffset;
	char m_iLineCap;
	int m_iClipPath;

	std::string m_sTitle;
};



class CSVGObject_Ellipse : public CSVGObject {
public:

	void Print(FILE *a);

	double m_fCX;
	double m_fCY;
	double m_fRX;
	double m_fRY;

	bool m_bFill;
	unsigned long m_iFillColor;
	double m_fFillOpacity;

	bool m_bStroke;
	double m_fStrokeWidth;
	unsigned long m_iStrokeColor;
	double m_fStrokeOpacity;
	bool m_bStrokeDash;
	std::vector<double> m_faStrokeDashArray;
	double m_fStrokeDashOffset;
	int m_iClipPath;

	std::string m_sTitle;
};



class CSVGObject_PolyLine : public CSVGObject {
public:

	void Print(FILE *a);

	std::vector<double> m_faPoints;

	double m_fStrokeWidth;
	unsigned long m_iStrokeColor;
	double m_fStrokeOpacity;
	bool m_bStrokeDash;
	std::vector<double> m_faStrokeDashArray;
	double m_fStrokeDashOffset;
	char m_iLineCap;
	char m_iLineJoin;
	int m_iClipPath;

	std::string m_sTitle;
};



class CSVGObject_Polygon : public CSVGObject {
public:

	void Print(FILE *a);

	std::vector<double> m_faPoints;

	bool m_bFill;
	unsigned long m_iFillColor;
	double m_fFillOpacity;

	bool m_bStroke;
	double m_fStrokeWidth;
	unsigned long m_iStrokeColor;
	double m_fStrokeOpacity;
	bool m_bStrokeDash;
	std::vector<double> m_faStrokeDashArray;
	double m_fStrokeDashOffset;
	char m_iLineJoin;
	int m_iClipPath;

	std::string m_sTitle;
};



class CSVGObject_ClipPolygon : public CSVGObject {
public:

	void Print(FILE *a) {
		UNUSED(a);
	}

	std::vector<double> m_faPoints;

	std::string m_sTitle;
};



class CSVGPathElement {
public:

	CSVGPathElement( int type, ... ) {
		va_list vl;
		m_iType = type;
		va_start(vl,type);
		switch(m_iType) {
			case SVG_PATHELEM_MOVETO:
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_MOVEREL:
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_LINETO:
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_LINEREL:
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_HORILINETO:
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_HORILINEREL:
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_VERTLINETO:
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_VERTLINEREL:
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_CLOSE:
				break;
			case SVG_PATHELEM_CUBICBEZIERTO:
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_CUBICBEZIERREL:
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_CUBICBEZIERCONTINUETO:
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_CUBICBEZIERCONTINUEREL:
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_QUADRATICBEZIERTO:
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_QUADRATICBEZIERREL:
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_QUADRATICBEZIERCONTINUETO:
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_QUADRATICBEZIERCONTINUEREL:
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_ARCTO:
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_iaValues.push_back( va_arg( vl, int ) );
				m_iaValues.push_back( va_arg( vl, int ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				break;
			case SVG_PATHELEM_ARCREL:
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_iaValues.push_back( va_arg( vl, int ) );
				m_iaValues.push_back( va_arg( vl, int ) );
				m_faValues.push_back( va_arg( vl, double ) );
				m_faValues.push_back( va_arg( vl, double ) );
				break;
		}
		va_end(vl);
	}


	void Print(FILE *a);

	int m_iType;
	std::vector<double> m_faValues;
	std::vector<int> m_iaValues;
};



class CSVGObject_Path : public CSVGObject {
public:

	void Print(FILE *a);

	std::vector<CSVGPathElement> m_oaPathElements;

	void AddMoveTo(double x, double y) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_MOVETO, x, y ) );
	}

	void AddMoveRel(double dx, double dy) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_MOVEREL, dx, dy ) );
	}

	void AddLineTo(double x, double y) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_LINETO, x, y ) );
	}

	void AddLineRel(double dx, double dy) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_LINEREL, dx, dy ) );
	}

	void AddHorizontalLineTo(double x) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_HORILINETO, x ) );
	}

	void AddHorizontalLineRel(double dx) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_HORILINEREL, dx ) );
	}

	void AddVerticalLineTo(double y) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_VERTLINETO, y ) );
	}

	void AddVerticalLineRel(double dy) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_VERTLINEREL, dy ) );
	}

	void AddClose() {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_CLOSE ) );
	}

	void AddCubicBezierTo(double x1, double y1, double x2, double y2, double x, double y) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_CUBICBEZIERTO, x1, y1, x2, y2, x, y ) );
	}

	void AddCubicBezierRel(double dx1, double dy1, double dx2, double dy2, double dx, double dy) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_CUBICBEZIERREL, dx1, dy1, dx2, dy2, dx, dy ) );
	}

	void AddCubicBezierContinueTo(double x2, double y2, double x, double y) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_CUBICBEZIERCONTINUETO, x2, y2, x, y ) );
	}

	void AddCubicBezierContinueRel(double dx2, double dy2, double dx, double dy) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_CUBICBEZIERCONTINUEREL, dx2, dy2, dx, dy ) );
	}

	void AddQuadraticBezierTo(double x1, double y1, double x, double y) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_QUADRATICBEZIERTO, x1, y1, x, y ) );
	}

	void AddQuadraticBezierRel(double dx1, double dy1, double dx, double dy) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_QUADRATICBEZIERREL, dx1, dy1, dx, dy ) );
	}

	void AddQuadraticBezierContinueTo(double x, double y) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_QUADRATICBEZIERCONTINUETO, x, y ) );
	}

	void AddQuadraticBezierContinueRel(double dx, double dy) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_QUADRATICBEZIERCONTINUEREL, dx, dy ) );
	}

	void AddArcTo(double rx, double ry, double rotx, int flong, int fclock, double x, double y) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_ARCTO, rx, ry, rotx, flong, fclock, x, y ) );
	}

	void AddArcRel(double rx, double ry, double rotx, int flong, int fclock, double dx, double dy) {
		m_oaPathElements.push_back( CSVGPathElement( SVG_PATHELEM_ARCREL, rx, ry, rotx, flong, fclock, dx, dy ) );
	}


	bool m_bFill;
	unsigned long m_iFillColor;
	double m_fFillOpacity;

	bool m_bStroke;
	double m_fStrokeWidth;
	unsigned long m_iStrokeColor;
	double m_fStrokeOpacity;
	bool m_bStrokeDash;
	std::vector<double> m_faStrokeDashArray;
	double m_fStrokeDashOffset;
	char m_iLineJoin;
	int m_iClipPath;

	std::string m_sTitle;
};



class CSVGObject_Text : public CSVGObject {
public:

	void Print(FILE *a);

	std::string m_sText;

	double m_fPosX;
	double m_fPosY;
	double m_fAngle;
	double m_fFontSize;

	bool m_bFill;
	unsigned long m_iFillColor;
	double m_fFillOpacity;

	bool m_bStroke;
	double m_fStrokeWidth;
	unsigned long m_iStrokeColor;
	double m_fStrokeOpacity;
	bool m_bStrokeDash;
	std::vector<double> m_faStrokeDashArray;
	double m_fStrokeDashOffset;
	char m_iLineJoin;

	std::string m_sFontFamily;
	char m_iFontWeight;
	char m_iFontStyle;
	char m_iFontStretch;
	char m_iFontAlign;
	char m_iFontBaseline;
	int m_iClipPath;

	std::string m_sTitle;
};



class CSVGWriter {
public:

	CSVGWriter();

	virtual ~CSVGWriter();

	void SetUserSize(double width, double height);

	void SetPixelSize(int width);

	bool WriteSVG(const char *s);

	CSVGObject_Path* AddPath();

	void AddRectangle(double x1, double y1, double x2, double y2);
	void AddRectangle(double x1, double y1, double x2, double y2, const char *title);
	void AddRectangle(double x1, double y1, double x2, double y2, double radius);
	void AddRectangle(double x1, double y1, double x2, double y2, double radius, const char *title);
	void AddRectangle(double x1, double y1, double x2, double y2, double radx, double rady);
	void AddRectangle(double x1, double y1, double x2, double y2, double radx, double rady, const char *title);

	void AddLine(double x1, double y1, double x2, double y2);
	void AddLine(double x1, double y1, double x2, double y2, const char *title);

	void AddCircle(double cx, double cy, double r);
	void AddCircle(double cx, double cy, double r, const char *title);

	void AddEllipse(double cx, double cy, double rx, double ry);
	void AddEllipse(double cx, double cy, double rx, double ry, const char *title);

	void AddPolyLine_Variadic(int n, ...);
	void AddPolyLine_Variadic(const char *title, int n, ...);
	void AddPolyLine_Array(int n, double *pts);

	void AddPolygon_Variadic(int n, ...);
	void AddPolygon_Variadic(const char *title, int n, ...);
	void AddPolygon_Array(int n, double *pts);
	void AddPolygon_Array(const std::vector<double> *vec);

	void AddText(double x, double y, const char *text);
	void AddText(double x, double y, const char *text, const char *title);
	void AddText(double x, double y, const char *text, double angle);
	void AddText(double x, double y, const char *text, double angle, const char *title);

	/********************************************************************************************/

	void SetClipPolygon(int n, ...);
	void ResetClipPolygon();
	void PushClipPolygon(int n, ...);
	void PopClipPolygon();

	void SetStroke(bool b) {
		m_iaStroke.resize(1);
		m_iaStroke[0] = b?1:0;
	}

	void PushStroke(bool b) {
		m_iaStroke.push_back(b?1:0);
	}

	void PopStroke() {
		if (m_iaStroke.size() > 1)
			m_iaStroke.pop_back();
	}

	void SetStrokeWidth(double w) {
		m_faStrokeWidth.resize(1);
		m_faStrokeWidth[0] = w;
	}

	void PushStrokeWidth(double w) {
		m_faStrokeWidth.push_back(w);
	}

	void PopStrokeWidth() {
		if (m_faStrokeWidth.size() > 1)
			m_faStrokeWidth.pop_back();
	}

	void SetStrokeColor(long c) {
		m_iaStrokeColor.resize(1);
		m_iaStrokeColor[0] = c;
	}

	void PushStrokeColor(long c) {
		m_iaStrokeColor.push_back(c);
	}

	void PopStrokeColor() {
		if (m_iaStrokeColor.size() > 1)
			m_iaStrokeColor.pop_back();
	}

	void SetStrokeOpacity(double o) {
		m_faStrokeOpacity.resize(1);
		m_faStrokeOpacity[0] = o;
	}

	void PushStrokeOpacity(double o) {
		m_faStrokeOpacity.push_back(o);
	}

	void PopStrokeOpacity() {
		if (m_faStrokeOpacity.size() > 1)
			m_faStrokeOpacity.pop_back();
	}

	void SetStrokeDash(bool b) {
		m_iaStrokeDash.resize(1);
		m_iaStrokeDash[0] = b?1:0;
	}

	void PushStrokeDash(bool b) {
		m_iaStrokeDash.push_back(b?1:0);
	}

	void PopStrokeDash() {
		if (m_iaStrokeDash.size() > 1)
			m_iaStrokeDash.pop_back();
	}

	void SetStrokeDashArray(int n, ...) {
		va_list vl;
		m_faaStrokeDashArray.resize(1);
		m_faaStrokeDashArray[0].clear();
		va_start(vl,n);
		for (int z=0;z<n;z++)
			m_faaStrokeDashArray[0].push_back( va_arg(vl,double) );
		va_end(vl);
	}

	void PushStrokeDashArray(int n, ...) {
		va_list vl;
		m_faaStrokeDashArray.push_back(std::vector<double>());
		va_start(vl,n);
		for (int z=0;z<n;z++)
			m_faaStrokeDashArray.back().push_back( va_arg(vl,double) );
		va_end(vl);
	}

	void PopStrokeDashArray() {
		m_faaStrokeDashArray.pop_back();
	}

	void SetStrokeDashOffset(double o) {
		m_faStrokeDashOffset.resize(1);
		m_faStrokeDashOffset[0] = o;
	}

	void PushStrokeDashOffset(double o) {
		m_faStrokeDashOffset.push_back(o);
	}

	void PopStrokeDashOffset() {
		if (m_faStrokeDashOffset.size() > 1)
			m_faStrokeDashOffset.pop_back();
	}

	void SetFill(bool b) {
		m_iaFill.resize(1);
		m_iaFill[0] = b?1:0;
	}

	void PushFill(bool b) {
		m_iaFill.push_back(b?1:0);
	}

	void PopFill() {
		if (m_iaFill.size() > 1)
			m_iaFill.pop_back();
	}

	void SetFillColor(long c) {
		m_iaFillColor.resize(1);
		m_iaFillColor[0] = c;
	}

	void PushFillColor(long c) {
		m_iaFillColor.push_back(c);
	}

	void PopFillColor() {
		if (m_iaFillColor.size() > 1)
			m_iaFillColor.pop_back();
	}

	void SetFillOpacity(double o) {
		m_faFillOpacity.resize(1);
		m_faFillOpacity[0] = o;
	}

	void PushFillOpacity(double o) {
		m_faFillOpacity.push_back(o);
	}

	void PopFillOpacity() {
		if (m_faFillOpacity.size() > 1)
			m_faFillOpacity.pop_back();
	}

	void SetLineCap(char i) {
		m_iaLineCap.resize(1);
		m_iaLineCap[0] = i;
	}

	void PushLineCap(char i) {
		m_iaLineCap.push_back(i);
	}

	void PopLineCap() {
		if (m_iaLineCap.size() > 1)
			m_iaLineCap.pop_back();
	}

	void SetLineJoin(char i) {
		m_iaLineJoin.resize(1);
		m_iaLineJoin[0] = i;
	}

	void PushLineJoin(char i) {
		m_iaLineJoin.push_back(i);
	}

	void PopLineJoin() {
		if (m_iaLineJoin.size() > 1)
			m_iaLineJoin.pop_back();
	}

	void SetFontFamily(const char *p) {
		m_saFontFamily.resize(1);
		m_saFontFamily[0] = p;
	}

	void PushFontFamily(const char *p) {
		m_saFontFamily.push_back(p);
	}

	void PopFontFamily() {
		if (m_saFontFamily.size() > 1)
			m_saFontFamily.pop_back();
	}

	void SetFontSize(double o) {
		m_faFontSize.resize(1);
		m_faFontSize[0] = o;
	}

	void PushFontSize(double o) {
		m_faFontSize.push_back(o);
	}

	void PopFontSize() {
		if (m_faFontSize.size() > 1)
			m_faFontSize.pop_back();
	}

	void SetFontWeight(char i) {
		m_iaFontWeight.resize(1);
		m_iaFontWeight[0] = i;
	}

	void PushFontWeight(char i) {
		m_iaFontWeight.push_back(i);
	}

	void PopFontWeight() {
		if (m_iaFontWeight.size() > 1)
			m_iaFontWeight.pop_back();
	}

	void SetFontStyle(char i) {
		m_iaFontStyle.resize(1);
		m_iaFontStyle[0] = i;
	}

	void PushFontStyle(char i) {
		m_iaFontStyle.push_back(i);
	}

	void PopFontStyle() {
		if (m_iaFontStyle.size() > 1)
			m_iaFontStyle.pop_back();
	}

	void SetFontStretch(char i) {
		m_iaFontStretch.resize(1);
		m_iaFontStretch[0] = i;
	}

	void PushFontStretch(char i) {
		m_iaFontStretch.push_back(i);
	}

	void PopFontStretch() {
		if (m_iaFontStretch.size() > 1)
			m_iaFontStretch.pop_back();
	}

	void SetFontAlign(char i) {
		m_iaFontAlign.resize(1);
		m_iaFontAlign[0] = i;
	}

	void PushFontAlign(char i) {
		m_iaFontAlign.push_back(i);
	}

	void PopFontAlign() {
		if (m_iaFontAlign.size() > 1)
			m_iaFontAlign.pop_back();
	}

	void SetFontBaseline(char i) {
		m_iaFontBaseline.resize(1);
		m_iaFontBaseline[0] = i;
	}

	void PushFontBaseline(char i) {
		m_iaFontBaseline.push_back(i);
	}

	void PopFontBaseline() {
		if (m_iaFontBaseline.size() > 1)
			m_iaFontBaseline.pop_back();
	}

	/********************************************************************************************/

	std::vector<CSVGObject*> m_oaObjects;

	double m_fUserWidth;
	double m_fUserHeight;

	int m_iPixelWidth;
	int m_iPixelHeight;

	/********************************************************************************************/

	std::vector<char> m_iaStroke;
	std::vector<double> m_faStrokeWidth;
	std::vector<unsigned long> m_iaStrokeColor;
	std::vector<double> m_faStrokeOpacity;
	std::vector<char> m_iaStrokeDash;
	std::vector<std::vector<double> > m_faaStrokeDashArray;
	std::vector<double> m_faStrokeDashOffset;

	std::vector<char> m_iaFill;
	std::vector<unsigned long> m_iaFillColor;
	std::vector<double> m_faFillOpacity;

	std::vector<char> m_iaLineCap;
	std::vector<char> m_iaLineJoin;

	std::vector<std::string> m_saFontFamily;
	std::vector<double> m_faFontSize;
	std::vector<char> m_iaFontWeight;
	std::vector<char> m_iaFontStyle;
	std::vector<char> m_iaFontStretch;
	std::vector<char> m_iaFontAlign;
	std::vector<char> m_iaFontBaseline;
	std::vector<int> m_iaClipPolygon;

	std::vector<CSVGObject_ClipPolygon*> m_oaClipPolygons;

};



#endif




