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


#ifndef GRACE_H
#define GRACE_H


// This must always be the first include directive
#include "config.h"

#include "xobject.h"
#include "xobarray.h"
#include "xdoublearray.h"


class CGraceDataset;


class CGraceGraph : public CxObject
{
public:
	CGraceGraph(); 
	~CGraceGraph(); 

	void WriteData(FILE *a, bool silent);
	void Write(FILE *a);
	CGraceDataset* CurrentDataset();
	CGraceDataset* Dataset(int i);

	CxObArray m_oaDatasets;
	CxObArray m_oaCustomLabelsX;
	CxObArray m_oaCustomLabelsY;
	CxObArray m_oaLines;

	double m_fYAxisBarWidth;
	bool m_bShowFrame;
	bool m_bShowXAxis;

	double m_fMinRangeX;
	double m_fMaxRangeX;
	double m_fMinRangeY;
	double m_fMaxRangeY;
	double m_fTickMajorX;
	int m_iTickPrecX;
	double m_fTickMajorY;
	int m_iTickPrecY;
	int m_iTickMinorX;
	int m_iTickMinorY;

	double m_fMinValX;
	double m_fMaxValX;
	double m_fMinValY;
	double m_fMaxValY;

	char *m_sTitle;
	char *m_sLabelX;
	char *m_sLabelY;
	char *m_sSubTitle;

	bool m_bTicksBothSidesX;
	bool m_bTicksBothSidesY;
	bool m_bTickLabelsBothSidesX;
	bool m_bTickLabelsBothSidesY;

	bool m_bInvertXAxis;
	bool m_bInvertYAxis;
	bool m_bInvert;

	int m_iColorIndex;
	int m_iNumber;

	double m_fViewportX1;
	double m_fViewportY1;
	double m_fViewportX2;
	double m_fViewportY2;

	bool m_bTicks;
	bool m_bTickLabels;
	bool m_bTickInX;
	bool m_bTickInY;
	bool m_bLegend;
	double m_fFrameWidth;
};


class CGraceColor : public CxObject
{
public:
	void SetName(const char *s);
	CGraceColor();
	~CGraceColor();

	unsigned char m_iColorRed;
	unsigned char m_iColorGreen;
	unsigned char m_iColorBlue;
	char *m_sName;
};


class CGraceDataset : public CxObject
{
public:
	void CopyFrom(CGraceDataset *d);
	int m_iNumber;
	void WriteSet(FILE *a);
	void WriteHeader(FILE *a);
	CGraceDataset();
	~CGraceDataset();

	CxDoubleArray m_faValues;
//	unsigned long m_iLineColor;
//	unsigned long m_iSymbColor;
	int m_iLineColorIndex;
	int m_iSymbColorIndex;
	char *m_sName;
	double m_fLineWidth;
	unsigned char m_iLineStyle;
	int m_iBegin;
	int m_iEnd;
	bool m_bFill;
	bool m_bSymbol;
	CGraceGraph *m_pGraph;
};


class CGraceLine : public CxObject
{
public:
	CGraceLine() { }
	~CGraceLine() { }
	void WriteLine(FILE *a, int graph);

	double m_fX1;
	double m_fY1;
	double m_fX2;
	double m_fY2;
	double m_fLineWidth;
	int m_iLineStyle;
	int m_iLineColorIndex;
};


class CGraceCustomLabel : public CxObject
{
public:
	void Write(FILE *a);
	CGraceCustomLabel() { }
	~CGraceCustomLabel() { }

	bool m_bX;
	int m_iNumber;
	bool m_bMajor;
	char *m_sText;
	double m_fValue;
};


class CGrace : public CxObject  
{
public:
	int AddColor(unsigned char r, unsigned char g, unsigned char b, const char *name);
	void WriteCSV(const char *s);
	void SetDatasetName(const char *s);
	void SetDatasetName(int set, const char *s);
	CGraceDataset* LastDataset();
	void SetViewport(double x1, double y1, double x2, double y2);
	void AddGraph();
	CGraceGraph* CurrentGraph();
	void AddLine(double x1, double y1, double x2, double y2, double width, int style);
	void AddLine(double x1, double y1, double x2, double y2, double width, int style, unsigned char r, unsigned char g, unsigned char b);
	void AddCustomLabelY(bool major, double val, const char *s);
	void AddCustomLabelX(bool major, double val, const char *s);
	void SetSetLineWidth(int set, double width);
	void SetSetLineColorLong(int set, unsigned long col);
	void SetSetLineColor(int set, unsigned char r, unsigned char g, unsigned char b);
	void SetSetLineColor(unsigned char r, unsigned char g, unsigned char b);
	void SetSetRange(int set, int start, int end);
	void DuplicateSet(int set);
	void SetRangeY(double mi, double ma);
	void SetRangeX(double mi, double ma);
	void SetLabelY(const char *s);
	void SetLabelX(const char *s);
	void SetSubTitle(const char *s);
	void SetTitle(const char *s);
	void FindMinMaxVal();
	void AddXYTupel(int set, double x, double y);
	void AddXYTupel(double x, double y);
	void AddDataset();
	void WriteAgr(const char *s, bool silent);

	CxObArray m_oaGraceGraphs;
	CxObArray m_oaGraceColors;

	void MakeTicks();

	unsigned long m_iBGColor;

	CGrace();
	~CGrace();

};

#endif 
