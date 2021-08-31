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


#ifndef MATRIXPLOT_H
#define MATRIXPLOT_H


// This must always be the first include directive
#include "config.h"

#include <vector>
#include <string>



class CMatrixPlot {
public:

	CMatrixPlot();

	~CMatrixPlot();

	void Init(int rows, int cols);

	void SetRowLabel(int row, const char *label);

	void SetColumnLabel(int column, const char *label);

	void SetLabelShow(bool top, bool bottom, bool left, bool right);

	void AddRowCategory(int i);

	void AddColumnCategory(int i);

	void AddRowSubCategory(int i);

	void AddColumnSubCategory(int i);

	void WriteSVGPlot(const char *s);

	unsigned long ColorFunction(double v);

	unsigned long ColorFunction2D(double x, double y);

	void SetPlotExponent(double e);

	void SetPlotExponent2(double e);

	void SetRangeZero(bool b);

	void SetRangeZero2(bool b);

	void SetValueLabel(const char *p);

	void SetValueLabel2(const char *p);

	void SetInvertColorScale(bool b);

	void SetInvertColorScale2(bool b);

	void Set2DMode(bool b);

	void SetRangeCut(double d);

	void SetRangeCut2(double d);

	void SetRange(double mi, double ma);

	void ResetRangeCut();

	void ResetRangeCut2();

	int m_iRows;
	int m_iCols;

	std::vector<std::string> m_saRowLabels;
	std::vector<std::string> m_saColumnLabels;

	std::vector<int> m_iaRowCategories;
	std::vector<int> m_iaColumnCategories;

	std::vector<int> m_iaRowSubCategories;
	std::vector<int> m_iaColumnSubCategories;

	std::vector<double> m_faBin;
	std::vector<double> m_faBin2;
	std::vector<char> m_iaActive;

	bool m_bShowLabelTop;
	bool m_bShowLabelBottom;
	bool m_bShowLabelLeft;
	bool m_bShowLabelRight;

	bool m_bInvertColorScale;
	bool m_bInvertColorScale2;

	double m_fRangeCut;
	double m_fRangeCut2;

	bool m_b2DMode;

	bool m_bRangeSet;
	double m_fRangeMin;
	double m_fRangeMax;

	double m_fPlotExponent;
	double m_fPlotExponent2;
	bool m_bRangeZero;
	bool m_bRangeZero2;
	std::string m_sValueLabel;
	std::string m_sValueLabel2;
};




#endif



