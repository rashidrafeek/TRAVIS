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


#ifndef SANKEY_H
#define SANKEY_H


// This must always be the first include directive
#include "config.h"

#include <vector>



class CSankeyNode {
public:
	std::vector<double> m_faValues;
	std::vector<double> m_faSubPos;
	std::vector<int> m_iaConnectIndex;
	double m_fSum;
	unsigned long m_iColor;

	double m_fPosStart;
	double m_fPosEnd;
};



class CSankeyDiagramEngine {
public:

	CSankeyDiagramEngine();

	void BuildSankeyDiagram( const char *s = NULL );

	void WriteSVG( const char *s );

	int m_iType; // 0 - Transfer, 1 - Circular

	int m_iNodeCount;

	double m_fTotalSum;

	std::vector<double> m_faMatrix;
	std::vector<unsigned long> m_iaColors;
	std::vector<int> m_iaGroupPermutation;


	std::vector<double> m_faMatrixSource;
	std::vector<double> m_faMatrixSink;
	std::vector<int> m_iaGroupPermutationSource;
	std::vector<int> m_iaGroupPermutationSink;
	std::vector<unsigned long> m_iaColorsSource;
	std::vector<unsigned long> m_iaColorsSink;
	int m_iSourceCount;
	int m_iSinkCount;
	double m_fTotalSumSource;
	double m_fTotalSumSink;


	bool m_bSmooth;
	bool m_bSmoothHSV;
	bool m_bCorrectWidth;
	double m_fPadding;
	double m_fZeroThreshold;
	bool m_bCorrectRatio;

private:
	void WriteSVG_Circular( const char *s );
	void WriteSVG_Block( const char *s );

	std::vector<CSankeyNode> m_oaNodes;

	std::vector<CSankeyNode> m_oaNodesSource;
	std::vector<CSankeyNode> m_oaNodesSink;
};



#endif





