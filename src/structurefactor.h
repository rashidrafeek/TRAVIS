/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Martin Thomas.

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


#ifndef STRUCTUREFACTOR_H
#define STRUCTUREFACTOR_H


// This must always be the first include directive
#include "config.h"

#include "xintarray.h"
#include "xobarray.h"
#include "xobject.h"
#include "xstring.h"


class CDF;
class CTimeStep;


class CIsotope: public CxObject
{
public:
	friend class CStructureFactor;
	
	CIsotope(const char *label, double neutronFactor, double cma1, double cma2, double cma3, double cma4, double cmb1, double cmb2, double cmb3, double cmb4, double cmc);
	~CIsotope();
	
	const char *label() { return _label; }
	double neutronFactor() { return _neutronFactor; }
	double xrayFactor(double q);
	
private:
	char *_label;
	double _neutronFactor;
	double _cma[4];
	double _cmb[4];
	double _cmc;
};

// class CSFacObservation: public CxObject
// {
// public:
// 	CSFacObservation(bool global = false);
// 	~CSFacObservation();
// 	
// 	void initialize();
// 	void process(CTimeStep *ts);
// 	void finalize();
// 	
// private:
// 	bool _global;
// 	char *_name;
// 	
// 	CxObArray _agList;
// 	CxIntArray _atomIndex;
// 	CxIntArray _isotopeList;
// 	CxObArray _isotopeTypeList;
// 	CxIntArray _isotopeTypeCount;
// 	
// 	double _rdfMax;
// 	int _rdfRes;
// 	double _sfacMax;
// 	int _sfacRes;
// 	CxObArray _rdfList;
// 	
// 	bool _normFFac;
// 	int _normFFacFormula;
// 	bool _sharpening;
// 	CIsotope *_sharpeningIsotope;
// };

class CStructureFactorGroup: public CxObject
{
public:
	friend class CStructureFactorCross;
	
	CStructureFactorGroup(bool global, CxObArray &isotopeAssignList);
	~CStructureFactorGroup();
	
	void initialize(double rdfMax, int rdfRes);
	void process(CTimeStep *ts);
	void finalize(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate, bool m_bLorch);
	void finalizeIntra(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate, bool m_bLorch);
	void finalizeInter(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate, bool m_bLorch);
	
	CxString &getName() { return m_name; }
	bool sepInterIntra() const { return m_sepInterIntra; }
	CxObArray &getIsotopeTypeList() { return m_isotopeTypeList; }
	CxIntArray &getIsotopeTypeTotalCount() { return m_isotopeTypeTotalCount; }
	
private:
	CxString m_name;
	bool m_sepInterIntra;
	bool m_global;
	
	CxObArray m_atomGroupList;
	CxIntArray m_atomIndexList;
	CxIntArray m_singleMolList;
	CxIntArray m_isotopeList;
	CxObArray m_isotopeTypeList;
	CxIntArray m_isotopeTypeCount;
	CxIntArray m_isotopeTypeTotalCount;
	CxObArray m_rdfList;
	CxObArray m_rdfIntraList;
	CxObArray m_rdfInterList;
//>>>OMP
#ifdef USE_OMP
	std::vector< std::vector<CDF*> > m_tempLists;
	std::vector< std::vector<CDF*> > m_tempInterLists;
	std::vector< std::vector<CDF*> > m_tempIntraLists;
#endif
//<<<OMP
};

class CStructureFactorCross: public CxObject
{
public:
	CStructureFactorCross(CStructureFactorGroup *sfacGroup1, CStructureFactorGroup *sfacGroup2);
	~CStructureFactorCross();
	
	void initialize(double rdfMax, int rdfRes);
	void process(CTimeStep *ts);
	void finalize(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate, bool m_bLorch);
	void finalizeIntra(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate, bool m_bLorch);
	void finalizeInter(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate, bool m_bLorch);
	
	CxString &getName() { return m_name; }
	bool sepInterIntra() const { return m_sepInterIntra; }
	void setSepInterIntra(bool sep) { m_sepInterIntra = sep; }
	
private:
	CxString m_name;
	bool m_sepInterIntra;
	
	CStructureFactorGroup *m_sfacGroup1;
	CStructureFactorGroup *m_sfacGroup2;
	CxObArray m_rdfList;
	CxObArray m_rdfIntraList;
	CxObArray m_rdfInterList;
};



class CStructureFactor: public CxObject
{
public:
	CStructureFactor();
	~CStructureFactor();
	
	void initialize();
	void process(CTimeStep *ts);
	void finalize();
	void HandleGroup( const char *name, CStructureFactorGroup *sfgroup, CStructureFactorCross *sfcross, int mode );
	
private:
	double m_rdfMax;
	int m_rdfRes;
	double m_sfacMax;
	int m_sfacRes;
	//int m_normalization;
	bool m_sharpening;
	CIsotope *m_sharpeningIsotope;
	bool m_saveIntermediate;

	double m_fNormalizationFactor;

//>>> SG 27.08.2019
	bool m_bLorch;
//<<< SG 27.08.2019
	
	CStructureFactorGroup *m_globalSFac;
	CxObArray m_sFacGroups;
	CxObArray m_sFacCrosses;
};

bool gatherStructureFactor();
bool initializeStructureFactor();
void processStructureFactor(CTimeStep *ts);
void finalizeStructureFactor();

// class CStructureFactor : public CxObject
// {
// public:
// 	void TransformRDF(CDF *pin, CDF *pout);
// 	void Finish();
// 	void ProcessStep(CTimeStep *ts);
// 	void Parse();
// 	void Create();
// 	CStructureFactor();
// 	~CStructureFactor();
// 
// 	bool m_bDumpTotalRDF;
// 	bool m_bDumpElementRDFs;
// 	bool m_bDumpElementSFac;
// 	CxObArray m_oaRDFs;
// 	CxObArray m_oaSFacs;
// 	int m_iRDFRes, m_iSQRes;
// 	double m_fRDFRange;
// 	double m_fSQRange;
// };


#endif
