/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Sascha Gehrke.

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


#ifndef VOID_H
#define VOID_H


// This must always be the first include directive
#include "config.h"

#include "domain.h"
#include "moltools.h"

#ifdef USE_OMP
#include <omp.h>
#endif


bool gatherVoid();
bool initializeVoid();
bool processVoid(CTimeStep* ts);
bool finalizeVoid();



class CSphere: public CxObject
{
public:
	CSphere();
	~CSphere();

	void FindCenter(CxDVector3*,CxDVector3*,CxDVector3*,CxDVector3*);
	bool GrowIn(std::vector<CSphere*>*,double, double); //, double*);
		
	double m_dRadius;
	CxDVector3 m_vCenter;
	std::vector<int> m_iaTouch;
	int dom;
	bool m_bActive;
private:

};



class CVoidObservation: public CObservation
{
public:
	CVoidObservation();
	~CVoidObservation();

	std::vector< CxIntArray* > m_oaDomains;	//  Contains one IntArray with Offsets per Domain 
	double m_dMinBox;		//  Minimum size of sub-box
	double m_dBox[3];		//  Current size of sub-box

	double m_dCutOff;		//  CutOff (square) Distance to see to voids as seperated
	double m_dAV;			//  AV CutOff which is needed by a void domain to be counted as "spheric"
	double m_dMinRadius;		//  CutOff minimum radius to accept spheric void
	double m_dTemperature;		//  Temperature of the Simulation
	double m_dScale;		//  Scaling factor for the atom radii
	int m_iSpotGrid;		//  Number of Gridpoints per Dimension
	int m_iDepth;			//  Correlation depth for the fluidity analysis
	std::vector<CSphere*> m_oaVoids;	//  Contains temporary CSpheres of the Voids
	std::vector< std::vector<int> > VoidKey;  // Contains an int array for each box with the Void positions in m_oaVoids
	std::vector< std::vector<CSphere*> > spheres;
	std::vector<CDomain*> m_oaVoroVoids;	//  Contains temporary VoronoiCells of the Voids
	std::vector<CDomain*> m_oaDomVoids;	//  Contains temporary VoronoiDomains of the Voids
	std::vector< std::vector<CxIntArray*> > m_oaACFSpots;	//  Contains InterDom a CxIntArray per Controlpoint in which the timesteps of change are saved ( compare HB analysis )  THIS ONE IS NOT(!!!!) DOMAIN DEPENDED!!!!!!!
	std::vector<double> m_oaVolume;		//  Contains the total effective void volumes
	long int m_iPics;		//  Number of analyzed steps
	std::vector<double> m_daMasse;		//  Geometrical averaged mass of all molecules
	std::vector<double> m_daColRad;		//  Geometrical averaged collision radius of the molecules
	double m_dAveRad;		//  Average radius of the void spheres

	bool m_bMakePics;		//  Generate xyz files of voids in first step?
	bool m_bNoH;			//  Exclude hydrogen?
	bool m_bSurfTens;		//  Calculate Surface Tension?
	bool m_bViscosity;		//  Calculate Viscosity?
	bool m_bVoidSphDF;		//  Perform Spheric DF calculation?
	bool m_bVoidDomDF;		//  Perform Void domain analysis?
	bool m_bVoidACF;  		//  Perform Void Fluidity analysis?
	bool m_bHoleDF;			//  Perform Hole Theory Analysis?
	bool m_bFull;			//  Analyze whole cell? Adds another Domain with everthing not defined ("+1")
	int m_iInterDom;		//  Number of Intervals of domain decomposition series
	std::vector<CDF*>  m_oaVoidSphDF;	//  One for each Domain (last for border region "+1" for Full "+2")
	std::vector<CDF*>  m_oaVoidDomDF;
	std::vector<CDF*>  m_oaVoidAVDF;
	std::vector<C2DF*> m_oaVoidAV2DF;
	std::vector<CDF*>  m_oaHoleDF;

	void initialize();
	void process(CTimeStep* ts);
	void finalize();

	bool Start_Domainize();
	bool Domainize(CxIntArray*,CDomain*,CDomain*);
	bool CheckIfVoid(std::vector<int>,CxDVector3*,double);
	void AddAll(int,CxIntArray*);
	bool checkAtomChoice(const char*,int,CxIntArray*);
	void Add(int,int,int,CxIntArray*);
	void FindXYZ(int,int*,int*);

	double time_tot;
	double time_dom;
	double time_voro1;
	double time_voro2;
	double time_grow;
	double time_crit;
private:

};
#endif
