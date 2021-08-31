/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020,  152 (16),  164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8),  2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Martin Brehm and Philipp di Dio.

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


/****************************************************************************

    Sources of Van-der-Waals radii:

    1) R. Scott Rowland, Robin Taylor: Intermolecular Nonbonded Contact Distances in Organic Crystal Structures: Comparison with Distances Expected from van der Waals Radii. In: J. Phys. Chem. 1996,  100, S. 7384-7391, doi:10.1021/jp953141+.
    2) A. Bondi: van der Waals Volumes and Radii. In: J. Phys. Chem. 1964, 68, S. 441-451, doi:10.1021/j100785a001.
    3) Manjeera Mantina, Adam C. Chamberlin, Rosendo Valero, Christopher J. Cramer, Donald G. Truhlar: Consistent van der Waals Radii for the Whole Main Group. In: J. Phys. Chem. A. 2009,  113, S. 5806-5812, doi:10.1021/jp8111556. 


    Coherent Neutron Scattering Cross Sections:

	Published in "Neutron News, Vol. 3, No. 3,  1992, pp. 29-37".
    Taken from "http://www.ncnr.nist.gov/resources/n-lengths/".



	Source of atomic masses: https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl

	The middle of the "standard atomic weight" intervals was chosen.

*****************************************************************************/


// This must always be the first include directive
#include "config.h"

#include "travis.h"


const char *GetRevisionInfo_elementdata(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_elementdata() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


void AddElementData()
{
	// Element, Ord. Number, Mass, Covalent Radius / pm, VdW Radius / pm, ( Coherent Neutron Scattering Cross Section / barn )

	// Lone Pair (for force field simulations)
	AddElement( "LP",  254,   3.00000,     1.0,   1.0 );


	// 1st Period
	AddElement( "H",     1,   1.00798,    37.0,  109.0,   1.7568 );
	AddElement( "D",     1,   2.01410,    37.0,  109.0 );
	AddElement( "He",    2,   4.00260,    32.0,  140.0 );


	// 2nd Period
	AddElement( "Li",    3,   6.9675,    134.0,  182.0 );
	AddElement( "Be",    4,   9.01218,    90.0,  153.0 );
	AddElement( "B",     5,  10.8135,     90.0,  192.0 );
	AddElement( "C",     6,  12.0106,     82.0,  170.0,   5.551 );
	AddElement( "N",     7,  14.0069,     77.0,  155.0,  11.01  );
	AddElement( "O",     8,  15.9994,     75.0,  152.0,   4.232 );
	AddElement( "F",     9,  18.99840,    73.0,  147.0,   4.017 );
	AddElement( "Ne",   10,  20.17976,    69.0,  154.0 );


	// 3rd Period
	AddElement( "Na",   11,  22.98977,    71.0,  227.0 );
	AddElement( "Mg",   12,  24.3055,    130.0,  173.0 );
	AddElement( "Al",   13,  26.98154,   154.0,  184.0 );
	AddElement( "Si",   14,  28.085,     118.0,  210.0 );
	AddElement( "P",    15,  30.97376,   111.0,  180.0,   3.307  );
	AddElement( "S",    16,  32.0675,    106.0,  180.0,   1.0186 );
	AddElement( "Cl",   17,  35.4515,    102.0,  175.0,  11.5257 );
	AddElement( "Ar",   18,  39.9481,     97.0,  188.0 );


	// 4th Period
	AddElement( "K",    19,  39.09831,   196.0,  275.0 );
	AddElement( "Ca",   20,  40.0784,    174.0,  231.0 );

	AddElement( "Sc",   21,  44.95591,   144.0,  211.0 );
	AddElement( "Ti",   22,  47.8671,    136.0,  176.0 );
	AddElement( "V",    23,  50.94151,   125.0,  161.0 );
	AddElement( "Cr",   24,  51.99616,   127.0,  154.0 );
	AddElement( "Mn",   25,  54.93804,   139.0,  152.0 );
	AddElement( "Fe",   26,  55.8452,    125.0,  151.0 );
	AddElement( "Co",   27,  58.93319,   126.0,  150.0 );
	AddElement( "Ni",   28,  58.69344,   121.0,  163.0 );
	AddElement( "Cu",   29,  63.5463,    138.0,  140.0 );
	AddElement( "Zn",   30,  65.382,     131.0,  139.0 );

	AddElement( "Ga",   31,  69.7231,    126.0,  187.0 );
	AddElement( "Ge",   32,  72.6308,    122.0,  211.0 );
	AddElement( "As",   33,  74.92160,   121.0,  185.0 );
	AddElement( "Se",   34,  78.9718,    116.0,  190.0 );
	AddElement( "Br",   35,  79.904,     114.0,  185.0 );
	AddElement( "Kr",   36,  83.7982,    110.0,  202.0 );


	// 5th Period
	AddElement( "Rb",   37,  85.46783,   211.0,  303.0 );
	AddElement( "Sr",   38,  87.621,     192.0,  249.0 );

	AddElement( "Y",    39,  88.90584,   162.0,  216.0 );
	AddElement( "Zr",   40,  91.2242,    148.0,  192.0 );
	AddElement( "Nb",   41,  92.90637,   137.0,  175.0 );
	AddElement( "Mo",   42,  95.951,     145.0,  167.0 );
	AddElement( "Tc",   43,  98.00,      131.0,  163.0 );
	AddElement( "Ru",   44,  101.072,    126.0,  161.0 );
	AddElement( "Rh",   45,  102.90550,  135.0,  161.0 );
	AddElement( "Pd",   46,  106.421,    131.0,  163.0 );
	AddElement( "Ag",   47,  107.86822,  153.0,  172.0 );
	AddElement( "Cd",   48,  112.4144,   148.0,  158.0 );

	AddElement( "In",   49,  114.8181,   144.0,  193.0 );
	AddElement( "Sn",   50,  118.7107,   141.0,  217.0 );
	AddElement( "Sb",   51,  121.7601,   138.0,  206.0 );
	AddElement( "Te",   52,  127.603,    135.0,  206.0 );
	AddElement( "I",    53,  126.90447,  133.0,  198.0 );
	AddElement( "Xe",   54,  131.2936,   130.0,  216.0 );


	// 6th Period
	AddElement( "Cs",   55,  132.90545,  225.0,  343.0 );
	AddElement( "Ba",   56,  137.3277,   198.0,  268.0 );

	AddElement( "La",   57,  138.90548,  169.0,  224.0 );

	AddElement( "Ce",   58,  140.1161,   204.0,  218.0 );
	AddElement( "Pr",   59,  140.90766,  203.0,  219.0 );
	AddElement( "Nd",   60,  144.2423,   201.0,  218.0 );
	AddElement( "Pm",   61,  145.00,     199.0,  220.0 );
	AddElement( "Sm",   62,  150.362,    198.0,  216.0 );
	AddElement( "Eu",   63,  151.9641,   198.0,  216.0 );
	AddElement( "Gd",   64,  157.253,    196.0,  216.0 );
	AddElement( "Tb",   65,  158.92535,  194.0,  213.0 );
	AddElement( "Dy",   66,  162.5001,   192.0,  214.0 );
	AddElement( "Ho",   67,  164.93033,  192.0,  211.0 );
	AddElement( "Er",   68,  167.2593,   189.0,  211.0 );
	AddElement( "Tm",   69,  168.93422,  190.0,  211.0 );
	AddElement( "Yb",   70,  173.0545,   187.0,  211.0 );
	AddElement( "Lu",   71,  174.96681,  187.0,  209.0 );

	AddElement( "Hf",   72,  178.492,    150.0,  191.0 );
	AddElement( "Ta",   73,  180.94788,  138.0,  175.0 );
	AddElement( "W",    74,  183.841,    146.0,  167.0 );
	AddElement( "Re",   75,  186.2071,   159.0,  164.0 );
	AddElement( "Os",   76,  190.233,    128.0,  162.0 );
	AddElement( "Ir",   77,  192.2173,   137.0,  163.0 );
	AddElement( "Pt",   78,  195.0849,   138.0,  175.0 );
	AddElement( "Au",   79,  196.96657,  144.0,  166.0 );
	AddElement( "Hg",   80,  200.5923,   149.0,  155.0 );

	AddElement( "Tl",   81,  204.384,    148.0,  196.0 );
	AddElement( "Pb",   82,  207.21,     146.0,  202.0 );
	AddElement( "Bi",   83,  208.98040,  146.0,  207.0 );
	AddElement( "Po",   84,  209.00,     140.0,  197.0 );
	AddElement( "At",   85,  210.00,     145.0,  202.0 );
	AddElement( "Rn",   86,  222.00,     145.0,  220.0 );


	/* 7th period */
	AddElement( "Fr",   87,  223.00,     260.0,  348.0 );
	AddElement( "Ra",   88,  226.00,     221.0,  283.0 );

	AddElement( "Ac",   89,  227.00,     215.0,  200.0 );

	AddElement( "Th",   90,  232.03774,  206.0,  200.0 );
	AddElement( "Pa",   91,  231.03588,  200.0,  196.0 );
	AddElement( "U",    92,  238.02891,  196.0,  186.0 );
	AddElement( "Np",   93,  237.00,     190.0,  186.0 );
	AddElement( "Pu",   94,  244.10,     187.0,  191.0 );
	AddElement( "Am",   95,  243.10,     180.0,  208.0 );
	AddElement( "Cm",   96,  247.10,     169.0,  209.0 );
	AddElement( "Bk",   97,  247.10,     160.0,  204.0 );
	AddElement( "Cf",   98,  251.10,     160.0,  223.0 );
	AddElement( "Es",   99,  254.10,     160.0,  223.0 );
	AddElement( "Fm",  100,  257.10,     160.0,  200.0 );
	AddElement( "Md",  101,  258.00,     160.0,  200.0 );
	AddElement( "No",  102,  259.00,     160.0,  200.0 );
	AddElement( "Lr",  103,  262.00,     160.0,  200.0 );
	AddElement( "Rf",  104,  267.00,     160.0,  200.0 );
	AddElement( "Db",  105,  268.00,     160.0,  200.0 );
	AddElement( "Sg",  106,  271.00,     160.0,  200.0 );
	AddElement( "Bh",  107,  272.00,     160.0,  200.0 );
	AddElement( "Hs",  108,  270.00,     160.0,  200.0 );
	AddElement( "Mt",  109,  276.00,     160.0,  200.0 );
	AddElement( "Ds",  110,  281.00,     160.0,  200.0 );
	AddElement( "Rg",  111,  280.00,     160.0,  200.0 );
	AddElement( "Cn",  112,  285.00,     160.0,  200.0 );
	AddElement( "Nh",  113,  284.00,     160.0,  200.0 );
	AddElement( "Fl",  114,  289.00,     160.0,  200.0 );
	AddElement( "Mc",  115,  288.00,     160.0,  200.0 );
	AddElement( "Lv",  116,  293.00,     160.0,  200.0 );
	AddElement( "Ts",  117,  292.00,     160.0,  200.0 );
	AddElement( "Og",  118,  294.00,     160.0,  200.0 );


	// Virtual Atom
	AddElement( "#",     0,    0.00,       0.0,    0.0 );

	// Colors for some common atoms. Other atoms have standard color.
	SetElementColor( "H",  255,  255,  255,  150,  150,  150 );
	SetElementColor( "D",  255,  255,  255,  150,  150,  150 );
//	SetElementColor( "C",  228,  113,   0,  180,  180,  180 ); // Orange
	SetElementColor( "B",  255,  181,  181,  200,  200,  200 );
	SetElementColor( "C",  128,  128,  128,  128,  128,  128 );
	SetElementColor( "N",    0,   0,  255,  140,  140,  140 );
	SetElementColor( "O",  255,   0,   0,  160,  160,  160 );
	SetElementColor( "F",    0,  255,   0,  160,  160,  160 );
	SetElementColor( "P",  228,  113,   0,  180,  180,  180 );
	SetElementColor( "S",  255,  255,   0,  200,  200,  200 );
	SetElementColor( "Cl",   0,  255,   0,  160,  160,  160 );
	SetElementColor( "Br",  166,  41,  41,  160,  160,  160 );
	SetElementColor( "I",  255,   0,  255,  160,  160,  160 );
}

