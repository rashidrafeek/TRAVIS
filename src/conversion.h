/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Martin Brehm and Martin Thomas.

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


#ifndef CONVERSION_H
#define CONVERSION_H


// This must always be the first include directive
#include "config.h"


// Physical constants
#define CONST_ATOMIC_MASS_UNIT (1.660538921E-27) // CODATA 2010, Atomic mass unit in kg
#define CONST_BOHR_MAGNETON (927.400968E-26) // CODATA 2010, Bohr magneton in A*m^2
#define CONST_BOHR_RADIUS (0.52917721092E-10) // CODATA 2010, Bohr radius in m
#define CONST_ELEMENTARY_CHARGE (1.602176565E-19) // CODATA 2010, Elementary charge in C
#define CONST_SPEED_OF_LIGHT (299792458.0) // CODATA 2010, Vacuum speed of light in m/s
#define CONST_AVOGADRO (6.02214085774E23) // Avogadro's number in 1/mol
#define CONST_EPSILON0 (8.8541878176E-12) // Vacuum permittivity in F/m
#define CONST_PLANCK (6.626070040E-34) // Planck constant
#define CONST_BOLTZMANN (1.38064852E-23) // Boltzmann Constant

// Length conversion
#define LEN_AU2PM (52.917721092) // Atomic units to pm: CONST_BOHR_RADIUS*1E12

// Time conversion
#define TIME_AU2FS (0.02418884326505) // Atomic units to fs

// Energy conversion (greetings to the Max Planck Institute :-D )
#define ENER_HARTREE2KJMOL (2625.4995)

// Dipole moment conversion
#define DIP_DEBYE2CM (3.335640952E-30) // Debye to C*m: 1E-21/CONST_SPEED_OF_LIGHT
#define DIP_EPM2DEBYE (0.04803204506) // Elementary charge*pm to Debye: CONST_ELEMENTARY_CHARGE*CONST_SPEED_OF_LIGHT*1E9

// Electric field strength conversion
#define EFIELD_AU2VPM (0.514220652) // Atomic unit of electric field "Eh / (e * a0)" to "V/pm"

// Total current conversion
#define CURR_AUFS2MBPM (9.142057808E-4) // Atomic unit of dipole moment/fs to Bohr magneton/pm: CONST_ELEMENTARY_CHARGE*CONST_BOHR_RADIUS*1E3/CONST_BOHR_MAGNETON
#define CURR_EMS2MBPM (1.727598547E-8) // Elementary charge*m/s to Bohr magneton/pm: CONST_ELEMENTARY_CHARGE*1E-12/CONST_BOHR_MAGNETON
#define CURR_MBPM2CMS (9.27400968114E-12) // Bohr magneton/pm to Coulomb*m/s: CONST_BOHR_MAGNETON/CONST_ELEMENTARY_CHARGE/1E-12*CONST_ELEMENTARY_CHARGE

// Magnetic moment conversion
#define MAG_AUPMFS2MB (9.142057808E-4) // Atomic unit of dipole moment*pm/fs to Bohr magneton: CONST_ELEMENTARY_CHARGE*CONST_BOHR_RADIUS*1E3/CONST_BOHR_MAGNETON
#define MAG_EPMMS2MB (1.727598547E-8) // Elementary charge*pm*m/s to Bohr magneton: CONST_ELEMENTARY_CHARGE*1E-12/CONST_BOHR_MAGNETON


#endif
