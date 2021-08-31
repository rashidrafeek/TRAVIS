/*****************************************************************************

    libTRAVIS - Class library for trajectory analysis and visualization

    Copyright (C) 2015-2021 Martin Thomas

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


#ifndef LT_CONSTANTS_H
#define LT_CONSTANTS_H


#include "../config.h"


#ifdef NEW_CHARGEVAR


namespace Libtravis {

namespace Travis {

namespace Constants {
	const double pi = 3.1415926535897931;
	
	//! Atomic mass contant (CODATA 2014) in kg
	const double atomicMassUnit = 1.660539040e-27;
	//! Bohr radius (CODATA 2014) in m
	const double bohrRadius = 0.52917721067e-10;
	//! Boltzmann constant (CODATA 2014) in J/K
	const double boltzmannConstant = 1.38064852e-23;
	//! Elementary charge (CODATA 2014) in C
	const double elementaryCharge = 1.6021766208e-19;
	//! Speed of light in vacuum (CODATA 2014) in m/s
	const double speedOfLight = 299792458;
	
	//! Length conversion: atomic units -> picometers (= Bohr radius * 1e12)
	const double lengthAu2Pm = 52.917721067;
	//! Electric dipole moment conversion: Debye -> C*m (= 1e-21 / speed of light)
	const double dipoleDebye2Cm = 3.33564095e-30;
}

}

}

#endif

#endif


