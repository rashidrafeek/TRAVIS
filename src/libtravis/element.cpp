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


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include "element.h"

#include <algorithm>
#include <iterator>
#include <string>

namespace Libtravis {

namespace Travis {

namespace Elements {

namespace {
const Element elements[] = {
	{ 1, "H", 1.01, 37.0, 110.0 },
	{ 6, "C", 12.01, 82.0, 170.0 },
	{ 7, "N", 14.01, 77.0, 155.0 },
	{ 8, "O", 16.00, 75.0, 152.0 }
};
}

const Element &findElement(unsigned int atomicNumber) {
	auto it = std::find_if(std::cbegin(elements), std::cend(elements), [&atomicNumber] (const Element &element) { return element.atomicNumber == atomicNumber; });
	if (it == std::cend(elements))
		throw UnknownElement();
	return *it;
}

const Element &findElement(const std::string &name) {
	auto it = std::find_if(std::cbegin(elements), std::cend(elements), [&name] (const Element &element) { return element.name == name; });
	if (it == std::cend(elements))
		throw UnknownElement();
	return *it;
}

}

}

}

#endif


