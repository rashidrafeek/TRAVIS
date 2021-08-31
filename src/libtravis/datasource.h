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


#ifndef LT_DATASOURCE_H
#define LT_DATASOURCE_H


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include <cstddef>

namespace Libtravis {

namespace Travis {

class DataSource
{
public:
	virtual ~DataSource() = default;
	
	bool hasNextStep() { return p_hasNextStep(); }
	bool readNextStep() { return p_readNextStep(); }
	bool skipSteps(std::size_t numSteps) { return p_skipSteps(numSteps); }
	void rewind() { return p_rewind(); }
	
private:
	virtual bool p_hasNextStep() = 0;
	virtual bool p_readNextStep() = 0;
	virtual bool p_skipSteps(std::size_t numSteps) = 0;
	virtual void p_rewind() = 0;
};

}

}

#endif

#endif


