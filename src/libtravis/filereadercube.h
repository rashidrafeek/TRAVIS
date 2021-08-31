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


#ifndef LT_FILEREADERCUBE_H
#define LT_FILEREADERCUBE_H


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include "cellinfo.h"
#include "datagrid.h"
#include "filereader.h"
#include "vector.h"

#include <cstddef>
#include <functional>
#include <ios>
#include <string>
#include <vector>

namespace Libtravis {

namespace Travis {

class FileReaderCube : public FileReader
{
public:
	FileReaderCube() : FileReader() {}
	FileReaderCube(const std::string &filename) : FileReader(filename) {}
	
	void setCommentTargetAccessor(const std::function<std::string &()> &commentTargetAccessor) { m_commentTargetAccessor = commentTargetAccessor; }
	void setCellInfoTargetAccessor(const std::function<CellInfo &()> &cellInfoTargetAccessor) { m_cellInfoTargetAccessor = cellInfoTargetAccessor; }
	void setAtomicNumbersTargetAccessor(const std::function<std::vector<unsigned int> &()> &atomicNumbersTargetAccessor) { m_atomicNumbersTargetAccessor = atomicNumbersTargetAccessor; }
	void setAtomPositionsTargetAccessor(const std::function<std::vector<Vector3> &()> &atomPositionsTargetAccessor) { m_atomPositionsTargetAccessor = atomPositionsTargetAccessor; }
	void setDataGridTargetAccessor(const std::function<DataGrid3D<double> &()> &dataGridTargetAccessor) { m_dataGridTargetAccessor = dataGridTargetAccessor; }
	
	void setCoordinatesConversionFactor(double factor) { m_coordinatesConversionFactor = factor; }
	void setFastStringConversion(bool fastStringConversion) { m_fastStringConversion = fastStringConversion; }
	
private:
	std::size_t m_currentStep = 0;
	std::vector<std::streampos> m_stepPos;
	
	std::function<std::string &()> m_commentTargetAccessor;
	std::function<CellInfo &()> m_cellInfoTargetAccessor;
	std::function<std::vector<unsigned int> &()> m_atomicNumbersTargetAccessor;
	std::function<std::vector<Vector3> &()> m_atomPositionsTargetAccessor;
	std::function<DataGrid3D<double> &()> m_dataGridTargetAccessor;
	
	double m_coordinatesConversionFactor = 1.0;
	bool m_fastStringConversion = true;
	
	virtual bool p_hasNextStep() override;
	virtual bool p_readNextStep() override;
	virtual bool p_skipSteps(std::size_t numSteps) override;
	virtual void p_rewind() override;
	
	static double p_fastStringToDouble(const char *string, char **stringEnd);
};

}

}

#endif

#endif


