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


#include "filereadercube.h"

#include "cellinfo.h"
#include "vector.h"

#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace Libtravis {

namespace Travis {

bool FileReaderCube::p_hasNextStep() {
	m_file.peek();
	return !m_file.eof();
}

bool FileReaderCube::p_readNextStep() {
	if (m_currentStep > 0 && m_stepPos.size() == m_currentStep - 1)
		m_stepPos.push_back(m_file.tellg());
	
	if (!hasNextStep())
		return false;
	
	std::string line;
	if (!std::getline(m_file, line)) {
		throw std::runtime_error("FileReaderCube: Could not read first comment line");
	}
	if (m_commentTargetAccessor)
		m_commentTargetAccessor() = std::move(line);
	if (!std::getline(m_file, line))
		throw std::runtime_error("FileReaderCube: Could not read second comment line");
	if (m_commentTargetAccessor) {
		m_commentTargetAccessor().append("\n");
		m_commentTargetAccessor().append(line);
	}
	
	if (!std::getline(m_file, line))
		throw std::runtime_error("FileReaderCube: Could not read line of atom number and cell origin");
	std::istringstream iss(line);
	std::size_t numAtoms;
	Vector3 origin;
	if (!(iss >> numAtoms >> origin[0] >> origin[1] >> origin[2]))
		throw std::runtime_error("FileReaderCube: Could not parse number of atoms and cell origin");
	
	std::size_t numPoints[3];
	Vector3 cellVectors[3];
	for (unsigned int i = 0; i < 3; ++i) {
		if (!std::getline(m_file, line))
			throw std::runtime_error("FileReaderCube: Could not read line of cell vector " + std::to_string(i + 1));
		iss.str(line);
		iss.seekg(0);
		if (!(iss >> numPoints[i] >> cellVectors[i][0] >> cellVectors[i][1] >> cellVectors[i][2]))
			throw std::runtime_error("FileReaderCube: Could not parse line of cell vector " + std::to_string(i + 1));
	}
	
	if (m_cellInfoTargetAccessor) {
		CellInfo cellInfo(cellVectors[0] * m_coordinatesConversionFactor * numPoints[0], cellVectors[1] * m_coordinatesConversionFactor * numPoints[1], cellVectors[2] * m_coordinatesConversionFactor * numPoints[2]);
		cellInfo.setOrigin(origin);
		m_cellInfoTargetAccessor() = std::move(cellInfo);
	}
	
	if (m_atomicNumbersTargetAccessor) {
		m_atomicNumbersTargetAccessor().clear();
		m_atomicNumbersTargetAccessor().reserve(numAtoms);
	}
	if (m_atomPositionsTargetAccessor) {
		m_atomPositionsTargetAccessor().clear();
		m_atomPositionsTargetAccessor().reserve(numAtoms);
	}
	
	for (std::size_t i = 0; i < numAtoms; ++i) {
		if (!std::getline(m_file, line))
			throw std::runtime_error("FileReaderCube: Could not read line of atom " + std::to_string(i + 1));
		iss.str(line);
		iss.seekg(0);
		unsigned int atomicNumber;
		double charge;
		double v[3];
		if (!(iss >> atomicNumber >> charge >> v[0] >> v[1] >> v[2]))
			throw std::runtime_error("FileReaderCube: Could not parse line of atom " + std::to_string(i + 1));
		
		if (m_atomicNumbersTargetAccessor) {
			m_atomicNumbersTargetAccessor().push_back(atomicNumber);
		}
		if (m_atomPositionsTargetAccessor) {
			m_atomPositionsTargetAccessor().emplace_back(v[0] * m_coordinatesConversionFactor, v[1] * m_coordinatesConversionFactor, v[2] * m_coordinatesConversionFactor);
		}
	}
	
	if (m_dataGridTargetAccessor) {
		DataGrid3D<double> &dataGridTarget = m_dataGridTargetAccessor();
		std::size_t totalCount = 0;
		const double conversionFactorCube = m_coordinatesConversionFactor * m_coordinatesConversionFactor * m_coordinatesConversionFactor;
		dataGridTarget.setSizes(numPoints[0], numPoints[1], numPoints[2]);
		for (std::size_t i = 0; i < numPoints[0]; ++i) {
			for (std::size_t j = 0; j < numPoints[1]; ++j) {
				for (std::size_t k = 0; k < numPoints[2] / 6; ++k) {
					if (!std::getline(m_file, line))
						throw std::runtime_error("FileReaderCube: Could not read line of grid data");
					const char *p = line.c_str();
					char *q;
					for (std::size_t l = 0; l < 6; ++l) {
						dataGridTarget.m_data[totalCount] = (m_fastStringConversion ? p_fastStringToDouble(p, &q) : std::strtod(p, &q)) / conversionFactorCube;
						if (p == q)
							throw std::runtime_error("FileReaderCube: Could not parse line of grid data");
						p = q;
						++totalCount;
					}
				}
				if (numPoints[2] % 6 != 0) {
					if (!std::getline(m_file, line))
						throw std::runtime_error("FileReaderCube: Could not read line of grid data");
					const char *p = line.c_str();
					char *q;
					for (std::size_t l = 0; l < numPoints[2] % 6; ++l) {
						dataGridTarget.m_data[totalCount] = (m_fastStringConversion ? p_fastStringToDouble(p, &q) : std::strtod(p, &q)) / conversionFactorCube;
						if (p == q)
							throw std::runtime_error("FileReaderCube: Could not parse line of grid data");
						p = q;
						dataGridTarget.m_data[totalCount] /= conversionFactorCube;
						++totalCount;
					}
				}
			}
		}
	} else {
		for (std::size_t i = 0; i < numPoints[0]; ++i) {
			for (std::size_t j = 0; j < numPoints[1]; ++j) {
				std::size_t numLines = numPoints[2] / 6;
				if (numPoints[2] % 6 != 0)
					++numLines;
				for (std::size_t k = 0; k < numLines; ++k)
					if (!m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n'))
						throw std::runtime_error("FileReaderCube: Could not skip line of grid data");
			}
		}
	}
	
	return true;
}

bool FileReaderCube::p_skipSteps(std::size_t numSteps) {
	if (m_stepPos.size() > m_currentStep + numSteps) {
		m_currentStep += numSteps;
		m_file.seekg(m_stepPos[m_currentStep - 1]);
	} else {
		if (m_stepPos.size() > m_currentStep) {
			numSteps -= m_stepPos.size() - m_currentStep - 1;
			m_currentStep += m_stepPos.size() - m_currentStep - 1;
			m_file.seekg(m_stepPos[m_currentStep - 1]);
		}
		for (std::size_t i = 0; i < numSteps; ++i) {
			if (m_currentStep > 0 && m_stepPos.size() == m_currentStep - 1)
				m_stepPos.push_back(m_file.tellg());
			
			if (!hasNextStep())
				return false;
			
			std::string line;
			if (!std::getline(m_file, line)) {
				throw std::runtime_error("FileReaderCube: Could not read first comment line");
			}
			if (!m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n'))
				throw std::runtime_error("FileReaderCube: Could not skip second comment line");
			
			if (!std::getline(m_file, line))
				throw std::runtime_error("FileReaderCube: Could not read line of atom number");
			std::istringstream iss(line);
			std::size_t numAtoms;
			if (!(iss >> numAtoms))
				throw std::runtime_error("FileReaderCube: Could not parse number of atoms");
			
			std::size_t numPoints[3];
			for (unsigned int j = 0; j < 3; ++j) {
				if (!std::getline(m_file, line))
					throw std::runtime_error("FileReaderCube: Could not read line of cell vector " + std::to_string(i + 1));
				iss.str(line);
				iss.seekg(0);
				if (!(iss >> numPoints[j]))
					throw std::runtime_error("FileReaderCube: Could not parse line of cell vector " + std::to_string(i + 1));
			}
			
			for (std::size_t j = 0; j < numAtoms; ++j) {
				if (!m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n'))
					throw std::runtime_error("FileReaderCube: Could not skip line of atom " + std::to_string(i + 1));
			}
			
			for (std::size_t j = 0; j < numPoints[0]; ++j) {
				for (std::size_t k = 0; k < numPoints[1]; ++k) {
					std::size_t numLines = numPoints[2] / 6;
					if (numPoints[2] % 6 != 0)
						++numLines;
					for (std::size_t l = 0; l < numLines; ++l)
						if (!m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n'))
							throw std::runtime_error("FileReaderCube: Could not skip line of grid data");
				}
			}
			
			++m_currentStep;
		}
	}
	
	return true;
}

void FileReaderCube::p_rewind() {
	rewindFile();
	m_currentStep = 0;
}

double FileReaderCube::p_fastStringToDouble(const char *string, char **stringEnd) {
	while (*string == ' ' || (*string >= 0x9 && *string <= 0xD))
		++string;
	
	bool negative = false;
	if (*string == '-') {
		negative = true;
		++string;
	} else if (*string == '+') {
		++string;
	}
	
	double result = 0.0;
	while (*string >= '0' && *string <= '9') {
		result = result * 10.0 + static_cast<double>(*string - '0');
		++string;
	}
	
	if (*string == '.') {
		++string;
		const char *begin = string;
		while (*string >= '0' && *string <= '9')
			++string;
		const char *end = string;
		double fracPart = 0.0;
		while (string != begin) {
			--string;
			fracPart = fracPart / 10.0 + static_cast<double>(*string - '0');
		}
		result += fracPart / 10.0;
		string = end;
	}
	
	if (*string == 'e' || *string == 'E') {
		++string;
		bool negExp = false;
		if (*string == '-') {
			negExp = true;
			++string;
		} else if (*string == '+') {
			++string;
		}
		unsigned int exponent = 0;
		while (*string >= '0' && *string <= '9') {
			exponent = exponent * 10 + static_cast<unsigned int>(*string - '0');
			++string;
		}
		double base = 10.0;
		double expFactor = 1.0;
		while (exponent != 0) {
			if (exponent % 2 == 1)
				expFactor *= base;
			base *= base;
			exponent /= 2;
		}
		if (negExp)
			expFactor = 1.0 / expFactor;
		result *= expFactor;
	}
	
	*stringEnd = const_cast<char *>(string);
	return negative ? -result : result;
}

}

}

#endif


