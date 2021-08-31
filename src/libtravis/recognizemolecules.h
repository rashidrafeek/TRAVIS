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


#ifndef LT_RECOGNIZEMOLECULES_H
#define LT_RECOGNIZEMOLECULES_H


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include "cellinfo.h"
#include "datacalculator.h"
#include "element.h"
#include "moleculeset.h"

#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Libtravis {

namespace Travis {

class RecognizeMolecules : public DataCalculator
{
public:
	void setAtomNamesAccessor(const std::function<std::vector<std::string> &()> &atomNamesAccessor) { m_atomNamesAccessor = atomNamesAccessor; }
	void setAtomicNumbersAccessor(const std::function<std::vector<unsigned int> &()> &atomicNumbersAccessor) { m_atomicNumbersAccessor = atomicNumbersAccessor; }
	void setAtomPositionsAccessor(const std::function<std::vector<Vector3> &()> &atomPositionsAccessor) { m_atomPositionsAccessor = atomPositionsAccessor; }
	void setCellInfoAccessor(const std::function<CellInfo &()> &cellInfoAccessor) { m_cellInfoAccessor = cellInfoAccessor; }
	void setMoleculeSetTargetAccessor(const std::function<MoleculeSet &()> &moleculeSetTargetAccessor) { m_moleculeSetTargetAccessor = moleculeSetTargetAccessor; }
	
	void addSpecialElement(const Element &element) {
		m_specialElements.push_back(element);
	}
	
	void setRadiusFactor(double factor) { m_radiusFactor = factor; }
	void setBondStretchingFactor(double factor) { m_bondStretchingFactor = factor; }
	void setNonHydrogenMassThreshold(double threshold) { m_nonHydrogenMassThreshold = threshold; }
	void setAtomCodeDifferenceThreshold(double threshold) { m_atomCodeDifferenceThreshold = threshold; }
	
	MoleculeSet recognizeMolecules(const std::vector<std::string> &atomNames, const std::vector<Vector3> &atomPositions, const CellInfo &cellInfo = CellInfo());
	MoleculeSet recognizeMolecules(const std::vector<unsigned int> &atomicNumbers, const std::vector<Vector3> &atomPositions, const CellInfo &cellInfo = CellInfo());
	
private:
	std::function<std::vector<std::string> &()> m_atomNamesAccessor;
	std::function<std::vector<unsigned int> &()> m_atomicNumbersAccessor;
	std::function<std::vector<Vector3> &()> m_atomPositionsAccessor;
	std::function<CellInfo &()> m_cellInfoAccessor;
	std::function<MoleculeSet &()> m_moleculeSetTargetAccessor;
	
	CellInfo m_cellInfo;
	std::vector<Element> m_specialElements;
	
	double m_radiusFactor = 1.0;
	double m_bondStretchingFactor = 1.15;
	double m_nonHydrogenMassThreshold = 4.5;
	double m_atomCodeDifferenceThreshold = 1e-3;
	
	virtual void p_calculateNextStep() override {
		if (!m_atomNamesAccessor)
			throw std::runtime_error("RecognizeMolecules: source for names not set");
		if (!m_atomPositionsAccessor)
			throw std::runtime_error("RecognizeMolecules: source for positions not set");
		if (!m_moleculeSetTargetAccessor)
			throw std::runtime_error("RecognizeMolecules: target for molecule set not set");
		
		if (!m_cellInfoAccessor)
			m_moleculeSetTargetAccessor() = recognizeMolecules(m_atomNamesAccessor(), m_atomPositionsAccessor());
		else
			m_moleculeSetTargetAccessor() = recognizeMolecules(m_atomNamesAccessor(), m_atomPositionsAccessor(), m_cellInfoAccessor());
	}
	
	void p_recognizeMolecules(MoleculeSet &moleculeSet, const std::vector<Vector3> &atomPositions);
	
	struct RecursionInfo
	{
		const std::vector<MoleculeSet::Atom> &atoms;
		const std::vector<Vector3> &atomPositions;
		std::vector<std::size_t> &atomMolIndex;
		std::vector<std::size_t> &helperMolecule;
		std::vector<std::pair<std::size_t, std::size_t>> &helperMoleculeBonds;
		std::vector<std::pair<std::size_t, std::size_t>> &helperMoleculeBondsCross;
	};
	
	std::unique_ptr<RecursionInfo> m_recursionInfo;
	
	void p_recognizeMoleculesRecursion(std::size_t beginAtom, std::size_t lastAtom, std::size_t beginAtomIndex);
	bool p_bondRange(const Vector3 &position1, const Vector3 &position2, const double radius1, const double radius2);
	
	struct RecursionInfoPath
	{
		const std::vector<std::vector<std::size_t>> &adjacencyList;
		std::vector<std::vector<std::size_t>> &ringList;
		std::size_t target;
		std::vector<std::size_t> &stack;
		std::vector<bool> &unused;
	};
	
	std::unique_ptr<RecursionInfoPath> m_recursionInfoPath;
	
	void p_findPathRecursion(std::size_t current);
};

}

}

#endif

#endif


