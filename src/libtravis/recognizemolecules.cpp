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


#include "recognizemolecules.h"

#include "element.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <stdexcept>

namespace Libtravis {

namespace Travis {

MoleculeSet RecognizeMolecules::recognizeMolecules(const std::vector<std::string> &atomNames, const std::vector<Vector3> &atomPositions, const CellInfo &cellInfo) {
	m_cellInfo = cellInfo;
	MoleculeSet moleculeSet;
	
	std::vector<std::size_t> sortPermutation(atomNames.size());
	std::iota(sortPermutation.begin(), sortPermutation.end(), 0);
	std::sort(sortPermutation.begin(), sortPermutation.end(), [&atomNames] (std::size_t a, std::size_t b) { return atomNames[a] < atomNames[b]; });
	
	// Create list of atom types and atoms
	std::size_t countType = 0;
	std::size_t last = 0;
	for (std::size_t i = 0; i < atomNames.size(); ++i) {
		const std::string &atomName = atomNames[sortPermutation[i]];
		if (atomName != atomNames[sortPermutation[last]] || i == 0) {
			auto it = std::find_if(m_specialElements.cbegin(), m_specialElements.cend(), [atomName] (const Element &element) { return element.name == atomName; });
			if (it != m_specialElements.cend()) {
				moleculeSet.m_atomTypes.emplace_back(MoleculeSet::AtomType(&moleculeSet, countType, atomName, *it, m_radiusFactor));
			} else {
				try {
					const Element &el = Elements::findElement(atomName);
					moleculeSet.m_atomTypes.emplace_back(MoleculeSet::AtomType(&moleculeSet, countType, atomName, el, m_radiusFactor));
				} catch(const Elements::UnknownElement &) {
					moleculeSet.m_atomTypes.emplace_back(MoleculeSet::AtomType(&moleculeSet, countType, atomName));
				}
			}
			++countType;
			last = i;
		}
		moleculeSet.m_atoms.emplace_back(MoleculeSet::Atom(&moleculeSet, sortPermutation[i], countType - 1));
	}
	std::sort(moleculeSet.m_atoms.begin(), moleculeSet.m_atoms.end(), [] (const MoleculeSet::Atom &a, const MoleculeSet::Atom &b) { return a.getId() < b.getId(); });
	
	p_recognizeMolecules(moleculeSet, atomPositions);
	
	return moleculeSet;
}

MoleculeSet RecognizeMolecules::recognizeMolecules(const std::vector<unsigned int> &atomicNumbers, const std::vector<Vector3> &atomPositions, const CellInfo &cellInfo) {
	m_cellInfo = cellInfo;
	MoleculeSet moleculeSet;
	
	std::vector<std::size_t> sortPermutation(atomicNumbers.size());
	std::iota(sortPermutation.begin(), sortPermutation.end(), 0);
	std::sort(sortPermutation.begin(), sortPermutation.end(), [&atomicNumbers] (std::size_t a, std::size_t b) { return atomicNumbers[a] < atomicNumbers[b]; });
	
	// Create list of atom types and atoms
	std::size_t countType = 0;
	std::size_t last = 0;
	for (std::size_t i = 0; i < atomicNumbers.size(); ++i) {
		unsigned int atomicNumber = atomicNumbers[sortPermutation[i]];
		if (atomicNumber != atomicNumbers[sortPermutation[last]] || i == 0) {
			try {
				const Element &el = Elements::findElement(atomicNumber);
				moleculeSet.m_atomTypes.emplace_back(MoleculeSet::AtomType(&moleculeSet, countType, el.name, el, m_radiusFactor));
			} catch(Elements::UnknownElement &) {
				throw std::runtime_error("RecognizeMolecules: Invalid atomic number");
			}
			++countType;
			last = i;
		}
		moleculeSet.m_atoms.emplace_back(MoleculeSet::Atom(&moleculeSet, sortPermutation[i], countType - 1));
	}
	std::sort(moleculeSet.m_atoms.begin(), moleculeSet.m_atoms.end(), [] (const MoleculeSet::Atom &a, const MoleculeSet::Atom &b) { return a.getId() < b.getId(); });
	
	p_recognizeMolecules(moleculeSet, atomPositions);
	
	return moleculeSet;
}

void RecognizeMolecules::p_recognizeMolecules(MoleculeSet &moleculeSet, const std::vector<Vector3> &atomPositions) {
	// Recursively find molecules
	std::vector<std::size_t> atomMolIndex(moleculeSet.m_atoms.size(), 0);
	std::size_t countMolType = 0;
	std::size_t countMol = 0;
	for (std::size_t i = 0; i < atomMolIndex.size(); ++i) {
		if (atomMolIndex[i] == 0) {
			atomMolIndex[i] = countMol + 1;
			
			std::vector<std::size_t> helperMolecule;
			std::vector<std::pair<std::size_t, std::size_t>> helperMoleculeBonds;
			std::vector<std::pair<std::size_t, std::size_t>> helperMoleculeBondsCross;
			
			helperMolecule.push_back(i);
			
			m_recursionInfo = std::make_unique<RecursionInfo>(RecursionInfo{moleculeSet.m_atoms, atomPositions, atomMolIndex, helperMolecule, helperMoleculeBonds, helperMoleculeBondsCross});
			
			p_recognizeMoleculesRecursion(i, i, 0);
			
			m_recursionInfo.reset();
			
			// Calculate initial atom codes
			std::vector<double> atomCodes(helperMolecule.size());
			for (std::size_t j = 0; j < helperMolecule.size(); ++j)
				atomCodes[j] = 10.0 * moleculeSet.m_atoms[helperMolecule[j]].getAtomType().getMass();
			for (std::pair<std::size_t, std::size_t> b: helperMoleculeBonds) {
				if (moleculeSet.m_atoms[helperMolecule[b.first]].getAtomType().getMass() > m_nonHydrogenMassThreshold)
					atomCodes[b.second] += 1.0;
				if (moleculeSet.m_atoms[helperMolecule[b.second]].getAtomType().getMass() > m_nonHydrogenMassThreshold)
					atomCodes[b.first] += 1.0;
			}
			for (std::pair<std::size_t, std::size_t> b: helperMoleculeBondsCross) {
				if (moleculeSet.m_atoms[helperMolecule[b.first]].getAtomType().getMass() > m_nonHydrogenMassThreshold)
					atomCodes[b.second] += 1.0;
				if (moleculeSet.m_atoms[helperMolecule[b.second]].getAtomType().getMass() > m_nonHydrogenMassThreshold)
					atomCodes[b.first] += 1.0;
			}
			
			// Iterate until number of different atom codes does not change anymore
			std::size_t count1 = 1, count2;
			std::vector<double> tempAtomCodes = atomCodes;
			std::sort(tempAtomCodes.begin(), tempAtomCodes.end());
			auto it = tempAtomCodes.cbegin();
			auto lastIt = tempAtomCodes.cbegin();
			while (++it != tempAtomCodes.cend()) {
				if (std::fabs(*it - *lastIt) > m_atomCodeDifferenceThreshold) {
					++count1;
					lastIt = it;
				}
			}
			do {
				for (std::size_t j = 0; j < atomCodes.size(); ++j)
					tempAtomCodes[j] = 5.0 * atomCodes[j];
				for (std::pair<std::size_t, std::size_t> b: helperMoleculeBonds) {
					tempAtomCodes[b.first] += atomCodes[b.second];
					tempAtomCodes[b.second] += atomCodes[b.first];
				}
				for (std::pair<std::size_t, std::size_t> b: helperMoleculeBondsCross) {
					tempAtomCodes[b.first] += atomCodes[b.second];
					tempAtomCodes[b.second] += atomCodes[b.first];
				}
				atomCodes = tempAtomCodes;
				std::sort(tempAtomCodes.begin(), tempAtomCodes.end());
				auto it2 = tempAtomCodes.cbegin();
				auto lastIt2 = tempAtomCodes.cbegin();
				count2 = count1;
				count1 = 1;
				while (++it2 != tempAtomCodes.cend()) {
					if (std::fabs(*it2 - *lastIt2) > m_atomCodeDifferenceThreshold) {
						++count1;
						lastIt2 = it2;
					}
				}
			} while (count2 != count1);
			
			// Count atom types and create sum formula
			std::vector<std::size_t> atomTypeCount(moleculeSet.m_atomTypes.size(), 0);
			std::string name;
			double mass = 0.0;
			
			std::vector<std::size_t> perm(helperMolecule.size());
			std::iota(perm.begin(), perm.end(), 0);
			std::sort(perm.begin(), perm.end(), [&moleculeSet, &helperMolecule] (std::size_t a, std::size_t b) { return moleculeSet.m_atoms[helperMolecule[a]].getAtomType().getId() < moleculeSet.m_atoms[helperMolecule[b]].getAtomType().getId(); });
			
			auto last = perm.begin();
			std::size_t lastAt = moleculeSet.m_atoms[helperMolecule[perm[0]]].getAtomType().getId();
			for (auto it2 = perm.begin(); it2 != perm.end(); ++it2) {
				std::size_t at = moleculeSet.m_atoms[helperMolecule[*it2]].getAtomType().getId();
				if (at != lastAt) {
					atomTypeCount[lastAt] = static_cast<std::size_t>(std::distance(last, it2));
					name.append(moleculeSet.m_atomTypes[lastAt].getName());
					if (std::distance(last, it2) > 1)
						name.append(std::to_string(std::distance(last, it2)));
					mass += moleculeSet.m_atomTypes[lastAt].getMass() * std::distance(last, it2);
					std::sort(last, it2, [&atomCodes] (std::size_t a, std::size_t b) { return atomCodes[a] > atomCodes[b]; });
					last = it2;
					lastAt = at;
				}
			}
			atomTypeCount[lastAt] = static_cast<std::size_t>(std::distance(last, perm.end()));
			name.append(moleculeSet.m_atomTypes[lastAt].getName());
			if (std::distance(last, perm.end()) > 1)
				name.append(std::to_string(std::distance(last, perm.end())));
			mass += moleculeSet.m_atomTypes[lastAt].getMass() * std::distance(last, perm.end());
			std::sort(last, perm.end(), [&atomCodes] (std::size_t a, std::size_t b) { return atomCodes[a] > atomCodes[b]; });
			
			// Create final lists of atoms and atom codes
			std::vector<double> atomCodeList(perm.size());
			std::vector<std::size_t> atomList(perm.size());
			for (std::size_t j = 0; j < perm.size(); ++j) {
				atomCodeList[j] = atomCodes[perm[j]];
				atomList[j] = helperMolecule[perm[j]];
			}
			
			// Check if molecule type is new
			auto it2 = std::find_if(moleculeSet.m_moleculeTypes.cbegin(), moleculeSet.m_moleculeTypes.cend(), [&name] (const MoleculeSet::MoleculeType &moleculeType) { return moleculeType.getName() == name; });
			bool exists = it2 != moleculeSet.m_moleculeTypes.cend();
			std::string name2 = name;
			if (exists) {
				char appendLetter = 'a';
				for(;;) {
					if (!std::equal(atomCodeList.cbegin(), atomCodeList.cend(), it2->m_atomCodes.cbegin(), [this] (double a, double b) { return !(std::fabs(a - b) > m_atomCodeDifferenceThreshold); })) {
						exists = false;
						name2 = name;
						name2.push_back(':');
						name2.push_back(appendLetter);
					}
					if (exists) {
						break;
					} else {
						it2 = std::find_if(moleculeSet.m_moleculeTypes.cbegin(), moleculeSet.m_moleculeTypes.cend(), [&name2] (const MoleculeSet::MoleculeType &moleculeType) { return moleculeType.getName() == name2; });
						if (it2 == moleculeSet.m_moleculeTypes.cend())
							break;
					}
					++appendLetter;
					if (appendLetter > 'z')
						throw std::runtime_error("RecognizeMolecules: Too many molecule types with the same name \"" + name + "\"");
				}
			}
			if (!exists) {
				// Create neighbor lists
				std::vector<std::size_t> invPerm(perm.size());
				for (std::size_t j = 0; j < perm.size(); ++j)
					invPerm[perm[j]] = j;
				
				std::vector<std::vector<std::size_t>> adjacencyList(helperMolecule.size());
				for (std::pair<std::size_t, std::size_t> b: helperMoleculeBonds) {
					adjacencyList[invPerm[b.first]].push_back(invPerm[b.second]);
					adjacencyList[invPerm[b.second]].push_back(invPerm[b.first]);
				}
				
				// Find all rings
				std::vector<std::vector<std::size_t>> ringList;
				while (!helperMoleculeBondsCross.empty()) {
					const std::pair<std::size_t, std::size_t> &b = helperMoleculeBondsCross.back();
					
					std::vector<std::size_t> stack(1, invPerm[b.first]);
					std::vector<bool> unused(helperMolecule.size(), true);
					unused[invPerm[b.first]] = false;
					m_recursionInfoPath = std::make_unique<RecursionInfoPath>(RecursionInfoPath{adjacencyList, ringList, invPerm[b.second], stack, unused});
					
					p_findPathRecursion(invPerm[b.first]);
					
					m_recursionInfoPath.reset();
					
					adjacencyList[invPerm[b.first]].push_back(invPerm[b.second]);
					adjacencyList[invPerm[b.second]].push_back(invPerm[b.first]);
					helperMoleculeBondsCross.pop_back();
				}
				
				// Erase rings that fully contain a smaller ring
				std::vector<std::size_t> eraseList;
				eraseList.reserve(ringList.size());
				for (std::size_t j = 0; j < ringList.size(); ++j) {
					for (std::size_t k = j + 1; k < ringList.size(); ++k) {
						if (std::find(eraseList.cbegin(), eraseList.cend(), j) != eraseList.cend())
							continue;
						if (std::find(eraseList.cbegin(), eraseList.cend(), k) != eraseList.cend())
							continue;
						bool contained = true;
						std::size_t j2;
						std::size_t k2;
						if (ringList[j].size() <= ringList[k].size()) {
							j2 = j;
							k2 = k;
						} else {
							j2 = k;
							k2 = j;
						}
						for (std::size_t l: ringList[j2]) {
							if (std::find(ringList[k2].cbegin(), ringList[k2].cend(), l) == ringList[k2].cend()) {
								contained = false;
								break;
							}
						}
						if (contained) {
							eraseList.push_back(k2);
						}
					}
				}
				
				std::sort(eraseList.begin(), eraseList.end());
				
				std::vector<std::size_t> permList;
				permList.reserve(ringList.size() - eraseList.size());
				for (std::size_t j = 0, k = 0; j < ringList.size(); ++j) {
					if (k < eraseList.size() && j == eraseList[k]) {
						++k;
					} else {
						permList.push_back(j);
					}
				}
				
				// Sort rings by size
				std::sort(permList.begin(), permList.end(), [&ringList] (std::size_t a, std::size_t b) { return ringList[a].size() < ringList[b].size(); });
				
				std::vector<std::vector<std::size_t>> ringListReduced(permList.size());
				for (std::size_t j = 0; j < permList.size(); ++j)
					ringListReduced[j] = std::move(ringList[permList[j]]);
				
				// Rotate each ring
				for (std::vector<std::size_t> &r: ringListReduced) {
					auto it3 = std::min_element(r.begin(), r.end());
					std::rotate(r.begin(), it3, r.end());
				}
				
				moleculeSet.m_moleculeTypes.emplace_back(MoleculeSet::MoleculeType(&moleculeSet, countMolType, name2, atomTypeCount, atomCodeList, adjacencyList, ringListReduced));
				moleculeSet.m_moleculeTypes.back().m_mass = mass;
				
				moleculeSet.m_molecules.emplace_back(MoleculeSet::Molecule(&moleculeSet, countMol, countMolType, atomList));
				
				++countMolType;
			} else {
				moleculeSet.m_molecules.emplace_back(MoleculeSet::Molecule(&moleculeSet, countMol, it2->getId(), atomList));
			}
			++countMol;
		}
	}
	
	// Sort molecules by molecular mass
	std::vector<std::size_t> perm(moleculeSet.m_moleculeTypes.size());
	std::iota(perm.begin(), perm.end(), 0);
	std::sort(perm.begin(), perm.end(), [&moleculeSet] (std::size_t a, std::size_t b) { return moleculeSet.m_moleculeTypes[a].getMass() > moleculeSet.m_moleculeTypes[b].getMass(); });
	
	for (std::size_t i = 0; i < moleculeSet.m_moleculeTypes.size(); ++i) {
		moleculeSet.m_moleculeTypes[i].m_id = perm[i];
	}
	for (std::size_t i = 0; i < moleculeSet.m_molecules.size(); ++i) {
		moleculeSet.m_molecules[i].m_moleculeTypeId = perm[moleculeSet.m_molecules[i].getMoleculeType().getId()];
	}
	std::sort(moleculeSet.m_moleculeTypes.begin(), moleculeSet.m_moleculeTypes.end(), [] (const MoleculeSet::MoleculeType &a, const MoleculeSet::MoleculeType &b) { return a.getId() < b.getId(); });
	std::sort(moleculeSet.m_molecules.begin(), moleculeSet.m_molecules.end(), [] (const MoleculeSet::Molecule &a, const MoleculeSet::Molecule &b) { return a.getMoleculeType().getId() < b.getMoleculeType().getId(); });
	
	for (std::size_t i = 0; i < moleculeSet.m_molecules.size(); ++i) {
		moleculeSet.m_molecules[i].m_id = i;
	}
}

void RecognizeMolecules::p_recognizeMoleculesRecursion(std::size_t beginAtom, std::size_t lastAtom, std::size_t beginAtomIndex) {
	for (std::size_t i = 0; i < m_recursionInfo->atomMolIndex.size(); ++i) {
		if (i == beginAtom || i == lastAtom)
			continue;
		if (m_recursionInfo->atomMolIndex[i] == 0) {
			if (p_bondRange(m_recursionInfo->atomPositions.at(beginAtom), m_recursionInfo->atomPositions.at(i), m_recursionInfo->atoms.at(beginAtom).getAtomType().getCovalentRadius(), m_recursionInfo->atoms.at(i).getAtomType().getCovalentRadius())) {
				m_recursionInfo->atomMolIndex[i] = m_recursionInfo->atomMolIndex[beginAtom];
				std::size_t index = m_recursionInfo->helperMolecule.size();
				m_recursionInfo->helperMolecule.push_back(i);
				m_recursionInfo->helperMoleculeBonds.emplace_back(beginAtomIndex, index);
				p_recognizeMoleculesRecursion(i, beginAtom, index);
			}
		} else if (m_recursionInfo->atomMolIndex[i] == m_recursionInfo->atomMolIndex[beginAtom]) {
			if (p_bondRange(m_recursionInfo->atomPositions.at(beginAtom), m_recursionInfo->atomPositions.at(i), m_recursionInfo->atoms.at(beginAtom).getAtomType().getCovalentRadius(), m_recursionInfo->atoms.at(i).getAtomType().getCovalentRadius())) {
				auto it = std::find(m_recursionInfo->helperMolecule.cbegin(), m_recursionInfo->helperMolecule.cend(), i);
				std::size_t index = static_cast<std::size_t>(std::distance(m_recursionInfo->helperMolecule.cbegin(), it));
				if (index > beginAtomIndex)
					m_recursionInfo->helperMoleculeBondsCross.emplace_back(beginAtomIndex, index);
			}
		}
	}
}

bool RecognizeMolecules::p_bondRange(const Vector3 &position1, const Vector3 &position2, const double radius1, const double radius2) {
	if (radius1 < m_radiusFactor || radius2 < m_radiusFactor)
		return false;
	
	Vector3 dist = position1 - position2;
	dist = m_cellInfo.foldVector(dist);
	
	const double distance = dist.norm();
	return distance < (radius1 + radius2) * m_bondStretchingFactor;
}

void RecognizeMolecules::p_findPathRecursion(std::size_t current) {
	for (std::size_t a: m_recursionInfoPath->adjacencyList[current]) {
		if (a == m_recursionInfoPath->target) {
			m_recursionInfoPath->stack.push_back(a);
			m_recursionInfoPath->ringList.push_back(m_recursionInfoPath->stack);
			m_recursionInfoPath->stack.pop_back();
		} else if (m_recursionInfoPath->unused[a]) {
			m_recursionInfoPath->stack.push_back(a);
			m_recursionInfoPath->unused[a] = false;
			p_findPathRecursion(a);
			m_recursionInfoPath->stack.pop_back();
			m_recursionInfoPath->unused[a] = true;
		}
	}
}

}

}

#endif


