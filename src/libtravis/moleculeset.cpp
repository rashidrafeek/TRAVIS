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


#include "moleculeset.h"

#include "vector.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace Libtravis {

namespace Travis {

namespace {
const std::string validAtomTypeCharacters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_#";
}

std::vector<std::size_t> MoleculeSet::AtomType::findAtomIds() const {
	std::vector<std::size_t> list;
	std::for_each(m_moleculeSet->m_atoms.cbegin(), m_moleculeSet->m_atoms.cend(), [&list, this] (const Atom &atom) { if (atom.getAtomType() == *this) list.push_back(atom.getId()); });
	return list;
}

std::vector<std::reference_wrapper<const MoleculeSet::Atom>> MoleculeSet::AtomType::findAtoms() const {
	std::vector<std::reference_wrapper<const Atom>> list;
	std::for_each(m_moleculeSet->m_atoms.cbegin(), m_moleculeSet->m_atoms.cend(), [&list, this] (const Atom &atom) { if (atom.getAtomType() == *this) list.push_back(atom); });
	return list;
}

MoleculeSet::AtomType::AtomType(const MoleculeSet *moleculeSet, std::size_t id, const std::string &name) : m_moleculeSet(moleculeSet), m_id(id), m_name(name) {
	if (m_name.find_first_not_of(validAtomTypeCharacters) != m_name.npos)
		throw std::runtime_error("AtomType: Invalid characters in name");
}

std::size_t MoleculeSet::MoleculeType::findLocalAtom(const std::string &localAtomName) const {
	std::size_t pos = localAtomName.find_first_not_of(' ');
	std::size_t end = localAtomName.find_first_not_of(validAtomTypeCharacters, pos);
	if (end == pos)
		throw std::runtime_error("MoleculeType::findLocalAtom: Expected atom type at beginning of \"" + localAtomName.substr(pos) + "\"");
	const AtomType &type = m_moleculeSet->findAtomType(localAtomName.substr(pos, end - pos));
	pos = end;
	
	pos = localAtomName.find_first_not_of(' ', pos);
	end = localAtomName.find_first_not_of("0123456789", pos);
	if (end == pos)
		throw std::runtime_error("MoleculeType::findLocalAtom: Expected number at beginning of \"" + localAtomName.substr(pos) + "\"");
	
	unsigned long num = std::stoul(localAtomName.substr(pos, end - pos));
	if (num == 0 || num > getAtomTypeCount()[type.getId()])
		throw std::runtime_error("MoleculeType::findLocalAtom: Invalid number");
	--num;
	
	std::size_t localAtomId = 0;
	for (std::size_t i = 0; i < type.getId(); ++i)
		localAtomId += getAtomTypeCount()[i];
	localAtomId += num;
	
	return localAtomId;
}

std::vector<std::size_t> MoleculeSet::MoleculeType::findLocalAtoms(const std::string &atomSelection) const {
	std::vector<std::size_t> list;
	
	if (atomSelection == "*") {
		list.resize(m_atomCodes.size());
		std::iota(list.begin(), list.end(), 0);
		return list;
	}
	
	const std::string typeChars = validAtomTypeCharacters;
	const std::string digits = "0123456789";
	
	std::size_t pos = 0;
	std::size_t type = 0;
	std::size_t baseId = 0;
	unsigned long num1 = 0, num2 = 0;
	
	enum class State { Initial, Type, Number, CommaAfterType, CommaAfterNumber, Range };
	State state = State::Initial;
	
	for(;;) {
		pos = atomSelection.find_first_not_of(' ', pos);
		if (pos >= atomSelection.size()) {
			if (state == State::CommaAfterType || state == State::CommaAfterNumber)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected number or atom type, but string ended");
			if (state == State::Range)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected number, but string ended");
			if (state == State::Type) {
				for (std::size_t i = baseId; i < baseId + m_atomTypeCount[type]; ++i)
					list.push_back(i);
			} else if (state == State::Number) {
				for (std::size_t i = baseId + num1; i <= baseId + num2; ++i)
					list.push_back(i);
			}
			break;
		}
		
		if (typeChars.find(atomSelection[pos]) != typeChars.npos) {
			if (state == State::Type)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected number or \",\", but found atom type at beginning of \"" + atomSelection.substr(pos) + "\"");
			if (state == State::Number)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected \"-\" or \",\", but found atom type at beginning of \"" + atomSelection.substr(pos) + "\"");
			if (state == State::Range)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected number, but found atom type at beginning of \"" + atomSelection.substr(pos) + "\"");
			if (state == State::CommaAfterType) {
				for (std::size_t i = baseId; i < baseId + m_atomTypeCount[type]; ++i)
					list.push_back(i);
			} else if (state == State::CommaAfterNumber) {
				for (std::size_t i = baseId + num1; i <= baseId + num2; ++i)
					list.push_back(i);
			}
			std::size_t end = atomSelection.find_first_not_of(typeChars, pos);
			type = m_moleculeSet->findAtomType(atomSelection.substr(pos, end - pos)).getId();
			baseId = 0;
			for (std::size_t i = 0; i < type; ++i)
				baseId += m_atomTypeCount[i];
			state = State::Type;
			pos = end;
		} else if (digits.find(atomSelection[pos]) != digits.npos) {
			if (state == State::Initial)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected atom type at beginning of \"" + atomSelection.substr(pos) + "\"");
			if (state == State::Number)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected \"-\" or \",\", but found number at beginning of \"" + atomSelection.substr(pos) + "\"");
			if (state == State::CommaAfterType)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected atom type, but found number at beginning of \"" + atomSelection.substr(pos) + "\"");
			if (state == State::CommaAfterNumber) {
				for (std::size_t i = baseId + num1; i <= baseId + num2; ++i)
					list.push_back(i);
			}
			std::size_t end = atomSelection.find_first_not_of(digits, pos);
			unsigned long num = std::stoul(atomSelection.substr(pos, end - pos));
			if (num == 0 || num > m_atomTypeCount[type])
				throw std::runtime_error("MoleculeType::findLocalAtoms: Invalid number");
			--num;
			if (state == State::Range) {
				num1 = num2;
				num2 = num;
			} else {
				num1 = num;
				num2 = num;
			}
			state = State::Number;
			pos = end;
		} else if (atomSelection[pos] == ',') {
			if (state == State::Initial)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected atom type at beginning of \"" + atomSelection.substr(pos) + "\"");
			if (state == State::CommaAfterType || state == State::CommaAfterNumber)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected number or atom type, but found \",\" at beginning of \"" + atomSelection.substr(pos) + "\"");
			if (state == State::Range)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected number, but found \",\" at beginning of \"" + atomSelection.substr(pos) + "\"");
			++pos;
			if (state == State::Type)
				state = State::CommaAfterType;
			else
				state = State::CommaAfterNumber;
		} else if (atomSelection[pos] == '-') {
			if (state == State::Initial)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected atom type at beginning of \"" + atomSelection.substr(pos) + "\"");
			if (state == State::CommaAfterType || state == State::CommaAfterNumber)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected number or atom type, but found \"-\" at beginning of \"" + atomSelection.substr(pos) + "\"");
			if (state == State::Type)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected number or \",\", but found \"-\" at beginning of \"" + atomSelection.substr(pos) + "\"");
			if (state == State::Range)
				throw std::runtime_error("MoleculeType::findLocalAtoms: Expected number, but found \"-\" at beginning of \"" + atomSelection.substr(pos) + "\"");
			++pos;
			state = State::Range;
		} else {
			throw std::runtime_error("MoleculeType::findLocalAtoms: Unknown character found at beginning of \"" + atomSelection.substr(pos) + "\"");
		}
	}
	
	std::sort(list.begin(), list.end());
	auto end = std::unique(list.begin(), list.end());
	list.erase(end, list.end());
	
	return list;
}

std::string MoleculeSet::MoleculeType::findLocalAtomName(std::size_t localAtomId) const {
	std::size_t type = 0, countId = 0;
	while (type < m_atomTypeCount.size() && localAtomId >= countId + m_atomTypeCount[type]) {
		countId += m_atomTypeCount[type];
		++type;
	}
	if (type >= m_atomTypeCount.size())
		throw std::runtime_error("MoleculeType::findLocalAtomName: Molecule type \"" + m_name + "\" has less than " + std::to_string(localAtomId) + " atoms");
	return m_moleculeSet->m_atomTypes[type].getName() + std::to_string(localAtomId - countId + 1);
}

const MoleculeSet::AtomType &MoleculeSet::MoleculeType::findLocalAtomType(std::size_t localAtomId) const {
	std::size_t type = 0, countId = 0;
	while (type < m_atomTypeCount.size() && localAtomId >= countId + m_atomTypeCount[type]) {
		countId += m_atomTypeCount[type];
		++type;
	}
	if (type >= m_atomTypeCount.size())
		throw std::runtime_error("MoleculeType::findLocalAtomName: Molecule type \"" + m_name + "\" has less than " + std::to_string(localAtomId) + " atoms");
	return m_moleculeSet->m_atomTypes[type];
}

std::vector<std::reference_wrapper<const MoleculeSet::Molecule>> MoleculeSet::MoleculeType::findMolecules() const {
	std::vector<std::reference_wrapper<const Molecule>> list;
	std::for_each(m_moleculeSet->m_molecules.cbegin(), m_moleculeSet->m_molecules.cend(), [&list, this] (const Molecule &molecule) { if (molecule.getMoleculeType() == *this) list.push_back(molecule); });
	return list;
}

std::vector<std::size_t> MoleculeSet::Molecule::findAtomIds(const std::vector<std::size_t> &localAtomIds) const {
	std::vector<std::size_t> list(localAtomIds.size());
	for (std::size_t i = 0; i < localAtomIds.size(); ++i)
		list[i] = m_atoms.at(localAtomIds[i]);
	return list;
}

std::vector<std::reference_wrapper<const MoleculeSet::Atom>> MoleculeSet::Molecule::findAtoms(const std::vector<std::size_t> &localAtomIds) const {
	std::vector<std::reference_wrapper<const Atom>> list;
	list.reserve(localAtomIds.size());
	for (auto &&i: localAtomIds)
		list.push_back(m_moleculeSet->m_atoms[m_atoms.at(i)]);
	return list;
}

std::vector<double> MoleculeSet::Molecule::getAtomMasses() const {
	std::vector<double> masses;
	masses.reserve(m_atoms.size());
	for (auto &&a: m_atoms)
		masses.push_back(m_moleculeSet->m_atoms[a].getAtomType().getMass());
	return masses;
}

const MoleculeSet::AtomType &MoleculeSet::findAtomType(const std::string &name) const {
	auto it = std::find_if(m_atomTypes.cbegin(), m_atomTypes.cend(), [&name] (const AtomType &atomType) { return atomType.getName() == name; });
	if (it == m_atomTypes.cend())
		throw std::runtime_error("MoleculeSet: Atom type \"" + name + "\" not known");
	return *it;
}

const MoleculeSet::MoleculeType &MoleculeSet::findMoleculeType(const std::string &name) const {
	auto it = std::find_if(m_moleculeTypes.cbegin(), m_moleculeTypes.cend(), [&name] (const MoleculeType &moleculeType) { return moleculeType.getName() == name; });
	if (it == m_moleculeTypes.cend())
		throw std::runtime_error("MoleculeSet: Molecule type \"" + name + "\" not known");
	return *it;
}

}

}

#endif


