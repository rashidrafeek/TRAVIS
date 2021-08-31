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


#ifndef LT_MOLECULESET_H
#define LT_MOLECULESET_H


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include "element.h"
#include "vector.h"

#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Libtravis {

namespace Travis {

/*!
 * \brief Results of molecule recognition
 * 
 * This class contains the results of the molecule recognition by RecognizeMolecules. If an object of this class is created by RecognizeMolecules, the described properties of the components are fulfilled.
 */
class MoleculeSet
{
	friend class RecognizeMolecules;
	
public:
	class Atom;
	/*!
	 * \brief An atom type
	 * 
	 * This class represents an atom type in a MoleculeSet. The IDs of atom types are assigned according to the alphabetical ordering of their names.
	 */ 
	class AtomType
	{
		friend class MoleculeSet;
		friend class RecognizeMolecules;
		
	public:
		//! Get the ID
		std::size_t getId() const { return m_id; }
		//! Get the name
		const std::string &getName() const { return m_name; }
		//! Get the mass in u
		double getMass() const { return m_mass; }
		//! Get the covalent radius in pm
		double getCovalentRadius() const { return m_covalentRadius; }
		//! Get the van der Waals radius in pm
		double getVdWRadius() const { return m_vdWRadius; }
		
		//! Find the IDs of all atoms of this atom type
		std::vector<std::size_t> findAtomIds() const;
		//! Find all atoms of this atom type
		std::vector<std::reference_wrapper<const Atom>> findAtoms() const;
		
		//! Returns true, if the atom types have the same ID
		bool operator==(const AtomType &other) const { return m_id == other.m_id; }
		//! Returns true, if the atom types do not have the same ID
		bool operator!=(const AtomType &other) const { return m_id != other.m_id; }
		
	private:
		const MoleculeSet *m_moleculeSet;
		std::size_t m_id;
		std::string m_name;
		double m_mass = 0.0;
		double m_covalentRadius = 0.0;
		double m_vdWRadius = 0.0;
		
		//! Construct an atom type with ID and name
		AtomType(const MoleculeSet *moleculeSet, std::size_t id, const std::string &name);
		//! Construct and atom type with ID, name, and an element to assign mass, covalent radius, and van der Waals radius
		AtomType(const MoleculeSet *moleculeSet, std::size_t id, const std::string &name, const Element &element, double radiusFactor) : AtomType(moleculeSet, id, name) {
			m_mass = element.mass;
			m_covalentRadius = element.covalentRadius * radiusFactor;
			m_vdWRadius = element.vdWRadius * radiusFactor;
		}
	};
	
	/*!
	 * \brief An atom
	 * 
	 * This class represents an atom in a MoleculeSet. The IDs are assigned according to the indices in the trajectory.
	 */
	class Atom
	{
		friend class MoleculeSet;
		friend class RecognizeMolecules;
		
	public:
		//! Get the ID
		std::size_t getId() const { return m_id; }
		//! Get the ID of the atom type
		std::size_t getAtomTypeId() const { return m_atomTypeId; }
		//! Get the atom type
		const AtomType &getAtomType() const { return m_moleculeSet->m_atomTypes[m_atomTypeId]; }
		
	private:
		const MoleculeSet *m_moleculeSet;
		std::size_t m_id;
		std::size_t m_atomTypeId;
		
		//! Construct an atom with ID and ID of atom type
		Atom(const MoleculeSet *moleculeSet, std::size_t id, std::size_t atomTypeId) : m_moleculeSet(moleculeSet), m_id(id), m_atomTypeId(atomTypeId) {}
	};
	
	class Molecule;
	/*!
	 * \brief A molecule type
	 * 
	 * This class represents a molecule type in a MoleculeSet. The IDs are assigned according to the ordering by molecular mass.
	 */
	class MoleculeType
	{
		friend class MoleculeSet;
		friend class RecognizeMolecules;
		
	public:
		//! Get the ID
		std::size_t getId() const { return m_id; }
		//! Get the name
		const std::string &getName() const { return m_name; }
		/*!
		 * \brief Get the number of atoms of different atom types
		 * 
		 * The returned array has one entry for each atom type in the MoleculeSet giving the number of atoms of this type in the molecule type.
		 */
		const std::vector<std::size_t> &getAtomTypeCount() const { return m_atomTypeCount; }
		/*!
		 * \brief Get the atom codes
		 * 
		 * The returned array has one entry for each atom in the molecule type. The atoms are orderd by their local ID, meaning that they are sorted by the respective atom types, and atoms of the same type are sorted by their atom code.
		 */
		const std::vector<double> &getAtomCodes() const { return m_atomCodes; }
		/*!
		 * \brief Get the bonds
		 * 
		 * The returned array contains one array for each atom of the molecule type. The atoms are sorted in the same way as in the atom codes. The inner arrays contain all the neighbors of the respective atom. The neighbors are given by the local ID.
		 */
		const std::vector<std::vector<std::size_t>> &getBonds() const { return m_bonds; }
		/*!
		 * \brief Get the rings
		 * 
		 * The returned array contains one array for each ring of the molecule type. The inner arrays contain the local IDs of the ring members.
		 */
		const std::vector<std::vector<std::size_t>> &getRings() const { return m_rings; }
		
		//! Get the molecular mass in u
		double getMass() const { return m_mass; }
		
		//! Get the number of atoms
		std::size_t getNumAtoms() const { return m_atomCodes.size(); }
		
		//! Find the local ID of an atom by its name
		std::size_t findLocalAtom(const std::string &localAtomName) const;
		//! Find the local IDs of an atom selection
		std::vector<std::size_t> findLocalAtoms(const std::string &atomSelection) const;
		
		//! Find the name of an atom by its local ID
		std::string findLocalAtomName(std::size_t localAtomId) const;
		//! Find the atom type of an atom by its local ID
		const AtomType &findLocalAtomType(std::size_t localAtomId) const;
		
		//! Find all molecules of this molecule type
		std::vector<std::reference_wrapper<const Molecule>> findMolecules() const;
		
		//! Returns true, if the molecule types have the same ID
		bool operator==(const MoleculeType &other) const { return m_id == other.m_id; }
		//! Returns true, if the molecule types do not have the same ID
		bool operator!=(const MoleculeType &other) const { return m_id == other.m_id; }
		
	private:
		const MoleculeSet *m_moleculeSet;
		std::size_t m_id;
		std::string m_name;
		std::vector<std::size_t> m_atomTypeCount;
		std::vector<double> m_atomCodes;
		std::vector<std::vector<std::size_t>> m_bonds;
		std::vector<std::vector<std::size_t>> m_rings;
		double m_mass = 0.0;
		
		//! Construct a molecule type with ID, name, number of atoms of different atom types, atom codes, bonds, and rings
		MoleculeType(const MoleculeSet *moleculeSet, std::size_t id, const std::string &name, const std::vector<std::size_t> &atomTypeCount, const std::vector<double> &atomCodes, const std::vector<std::vector<std::size_t>> &bonds, const std::vector<std::vector<std::size_t>> &rings) : m_moleculeSet(moleculeSet), m_id(id), m_name(name), m_atomTypeCount(atomTypeCount), m_atomCodes(atomCodes), m_bonds(bonds), m_rings(rings) {}
	};
	
	/*!
	 * \brief A molecule
	 * 
	 * This class represents a molecule in a MoleculeSet. The IDs are assigned according to the molecule type. The order of molecules of the same type is random.
	 */
	class Molecule
	{
		friend class MoleculeSet;
		friend class RecognizeMolecules;
		
	public:
		//! Get the ID
		std::size_t getId() const { return m_id; }
		//! Get the ID of the molecule type
		std::size_t getMoleculeTypeId() const { return m_moleculeTypeId; }
		//! Get the molecule type
		const MoleculeType &getMoleculeType() const { return m_moleculeSet->m_moleculeTypes[m_moleculeTypeId]; }
		/*!
		 * \brief Get the atom IDs
		 * 
		 * The returned array containes the ID of each Atom in the molecule. The ordering is according to the local ID in the molecule type, so this array translates local IDs from the molecule type to IDs of atoms.
		 */
		const std::vector<std::size_t> &getAtomIds() const { return m_atoms; }
		
		//! Get the number of atoms
		std::size_t getNumAtoms() const { return m_atoms.size(); }
		
		//! Find the ID of an atom by its local ID
		std::size_t findAtomId(std::size_t localAtomId) const { return m_atoms.at(localAtomId); }
		//! Find the ID of an atom by its name
		std::size_t findAtomId(const std::string &localAtomName) const { return m_atoms[getMoleculeType().findLocalAtom(localAtomName)]; }
		//! Find an atom by its local ID
		const Atom &findAtom(std::size_t localAtomId) const { return m_moleculeSet->m_atoms[findAtomId(localAtomId)]; }
		//! Find an atom by its name
		const Atom &findAtom(const std::string &localAtomName) const { return m_moleculeSet->m_atoms[findAtomId(localAtomName)]; }
		
		//! Find the IDs of atoms by their local IDs
		std::vector<std::size_t> findAtomIds(const std::vector<std::size_t> &localAtomIds) const;
		//! Find the IDs of atoms by an atom selection
		std::vector<std::size_t> findAtomIds(const std::string &atomSelection) const { return findAtomIds(m_moleculeSet->m_moleculeTypes[m_moleculeTypeId].findLocalAtoms(atomSelection)); }
		//! Find atoms by their local IDs
		std::vector<std::reference_wrapper<const Atom>> findAtoms(const std::vector<std::size_t> &localAtomIds) const;
		//! Find atoms by an atom selection
		std::vector<std::reference_wrapper<const Atom>> findAtoms(const std::string &atomSelection) const { return findAtoms(m_moleculeSet->m_moleculeTypes[m_moleculeTypeId].findLocalAtoms(atomSelection)); }
		
		//! Get masses of all atoms
		std::vector<double> getAtomMasses() const;
		
	private:
		const MoleculeSet *m_moleculeSet;
		std::size_t m_id;
		std::size_t m_moleculeTypeId;
		std::vector<std::size_t> m_atoms;
		
		//! Construct a molecule with ID, ID of molecule type, and list of atom IDs
		Molecule(const MoleculeSet *moleculeSet, std::size_t id, std::size_t moleculeTypeId, const std::vector<std::size_t> &atoms) : m_moleculeSet(moleculeSet), m_id(id), m_moleculeTypeId(moleculeTypeId), m_atoms(atoms) {}
	};
	
	MoleculeSet() = default;
	MoleculeSet(const MoleculeSet &other) : m_atomTypes(other.m_atomTypes), m_atoms(other.m_atoms), m_moleculeTypes(other.m_moleculeTypes), m_molecules(other.m_molecules) {
		p_updatePointers();
	}
	MoleculeSet(MoleculeSet &&other) noexcept {
		swap(*this, other);
	}
	MoleculeSet &operator=(MoleculeSet other) noexcept {
		swap(*this, other);
		return *this;
	}
	friend void swap(MoleculeSet &first, MoleculeSet &second) noexcept {
		using std::swap;
		swap(first.m_atoms, second.m_atoms);
		swap(first.m_atomTypes, second.m_atomTypes);
		swap(first.m_molecules, second.m_molecules);
		swap(first.m_moleculeTypes, second.m_moleculeTypes);
		first.p_updatePointers();
		second.p_updatePointers();
	}
	
	/*!
	 * \brief Get the array of atom types
	 * 
	 * The atom types are sorted by their ID.
	 */
	const std::vector<AtomType> &getAtomTypes() const { return m_atomTypes; }
	/*!
	 * \brief Get the array of atoms
	 * 
	 * The atoms are sorted by their ID.
	 */
	const std::vector<Atom> &getAtoms() const { return m_atoms; }
	/*! 
	 * \brief Get the array of molecule types
	 * 
	 * The molecule types are sorted by their ID
	 */
	const std::vector<MoleculeType> &getMoleculeTypes() const { return m_moleculeTypes; }
	/*! 
	 * \brief Get the array of molecules
	 * 
	 * The molecules are sorted by their ID
	 */
	const std::vector<Molecule> &getMolecules() const { return m_molecules; }
	
	//! Get a specific atom type
	const AtomType &getAtomType(std::size_t id) const { return m_atomTypes.at(id); }
	//! Get a specific atom
	const Atom &getAtom(std::size_t id) const { return m_atoms.at(id); }
	//! Get a specific molecule type
	const MoleculeType &getMoleculeType(std::size_t id) const { return m_moleculeTypes.at(id); }
	//! Get a specific molecule
	const Molecule &getMolecule(std::size_t id) const { return m_molecules.at(id); }
	
	//! Get the number of atom types
	std::size_t getNumAtomTypes() const { return m_atomTypes.size(); }
	//! Get the number of atoms
	std::size_t getNumAtoms() const { return m_atoms.size(); }
	//! Get the number of molecule types
	std::size_t getNumMoleculeTypes() const { return m_moleculeTypes.size(); }
	//! Get the number of molecules
	std::size_t getNumMolecules() const { return m_molecules.size(); }
	
	//! Find the ID of an atom type by its name
	const AtomType &findAtomType(const std::string &name) const;
	//! Find the ID of a molecule type by its name
	const MoleculeType &findMoleculeType(const std::string &name) const;
	
	//! Find IDs of all atoms of the atom type
	std::vector<std::size_t> findAtomIds(const AtomType &atomType) const { return atomType.findAtomIds(); }
	//! Find IDs of all atoms by the ID of the atom type
	std::vector<std::size_t> findAtomIds(std::size_t atomTypeId) const { return m_atomTypes.at(atomTypeId).findAtomIds(); }
	//! Find IDs of all atoms by the name of the atom type
	std::vector<std::size_t> findAtomIds(const std::string &name) const { return findAtomType(name).findAtomIds(); }
	//! Find all atoms of the atom type
	std::vector<std::reference_wrapper<const Atom>> findAtoms(const AtomType &atomType) const { return atomType.findAtoms(); }
	//! Find all atoms by the ID of the atom type
	std::vector<std::reference_wrapper<const Atom>> findAtoms(std::size_t atomTypeId) const { return m_atomTypes.at(atomTypeId).findAtoms(); }
	//! Find all atoms by the name of the atom type
	std::vector<std::reference_wrapper<const Atom>> findAtoms(const std::string &name) const { return findAtomType(name).findAtoms(); }
	//! Find all molecules of the molecule type
	std::vector<std::reference_wrapper<const Molecule>> findMolecules(const MoleculeType &moleculeType) const { return moleculeType.findMolecules(); }
	//! Find all molecules by the ID of the molecule type
	std::vector<std::reference_wrapper<const Molecule>> findMolecules(std::size_t moleculeTypeId) const { return m_moleculeTypes.at(moleculeTypeId).findMolecules(); }
	//! Find all molecules by the name of the molecule type
	std::vector<std::reference_wrapper<const Molecule>> findMolecules(const std::string &name) const { return findMoleculeType(name).findMolecules(); }
	
private:
	std::vector<AtomType> m_atomTypes;
	std::vector<Atom> m_atoms;
	std::vector<MoleculeType> m_moleculeTypes;
	std::vector<Molecule> m_molecules;
	
	void p_updatePointers() noexcept {
		for (auto &&at: m_atomTypes)
			at.m_moleculeSet = this;
		for (auto &&a: m_atoms)
			a.m_moleculeSet = this;
		for (auto &&mt: m_moleculeTypes)
			mt.m_moleculeSet = this;
		for (auto &&m: m_molecules)
			m.m_moleculeSet = this;
	}
};

}

}

#endif

#endif


