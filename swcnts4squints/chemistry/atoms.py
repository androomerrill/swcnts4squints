# -*- coding: utf-8 -*-
"""
===================================================================
Abstract representation of Atoms (:mod:`tasr.tools.chemistry.atoms`)
===================================================================

.. currentmodule:: tasr.tools.chemistry.atoms

.. autosummary::
   :toctree: generated/

   Atoms

"""
from __future__ import division, absolute_import, print_function
__docformat__ = 'restructuredtext'

import copy
import math

from collections import MutableSequence

import numpy as np

from .atom import Atom


class Atoms(MutableSequence):

    """Class for creating abstract represention of a collection of Atoms.

    .. versionadded:: 0.3.14

    Parameters
    ----------
    atoms : {None, sequence}, optional
        if not ``None``, then a list of ``Atom`` objects
    copylist : bool, optional
        perform shallow copy of atoms list
    deepcopy : bool, optional
        perform deepcopy of atoms list

    """

    def __init__(self, atoms=None, copylist=True, deepcopy=False):
        self._atoms = []
        self._masses = []
        self._charges = []
        self._coords = []

        self._atomtypes = {}

        self._property_lists = \
            {'m': self._masses, 'q': self._charges, 'r': self._coords}
        if atoms is not None:
            if isinstance(atoms, list):
                if copylist and not deepcopy:
                    self._atoms = atoms[:]
                elif deepcopy:
                    self._atoms = copy.deepcopy(atoms)
                else:
                    self._atoms = atoms

                for atom in self._atoms:
                    self._check_type(atom)
                    for p, plist in self._property_lists.iteritems():
                        plist.append(getattr(atom, p))
            #elif isinstance(atoms, Atoms):
            #    if copylist and not deepcopy:
            #        self._atoms.extend(atoms.atomlist[:])
            #    elif deepcopy:
            #        self._atoms.extend(copy.deepcopy(atoms.atomlist))
            #    else:
            #        self._atoms.extend(atoms.atomlist)
            else:
                raise TypeError('Atoms(atoms={!r}) '.format(atoms) +
                                'is not a valid Atoms constructor argument.\n'
                                'atoms must be None, a list of Atom object '
                                'instances, or an Atoms object\ninstance.')

    #@property
    #def atoms(self):
    #    return self._atoms

    @property
    def atomlist(self):
        """Return the list of Atom objects"""
        return self._atoms

    @property
    def atomtypes(self):
        """Return the atom type dictionary."""
        for atom in self._atoms:
            if atom.atomtype not in self._atomtypes:
                self._atomtypes[atom.atomtype] = {}
                self._atomtypes[atom.atomtype]['mass'] = atom.m
                self._atomtypes[atom.atomtype]['q'] = atom.q
        return self._atomtypes

    def add_atomtype(self, atom):
        """Add atom type to atom type dictionary.

        Parameters
        ----------
        atom : Atom
            an intance of an Atom object

        """
        if atom.atomtype not in self._atomtypes:
            self._atomtypes[atom.atomtype] = {}
            self._atomtypes[atom.atomtype]['mass'] = atom.m
            self._atomtypes[atom.atomtype]['q'] = atom.q

    def add_atomtypes(self, atomtypes=[]):
        """Add atom type in atomtypes list to atomtype dictionary.
        Parameters
        ----------
        atomtypes : sequence
            a list of Atom object instances

        """
        for atom in atomtypes:
            self.add_atomtype(atom)

    @property
    def Natoms(self):
        """Number of atoms in Atoms."""
        return len(self._atoms)

    @property
    def Ntypes(self):
        """Number of atom types in Atoms."""
        return len(self.atomtypes.keys())

    @property
    def masses(self):
        """Return the list of Atom masses"""
        #self._masses_array = np.asarray(self._masses)
        return self._masses

    @property
    def m(self):
        """Return the total mass of all atoms"""
        return math.fsum(self._masses)

    @property
    def charges(self):
        """Return list of Atom charges"""
        #self._charges_array = np.asarray(self._charges)
        return np.asarray(self._charges)

    @property
    def q(self):
        """Return the total net charge of all atoms."""
        return np.asarray(self._charges).sum()

    @property
    def coords(self):
        """Return array of Atom position array"""
        #self._coords_array = np.asarray(self._coords)
        return np.asarray(self._coords)

    @property
    def r(self):
        """Return the center of mass cooridates of all atoms."""
        return self.CM()

    #@property
    #def properties(self):
    #    return self._property_lists

    def _check_type(self, value):
        """Check that value is instance of ``Atom`` class.

        Parameters
        ----------
        value : Atom
            value to type check

        Raises
        ------
        TypeError
            if ``value`` is not instance of ``Atom``

        """
        if not isinstance(value, Atom):
            raise TypeError('{} is not an Atom.'.format(value))

    def __len__(self):
        """Concrete implementation of @abstractmethod.

        Returns
        -------
        int
            length of ``self.atoms`` list.

        """
        return len(self._atoms)

    def __delitem__(self, index):
        """Concrete implementation of @abstractmethod.

        delete list element ``self.atoms[index]`` and delete all elements
        from atom properties lists ``self.masses[index]``,
        ``self.charges[index]``, and ``self.coords[index]``

        Parameters
        ----------
        index : int
            index of target list element

        """
        del self._atoms[index]
        for plist in self._property_lists.itervalues():
            del plist[index]

    def __getitem__(self, index):
        """Concrete implementation of @abstractmethod.

        get Atom object instance at list element ``self.atoms[index]``

        Parameters
        ----------
        index : int
            index of target list element

        Returns
        -------
        Atom
            Atom object instance at target ``self.atoms[index]``

        """
        return self._atoms[index]

    def __setitem__(self, index, atom):
        """Concrete implementation of @abstractmethod.

        set target list element ``self.atoms[index] = atom``

        Also set element of all atom properties lists (``self.masses[index]``,
        ``self.charges[index]``, and ``self.coords[index]``) to atom instance
        properties (``atom.m``, ``atom.q``, ``atom.r``), respectively.


        Parameters
        ----------
        index : int
            index of target list element
        atom : Atom
            Atom object instance to set target list element to

        """
        self._check_type(atom)
        self._atoms[index] = atom
        for p, plist in self._property_lists.iteritems():
            plist[index] = getattr(atom, p)

    def insert(self, index, atom):
        """Concrete implementation of @abstractmethod.

        insert Atom instance at target list ``index``

        Also insert Atom instance properties at the given target list index
        for all Atom properties in ``self._property_lists.keys()``
        into their respective target lists of Atoms properties
        ``self._property_lists.values()``.


        Parameters
        ----------
        index : int
            index of target list element
        atom : Atom
            Atom object instance to set target list element to

        """
        self._check_type(atom)
        self._atoms.insert(index, atom)
        for p, plist in self._property_lists.iteritems():
            plist.insert(index, getattr(atom, p))

    def CM(self):
        """Compute the center-of-mass coordinates of Atoms.

        Returns
        -------
        ndarray
            3-element ndarray specifying center-of-mass coordinates of Atoms.

        """
        mass = np.asarray([self._masses])
        coords = np.asarray(self._coords)
        mr = mass.T * coords
        return np.sum(mr, axis=0) / np.sum(mass)

    def center_CM(self, r_indices=[0, 1, 2]):
        """Center atoms on CM coordinates

        Parameters
        ----------
        r_indices : sequence, optional
            list of component indices of CM position vector to translate
            and center on the origin.

        """
        dr = -self.CM()
        self.translate(dr, r_indices=r_indices)

    def fix_minus_zero_coords(self, epsilon=1.0e-10):
        """Set really really small negative coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon zero so we don't end up with -0.00000
        coordinates in structure data output.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any component of an
            atom's position

        """
        for atom in self._atoms:
            r = atom.r.tolist()
            for i, ri in enumerate(r[:]):
                if abs(ri) < epsilon:
                    r[i] = 0.0
            atom.x, atom.y, atom.z = r

    def rotate(self, R_matrix):
        """Rotate atom coordinates using rotation matrix ``R_matrix``.

        Parameters
        ----------
        R_matrix : array_like
            3x3 array representation of ``3D`` rotation matrix.

        """
        for atom in self._atoms:
            atom.r = np.dot(R_matrix, atom.r.T).T

    def translate(self, dr, r_indices=[0, 1, 2]):
        """Translate atom coordinates.

        Parameters
        ----------
        dr : array_like
            array representation of displacement vector to translate
            components of atoms position vector through
        r_indices : sequence, optional
            list of component indices of position vector that
            ``dr`` vector acts upon

        """
        for atom in self._atoms:
            r = atom.r.tolist()
            for i, ri in enumerate(r[:]):
                if i in r_indices:
                    r[i] += dr[i]
            atom.x, atom.y, atom.z = r

    def clip_bounds(self, min_limit=None, max_limit=None,
                    abs_limit=None, r_indices=[0, 1, 2]):
        """Remove atoms outside the given limits along given dimension.

        Parameters
        ----------
        min_limit : {None, float}, optional
        max_limit : {None, float}, optional
        abs_limit : {None, float}, optional
        r_indices : sequence, optional

        """
        for atom in self._atoms[:]:
            r = atom.r.tolist()
            for i, ri in enumerate(r[:]):
                if i in r_indices:
                    if abs_limit is not None:
                        if abs(ri) <= abs_limit:
                            continue
                        else:
                            del self._atoms[self._atoms.index(atom)]
                            break
                    else:
                        if min_limit is not None:
                            if ri >= min_limit:
                                continue
                            else:
                                del self._atoms[self._atoms.index(atom)]
                                break
                        if max_limit is not None:
                            if ri <= max_limit:
                                continue
                            else:
                                del self._atoms[self._atoms.index(atom)]
                                break
