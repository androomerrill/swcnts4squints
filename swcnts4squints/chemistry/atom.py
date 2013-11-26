# -*- coding: utf-8 -*-
"""
==================================================================
Abstract representation of Atom (:mod:`tasr.tools.chemistry.atom`)
==================================================================

.. currentmodule:: tasr.tools.chemistry.atom

.. autosummary::
   :toctree: generated/

   Atom

"""
from __future__ import division, absolute_import, print_function
__docformat__ = 'restructuredtext'

import numpy as np

from ..refdata import atomic_numbers, atomic_masses, element_symbols


class AtomAttributes(object):
    pass


class Atom(object):

    """Class for creating abstract object representing an Atom.

    .. versionadded:: 0.3.11

    Parameters
    ----------
    element : {str, int}
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    atomID : int, optional
        atom ID, a LAMMPS atom attribute
    moleculeID : int, optional
        molecule ID, a LAMMPS atom attribute
    q : {int, float}, optional
        net charge of Atom
    atomtype : int, optional
        atom type, a LAMMPS atom attribute

    """

    def __init__(self, element, atomID=0, moleculeID=0, q=0, atomtype=1):
        self._check_type(element, (int, float, str))
        self._element = element
        self._check_type(q, (int, float))
        #self._q = np.asarray([q], dtype=float)
        self._q = q
        self._atomID = atomID
        self._moleculeID = moleculeID
        self._atomtype = atomtype
        self._symbol = None
        self._Z = None
        self._m = None
        self._idx = None
        self._r = np.zeros(3, dtype=float)

        self._properties = {'symbol': self._symbol, 'Z': self._Z,
                            'm': self._m, 'q': self._q, 'r': self._r}

        if isinstance(element, (int, float)):
            self._Z = int(element)
            self._idx = self._Z - 1
            try:
                self._symbol = element_symbols[self._idx]
                self._m = atomic_masses[self._symbol]
            except KeyError as e:
                print(e)
                print('unrecognized element number: {}'.format(element))
        elif isinstance(element, str):
            self._symbol = element
            try:
                self._Z = atomic_numbers[self._symbol]
                self._idx = self._Z - 1
                self._m = atomic_masses[self._symbol]
            except KeyError as e:
                print(e)
                print('Unrecognized atomic symbol: {}'.format(element))

    def __str__(self):
        """Return string representation of atom."""
        return self._symbol

    def _check_type(self, value, valid_types):
        if not isinstance(value, valid_types):
            raise TypeError('{} not valid type.\n'.format(value) +
                            '(Valid Types: {})'.format(valid_types))

    @property
    def m(self):
        """Return atomic mass of atom."""
        return self._m

    @property
    def Z(self):
        """Return atomic number of atom."""
        return self._Z

    @property
    def symbol(self):
        """Return element symbol of atom."""
        return self._symbol

    @property
    def x(self):
        """Return the x-coordinate of atom in units of **Angstroms**.

        Returns
        -------
        float
            x-coordinate of atom in units of **Angstroms**.

        """
        return self._r[0]

    @x.setter
    def x(self, value):
        """Set x-coordinate.

        Parameters
        ----------
        value : float
            x-coordinate in units of **Angstroms**

        """
        self._check_type(value, (int, float))
        self._r[0] = float(value)

    @property
    def y(self):
        """Return the y-coordinate of atom in units of **Angstroms**.

        Returns
        -------
        float
            y-coordinate of atom in units of **Angstroms**.

        """
        return self._r[1]

    @y.setter
    def y(self, value):
        """Set y-coordinate of atom.

        Parameters
        ----------
        value : float
            y-coordinate in units of **Angstroms**

        """
        self._check_type(value, (int, float))
        self._r[1] = float(value)

    @property
    def z(self):
        """Return the z-coordinate of atom in units of **Angstroms**.

        Returns
        -------
        float
            z-coordinate of atom in units of **Angstroms**.

        """
        return self._r[2]

    @z.setter
    def z(self, value):
        """Set z-coordinate of atom.

        Parameters
        ----------
        value : float
            z-coordinate in units of **Angstroms**

        """
        self._check_type(value, (int, float))
        self._r[2] = float(value)

    @property
    def r(self):
        """Return x-, y-, z- coordinates of atom in units of **Angstroms**.

        Returns
        -------
        ndarray
            3-element ndarray of [x, y, z] coordinates of atom.

        """
        return self._r

    @r.setter
    def r(self, value):
        """Set x, y, z coordinates of atom.

        Parameters
        ----------
        value : array_like
            3-element array of x-, y-, z-coordinates in units of **Angstroms**.

        """
        self._check_type(value, np.ndarray)
        for i, v in enumerate(value):
            self._r[i] = v

    @property
    def q(self):
        """Return net charge of atom."""
        return self._q

    @q.setter
    def q(self, value):
        """Set net charge of atom.

        Parameters
        ----------
        value : {int, float}
            net charge on atom as a multiple of the ``elementary_charge`` *e*.

        """
        self._check_type(value, (int, float))
        self._q = value

    @property
    def atomID(self):
        """Return atom ID of atom."""
        return self._atomID

    @atomID.setter
    def atomID(self, value):
        """Set atom ID of atom.

        Parameters
        ----------
        value : int
            atom ID

        """
        self._check_type(value, (int, float))
        self._atomID = int(value)

    @property
    def moleculeID(self):
        """Return molecule ID of atom."""
        return self._moleculeID

    @moleculeID.setter
    def moleculeID(self, value):
        """Set molecule ID of atom.

        Parameters
        ----------
        value : int
            molecule ID

        """
        self._check_type(value, (int, float))
        self._moleculeID = int(value)

    @property
    def atomtype(self):
        """Return atom type of atom."""
        return self._atomtype

    @atomtype.setter
    def atomtype(self, value):
        """Set atom type of atom.

        Parameters
        ----------
        value : int
            atom type

        """
        self._check_type(value, (int, float))
        self._atomtype = int(value)

    def fix_minus_zero_coords(self, epsilon=1.0e-10):
        """Set really really small negative coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon zero so we don't end up with -0.00000
        coordinates in structure data output.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any component of atom's
            position

        """
        r = self._r.tolist()
        for i, ri in enumerate(r[:]):
            if abs(ri) < epsilon:
                r[i] = 0.0
        self._r[0], self._r[1], self._r[2] = r
