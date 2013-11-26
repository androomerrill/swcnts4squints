# -*- coding: utf-8 -*-
"""
=============================================================
swcnts4squints model (:mod:`swcnts4squints.s4s_model`)
=============================================================

.. currentmodule:: swcnts4squints.s4s_model

.. autosummary::
   :toctree: generated/

   Nanotube

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext'

from math import pi

#from collections import OrderedDict

import numpy as np

#from ..arrayfuncs import rotation_matrix
#from ..chemistry import Atom, Atoms
#from ..refdata import ccbond
#from .structure_io import XYZWriter

from .nanogen import Nanotube

__all__ = ['S4SModel']


class S4SModel(Nanotube):

    def __init__(self):

        self.observers = []

        super(S4SModel, self).__init__(n=10, m=10, verbose=True)
        print(self.Ch)

    def init(self):
        self._nzcells = 10
        print(self.Ch)

    def _compute_tube_params(self):

        self.notify_observers()

    @property
    def n(self):
        """Chiral index :math:`n`"""
        return self._n

    @n.setter
    def n(self, value):
        self._n = int(value)
        self._compute_tube_params()

    @property
    def m(self):
        """Chiral index :math:`m`"""
        return self._m

    @m.setter
    def m(self, value):
        self._m = int(value)
        self._compute_tube_params()

    @property
    def bond(self):
        """Bond length in **Angstroms**."""
        return self._bond

    @bond.setter
    def bond(self, value):
        self._bond = value

    @property
    def bundle_density(self):
        return self._bundle_density

    def compute_bundle_density(self):
        a_CC = self._ccbond
        m_C = self.cgs_mass_C
        n = self._n
        m = self._m
        d_vdw = None
        if self._n == self._m:
            d_vdw = 3.38
        elif (self._m == 0) or (self._n == 0):
            d_vdw = 3.41
        else:
            d_vdw = 3.39
        bundle_density = 8 * pi**2 * m_C * np.sqrt(n**2 + m**2 + n*m) / \
            (9 * np.sqrt(3) * (a_CC * 1e-8)**3 *
                (np.sqrt(n**2 + m**2 + n*m) +
                    pi * d_vdw / (np.sqrt(3) * a_CC))**2)
        return bundle_density

    def _update_bond_length_dependents(self):
        self._Ch = self.compute_Ch()
        self._dt = self._Ch / pi
        self._bundle_density = self.compute_bundle_density()
        self._T = self.compute_T()
        self.notify_observers()

    def _update_Ch_dependents(self):
        self._d = self.compute_d()
        self._dR = self.compute_dR()
        self._Ch = self.compute_Ch()
        self._dt = self._Ch / pi
        self._T = self.compute_T()
        self._N = self._Nhexs_per_cell = self.compute_Nhexs_per_cell()
        self._Natoms = self.compute_Natoms()
        self._bundle_density = self.compute_bundle_density()
        self._t1 = self.compute_t1()
        self._t2 = self.compute_t2()
        self._M, self._p, self._q = self.compute_R()
        self._chiral_angle = self.compute_chiral_angle()

    def _update_Nz_unit_cells(self, Ncells=None):
        if Ncells is None:
            #self._Nz_unit_cells = int(ceil(10 * self._L / self._T))
            self._nzcells = 10 * self._L / self._T
        else:
            self._nzcells = Ncells
        #self._update_Natoms_per_tube()
        #self._update_tube_mass()

        #self._update_L()

    #def _update_tube_mass(self):
    #    self._tube_mass = self.cgs_mass_C * \
    #        self._Natoms_per_cell * self._Nz_unit_cells

    def register_observer(self, observer):
        self.observers.append(observer)

    def remove_observer(self, observer):
        self.observers.remove(observer)

    def notify_observers(self):
        for observer in self.observers[:]:
            observer.update_app_view()
