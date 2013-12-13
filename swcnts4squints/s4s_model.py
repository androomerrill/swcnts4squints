# -*- coding: utf-8 -*-
"""
=============================================================
swcnts4squints model (:mod:`swcnts4squints.s4s_model`)
=============================================================

.. currentmodule:: swcnts4squints.s4s_model

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext'

from sknano.nanogen import NanotubeBundle

__all__ = ['S4SModel']


class S4SModel(NanotubeBundle):

    def __init__(self):

        self.observers = []

        super(S4SModel, self).__init__(n=10, m=10, nzcells=10, verbose=True)

        self._fix_nzcells = True
        self._fix_tube_length = False

    def init(self):
        #self._nzcells = 10
        self.notify_observers()

    def _compute_bundle_params(self):
        super(S4SModel, self).compute_bundle_params()
        self.notify_observers()

    @property
    def fix_nzcells(self):
        return self._fix_nzcells

    @fix_nzcells.setter
    def fix_nzcells(self, value):
        self._fix_nzcells = value

    @property
    def fix_tube_length(self):
        return self._fix_tube_length

    @fix_tube_length.setter
    def fix_tube_length(self, value):
        self._fix_tube_length = value

    @property
    def n(self):
        """Chiral index :math:`n`"""
        return self._n

    @n.setter
    def n(self, value):
        self._n = int(value)
        self._compute_bundle_params()

    @property
    def m(self):
        """Chiral index :math:`m`"""
        return self._m

    @m.setter
    def m(self, value):
        self._m = int(value)
        self._compute_bundle_params()

    @property
    def bond(self):
        return self._bond

    @bond.setter
    def bond(self, value):
        self._bond = value
        self._compute_bundle_params()

    @property
    def tube_length(self):
        return self._tube_length

    @tube_length.setter
    def tube_length(self, value):
        self._tube_length = value
        self._compute_bundle_params()

    @property
    def nzcells(self):
        return self._nzcells

    @nzcells.setter
    def nzcells(self, value):
        self._nzcells = value
        self._compute_bundle_params()

    @property
    def Ntubes(self):
        return self._Ntubes

    @Ntubes.setter
    def Ntubes(self, value=int):
        self._Ntubes = value
        self._nxcells = value
        self._nycells = 1
        self._compute_bundle_params()

    #def _update_bond_length_dependents(self):
    #    self._Ch = self.compute_Ch()
    #    self._dt = self._Ch / pi
    #    self._bundle_density = self.compute_bundle_density()
    #    self._T = self.compute_T()
    #    self.notify_observers()

    #def _update_Ch_dependents(self):
    #    self._d = self.compute_d()
    #    self._dR = self.compute_dR()
    #    self._Ch = self.compute_Ch()
    #    self._dt = self._Ch / pi
    #    self._T = self.compute_T()
    #    self._N = self._Nhexs_per_cell = self.compute_Nhexs_per_cell()
    #    self._Natoms = self.compute_Natoms()
    #    self._bundle_density = self.compute_bundle_density()
    #    self._t1 = self.compute_t1()
    #    self._t2 = self.compute_t2()
    #    self._M, self._p, self._q = self.compute_R()
    #    self._chiral_angle = self.compute_chiral_angle()

    #def _update_Nz_unit_cells(self, Ncells=None):
    #    if Ncells is None:
    #        #self._Nz_unit_cells = int(ceil(10 * self._L / self._T))
    #        self._nzcells = 10 * self._L / self._T
    #    else:
    #        self._nzcells = Ncells
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
