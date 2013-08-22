from __future__ import division

import sys

from fractions import gcd
from math import acos, degrees, floor, pi, sqrt

import numpy as np

if sys.platform == 'linux2':
    from glpk.glpkpi import glp_create_prob, glp_load_matrix, \
            glp_set_obj_dir, glp_add_rows, glp_set_row_name, \
            glp_set_row_bnds, glp_add_cols, glp_set_col_name, \
            glp_set_col_bnds, glp_set_obj_coef, glp_set_col_kind, \
            glp_intopt, glp_simplex, glp_mip_col_val, glp_mip_obj_val, \
            GLP_MIN, GLP_MAX, GLP_UP, GLP_LO, GLP_DB, GLP_FX, \
            GLP_IV, GLP_MPS_FILE, intArray, doubleArray

#from moldycode.misc.constants import elements, mass_units
#from moldycode.tools.nanogen.nanotube import Nanotube
mass_C = 12.0107
gram_per_kilogram = 1e3
kilogram_per_Dalton = 1.660538782e-27
gram_per_Dalton = kilogram_per_Dalton * gram_per_kilogram
cm_per_angstrom = 1e-8


class S4SModel(object):

    def __init__(self):
        self.observers = []

        self.cgs_mass_C = mass_C * gram_per_Dalton

        self._ccbond = 1.421
        self._n = 10
        self._m = 10
        self._d = self.compute_d()
        self._dR = self.compute_dR()

        self._Ch = self.compute_Ch()
        self._T = self.compute_T()

        #self._Ch = None
        #self._T = None

        self._L = None
        self._dt = None
        self._chiral_angle = None
        self._t1 = None
        self._t2 = None

        # swcnt parameters
        self._Nz_unit_cells = None
        self._Nhexs_per_cell = None
        self._Natoms_per_cell = None
        self._Natoms_per_tube = None
        self._tube_mass = None

        self._p = None
        self._q = None
        self._M = None
        self._R = None
        self._N = None

        # bundle parameters
        self._bundle_density = None


    def init(self):
        #self._n = 10
        #self._m = 0
        #self._ccbond = 1.421
        #self._L = 10.
        self._Nz_unit_cells = 1
        self._update_L()
        self._update_params()

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, value):
        self._n = value
        self._update_Ch_dependents()
        self._update_Nz_unit_cells(Ncells=1)
        self._update_L()
        self._update_params()

    @property
    def m(self):
        return self._m

    @m.setter
    def m(self, value):
        self._m = value
        self._update_Ch_dependents()
        self._update_Nz_unit_cells(Ncells=1)
        self._update_L()
        self._update_params()

    @property
    def ccbond(self):
        return self._ccbond

    @ccbond.setter
    def ccbond(self, value):
        self._ccbond = value
        self._update_bond_length_dependents()
        self.notify_observers()

    @property
    def L(self):
        return self._L

    @L.setter
    def L(self, value):
        self._L = value
        self._update_Nz_unit_cells()
        self.notify_observers()

    @property
    def Nz_unit_cells(self):
        return self._Nz_unit_cells

    @Nz_unit_cells.setter
    def Nz_unit_cells(self, value):
        self._Nz_unit_cells = value
        self._update_Natoms_per_tube()
        self._update_L()
        self._update_tube_mass()
        self.notify_observers()

    @property
    def Ch(self):
        return self._Ch

    @property
    def Natoms_per_cell(self):
        return self._Natoms_per_cell

    @property
    def Natoms_per_tube(self):
        return self._Natoms_per_tube

    @property
    def chiral_angle(self):
        return self._chiral_angle

    @property
    def T(self):
        return self._T

    @property
    def dt(self):
        return self._dt

    @property
    def dR(self):
        return self._dR

    @property
    def d(self):
        return self._d

    @property
    def t1(self):
        return self._t1

    @property
    def t2(self):
        return self._t2

    @property
    def p(self):
        return self._p

    @property
    def q(self):
        return self._q

    @property
    def M(self):
        return self._M

    @property
    def tube_mass(self):
        return self._tube_mass

    @property
    def bundle_density(self):
        return self._bundle_density

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
        self._Nhexs_per_cell = self.compute_Nhexs_per_cell()
        self._Natoms_per_cell = self.compute_Natoms_per_cell()
        self._bundle_density = self.compute_bundle_density()
        self._t1 = self.compute_t1()
        self._t2 = self.compute_t2()
        self._M, self._p, self._q = self.compute_R()
        self._chiral_angle = self.compute_chiral_angle()

    def _update_L_dependents(self):
        pass

    def _update_Nz_unit_cells_dependents(self):
        pass

    def _update_Nz_unit_cells(self, Ncells=None):
        if Ncells is None:
            #self._Nz_unit_cells = int(ceil(10 * self._L / self._T))
            self._Nz_unit_cells = 10 * self._L / self._T
        else:
            self._Nz_unit_cells = Ncells
        self._update_Natoms_per_tube()
        self._update_tube_mass()
        #self._update_L()

    def _update_L(self):
        self._L = self._Nz_unit_cells * self._T / 10

    def _update_Natoms_per_tube(self):
        self._Natoms_per_tube = int(floor(self.Natoms_per_cell * self._Nz_unit_cells))

    def _update_tube_mass(self):
        self._tube_mass = self.cgs_mass_C * self._Natoms_per_cell * self._Nz_unit_cells

    def _update_params(self):
        #self._d = gcd(self._n, self._m)
        #self._dR = Nanotube.compute_dR(n=self._n, m=self._m)
        #self._Natoms_per_cell = Nanotube.compute_Natoms_per_cell(self._n, self._m)
        #self._tube_mass = elements['m']['C'] * self._Natoms_per_cell * mass_units['eV']
        #self._t1 = int((2 * self._m + self._n) / self._dR)
        #self._t2 = -int((2 * self._n + self._m) / self._dR)
        #self._Ch = Nanotube.compute_Ch(n=self._n, m=self._m, a=self._ccbond)
        #self._dt = self._Ch / pi
        #self._chiral_angle = Nanotube.compute_chiral_angle(n=self._n, m=self._m)
        #self._T = Nanotube.compute_T(Ch=self._Ch, dR=self._dR)
        self._d = self.compute_d()
        self._dR = self.compute_dR()
        self._Nhexs_per_cell = self.compute_Nhexs_per_cell()
        self._Natoms_per_cell = self.compute_Natoms_per_cell()
        self._bundle_density = self.compute_bundle_density()
        self._t1 = self.compute_t1()
        self._t2 = self.compute_t2()
        self._M, self._p, self._q = self.compute_R()
        self._Ch = self.compute_Ch()
        self._dt = self._Ch / pi
        self._chiral_angle = self.compute_chiral_angle()
        self._T = self.compute_T()
        #self._update_Nz_unit_cells()
        self.notify_observers()

    def compute_d(self):
        n = self._n
        m = self._m
        return gcd(n, m)

    def compute_dR(self):
        n = self._n
        m = self._m
        return gcd(2 * m + n, 2 * n + m)

    def compute_t1(self):
        n = self._n
        m = self._m
        dR = self._dR
        return int((2 * m + n) / dR)

    def compute_t2(self):
        n = self._n
        m = self._m
        dR = self._dR
        return -int((2 * n + m) / dR)

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
        bundle_density = 8 * pi**2 * m_C * sqrt(n**2 + m**2 + n*m) / \
                (9 * sqrt(3) * (a_CC * cm_per_angstrom)**3 * \
                (sqrt(n**2 + m**2 + n*m) + pi * d_vdw / (sqrt(3) * a_CC))**2)
        return bundle_density

    def compute_Ch(self):
        n = self._n
        m = self._m
        a_CC = self._ccbond
        return a_CC * sqrt(3 * (n**2 + m**2 + n*m))

    def compute_chiral_angle(self):
        n = self._n
        m = self._m
        return degrees(acos((2 * n + m) / (2 * sqrt(n**2 + m**2 + n*m))))

    def compute_T(self):
        Ch = self._Ch
        dR = self._dR
        return sqrt(3) * Ch / dR

    def compute_Natoms_per_cell(self):
        n = self._n
        m = self._m
        dR = self._dR
        return int(4 * (n**2 + m**2 + n*m) / dR)

    def compute_Nhexs_per_cell(self):
        n = self._n
        m = self._m
        dR = self._dR
        return int(2 * (n**2 + m**2 + n*m) / dR)

    def compute_R(self):
        if sys.platform == 'linux2':
            size = 1000+1
            ia = intArray(size)
            ja = intArray(size)
            ar = doubleArray(size)
            lp = glp_create_prob()
            glp_set_obj_dir(lp, GLP_MIN)
            glp_add_rows(lp, 2)
            glp_set_row_name(lp, 1, "N")
            glp_set_row_bnds(lp, 1, GLP_DB, 1, self._Nhexs_per_cell)
            glp_set_row_name(lp, 2, "1")
            glp_set_row_bnds(lp, 2, GLP_FX, 1, 1)
            glp_add_cols(lp, 2)
            glp_set_col_name(lp, 1, "p")
            glp_set_col_bnds(lp, 1, GLP_LO, 1, 0)
            glp_set_obj_coef(lp, 1, self._m)
            glp_set_col_name(lp, 2, "q")
            glp_set_col_bnds(lp, 2, GLP_UP, 0, 0)
            glp_set_obj_coef(lp, 2, -self._n)
            glp_set_col_kind(lp, 1, GLP_IV)
            glp_set_col_kind(lp, 2, GLP_IV)
            ia[1]=1; ja[1]=1; ar[1]=self._m
            ia[2]=1; ja[2]=2; ar[2]=-self._n
            ia[3]=2; ja[3]=1; ar[3]=-self._t2
            ia[4]=2; ja[4]=2; ar[4]=self._t1

            glp_load_matrix(lp, 4, ia, ja, ar)
            glp_simplex(lp, None)
            glp_intopt(lp, None)
            M = int(glp_mip_obj_val(lp))
            p = int(glp_mip_col_val(lp, 1))
            q = int(glp_mip_col_val(lp, 2))
            del lp
            #return (M, p, q)
            return (np.nan, np.nan, np.nan)
        else:
            return (None, None, None)

    def register_observer(self, observer):
        self.observers.append(observer)

    def remove_observer(self, observer):
        self.observers.remove(observer)

    def notify_observers(self):
        for observer in self.observers[:]:
            observer.update_app_view()
