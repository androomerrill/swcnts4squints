# -*- coding: utf-8 -*-
"""
=============================================================
swcnts4squints model (:mod:`swcnts4squints.s4s_model`)
=============================================================

.. currentmodule:: swcnts4squints.s4s_model

.. autosummary::
   :toctree: generated/

   S4SModel

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext'

from fractions import gcd
from math import pi

from collections import OrderedDict

import numpy as np

#from ..arrayfuncs import rotation_matrix
#from ..chemistry import Atom, Atoms
#from ..refdata import ccbond
#from .structure_io import XYZWriter

param_units = {}
param_units['dt'] = \
    param_units['rt'] = \
    param_units['Ch'] = \
    param_units['T'] = \
    param_units['bond'] = u' \u212B'
param_units['chiral_angle'] = u'\u00b0'

param_symbols = {}
param_symbols['dt'] = u'd\u209C'
param_symbols['rt'] = u'r\u209C'
param_symbols['Ch'] = u'C\u2095'
param_symbols['t1'] = u't\u2081'
param_symbols['t2'] = u't\u2082'
param_symbols['chiral_angle'] = u'\u03b8\u1d04'

param_strfmt = {}
param_strfmt['Ch'] = \
    param_strfmt['T'] = \
    param_strfmt['dt'] = \
    param_strfmt['rt'] = \
    param_strfmt['chiral_angle'] = '{:.2f}'
param_strfmt['bond'] = '{:.3f}'

mass_C = 12.0107
gram_per_kilogram = 1e3
kilogram_per_Dalton = 1.660538782e-27
gram_per_Dalton = kilogram_per_Dalton * gram_per_kilogram
cm_per_angstrom = 1e-8
ccbond = 1.4210

__all__ = ['Nanotube', 'param_units', 'param_symbols', 'param_strfmt']


class Nanotube(object):

    u"""Class for creating interactive Nanotube objects.

    .. versionadded:: 0.3.10

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    bond : float, optional
        bond length between nearest neighbor atoms.
        Must be in units of **Angstroms**. Default value is
        the carbon-carbon bond length in graphite:
        :math:`\\mathrm{a}_{\\mathrm{CC}} = 1.421` \u212b ([SoFaCNTs]_)
    verbose : bool, optional
        verbose output

    References
    ----------
    .. [SoFaCNTs] Science of Fullerenes and Carbon Nanotubes,
       M. Dresselhaus, G. Dresselhaus, P. Eklund, 1996, p. 760.

    Examples
    --------

    """

    def __init__(self, n=int, m=int, bond=ccbond, nxcells=1, nycells=1,
                 nzcells=1, verbose=False):

        self._params = OrderedDict()

        # add each parameter in the order I want them to appear in
        # verbose output mode
        self._params['n'] = {}
        self._params['m'] = {}
        self._params['t1'] = {}
        self._params['t2'] = {}
        self._params['d'] = {}
        self._params['dR'] = {}
        self._params['N'] = {}
        self._params['M'] = {}
        self._params['R'] = {}
        self._params['bond'] = {}
        self._params['Ch'] = {}
        self._params['T'] = {}
        self._params['dt'] = {}
        self._params['rt'] = {}
        self._params['chiral_angle'] = {}
        self._params['nxcells'] = {}
        self._params['nycells'] = {}
        self._params['nzcells'] = {}

        self._n = int(n)
        self._m = int(m)
        self._bond = float(bond)
        self._verbose = verbose

        self._d = None
        self._dR = None
        self._t1 = None
        self._t2 = None
        self._Ch = None
        self._T = None
        self._dt = None
        self._rt = None
        self._chiral_angle = None
        self._N = None
        self._M = None
        self._R = None
        self._p = None
        self._q = None
        self._Natoms = None
        self._Natoms_per_tube = None

        self.observers = []

        self.cgs_mass_C = mass_C * gram_per_Dalton

        # swcnt parameters
        self._nxcells = int(nxcells)
        self._nycells = int(nycells)
        self._nzcells = int(nzcells)
        self._tube_length = None
        self._tube_mass = None

        # bundle parameters
        self._bundle_density = None

        for k, v in self.__dict__.iteritems():
            p = k.strip('_')
            if p in self._params.keys():
                self._params[p]['units'] = param_units.get(p)
                self._params[p]['strfmt'] = param_strfmt.get(p)
                if param_symbols.get(p) is not None:
                    self._params[p]['var'] = param_symbols[p]
                else:
                    self._params[p]['var'] = p

        self._compute_tube_params()

    def _compute_tube_params(self):
        self._d = self.compute_d(n=self._n, m=self._m)
        self._dR = self.compute_dR(n=self._n, m=self._m)
        self._t1 = self.compute_t1(n=self._n, m=self._m)
        self._t2 = self.compute_t2(n=self._n, m=self._m)

        #Compute geometric properties
        self._Ch = self.compute_Ch(n=self._n, m=self._m, bond=self._bond)
        self._T = self.compute_T(n=self._n, m=self._m, bond=self._bond)
        self._dt = self.compute_dt(n=self._n, m=self._m, bond=self._bond)
        self._rt = self.compute_rt(n=self._n, m=self._m, bond=self._bond)
        self._chiral_angle = self.compute_chiral_angle(n=self._n, m=self._m)
        self._N = self.compute_N(n=self._n, m=self._m)
        self._R = self.compute_R(n=self._n, m=self._m)
        self._M = self.compute_M(n=self._n, m=self._m)
        self._Natoms = self.compute_Natoms(n=self._n, m=self._m)

        for k, v in self.__dict__.iteritems():
            p = k.strip('_')
            if p in self._params.keys():
                self._params[p]['val'] = v

        if self._verbose:
            for p, pdict in self._params.iteritems():
                pvar = pdict['var']
                pval = pdict['val']
                punits = pdict['units']
                pstrfmt = pdict['strfmt']
                if punits is not None:
                    if pstrfmt is not None:
                        print(u"{}: {}{}".format(pvar,
                                                 pstrfmt.format(pval),
                                                 punits))
                    else:
                        print(u"{}: {}{}".format(pvar, pval, punits))
                else:
                    print(u"{}: {}".format(pvar, pval))

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
    def nzcells(self):
        return self._nzcells

    @nzcells.setter
    def nzcells(self, value):
        self._nzcells = value
        self._update_Natoms_per_tube()
        self._update_L()
        self._update_tube_mass()
        self.notify_observers()

    @property
    def Natoms_per_tube(self):
        return self._Natoms_per_tube

    @property
    def tube_mass(self):
        return self._tube_mass

    @property
    def bundle_density(self):
        return self._bundle_density

    @property
    def t1(self):
        """:math:`t_{1} = \\frac{2m + n}{d_{R}}`
        where :math:`d_{R} = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_{1}`:

        .. math::

           \\mathbf{T} = t_{1}\\mathbf{a}_{1} + t_{2}\\mathbf{a}_{2}

        """
        return self._t1

    @classmethod
    def compute_t1(cls, n=int, m=int):
        """Compute :math:`t_{1} = \\frac{2m + n}{d_{R}}`
        where :math:`d_{R} = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_{1}`:

        .. math::

           \\mathbf{T} = t_{1}\\mathbf{a}_{1} + t_{2}\\mathbf{a}_{2}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.

        Returns
        -------
        int
            :math:`t_{1}`

        """
        dR = Nanotube.compute_dR(n=n, m=m)
        return int((2 * m + n) / dR)

    @property
    def t2(self):
        """:math:`t_{2} = -\\frac{2n + m}{d_{R}}`
        where :math:`d_{R} = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_{2}`:

        .. math::

           \\mathbf{T} = t_{1}\\mathbf{a}_{1} + t_{2}\\mathbf{a}_{2}

        """
        return self._t2

    @classmethod
    def compute_t2(cls, n=int, m=int):
        """Compute :math:`t_{2} = -\\frac{2n + m}{d_{R}}`
        where :math:`d_{R} = \\gcd{(2n + m, 2m + n)}`.

        The component of the translation vector :math:`\\mathbf{T}`
        along :math:`\\mathbf{a}_{2}`:

        .. math::

           \\mathbf{T} = t_{1}\\mathbf{a}_{1} + t_{2}\\mathbf{a}_{2}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.

        Returns
        -------
        int
            :math:`t_{2}`

        """
        dR = Nanotube.compute_dR(n=n, m=m)
        return -int((2 * n + m) / dR)

    @property
    def N(self):
        """Number of hexagons per nanotube unit cell :math:`N`."""
        return self._N

    @classmethod
    def compute_N(cls, n=int, m=int):
        """Compute :math:`N = \\frac{2(n^2+m^2+nm)}{d_{R}}`.
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.

        Returns
        -------
        int
            Number of hexagons per nanotube unit cell :math:`N`.

        """
        dR = Nanotube.compute_dR(n=n, m=m)
        return int(2 * (n**2 + m**2 + n * m) / dR)

    @property
    def Natoms(self):
        """Number of atoms per nanotube **unit cell** :math:`2N`."""

        return self._Natoms

    @classmethod
    def compute_Natoms(cls, n=int, m=int):
        """Compute :math:`N_{\mathrm{atoms/cell}} = 2N`.

        .. math::

           N_{\\mathrm{atoms}} = 2N = \\frac{4(n^2 + m^2 + nm)}{d_{R}}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.

        Returns
        -------
        int
            Number of atoms per nanotube **unit cell** :math:`2N`.

        """
        N = Nanotube.compute_N(n=n, m=m)
        return 2 * N

    @property
    def R(self):
        """Symmetry vector :math:`\\mathbf{R} = (p, q)`"""
        return self._R

    @classmethod
    def compute_R(cls, n=int, m=int, bond=None, magnitude=False):
        """Compute symmetry vector :math:`\\mathbf{R} = (p, q)`

        .. math::

           \\mathbf{R} = p\\mathbf{a}_{1} + q\\mathbf{a}_{2}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.
        magnitude : bool, optional
            .. versionadded:: 0.3.11
            if ``True``, return the magnitude of R

        Returns
        -------
        (p, q) : tuple
            2-tuple of ints which are the components of R vector
        float
            magnitude of R if ``magnitude`` is `True`

        """
        t1 = Nanotube.compute_t1(n=n, m=m)
        t2 = Nanotube.compute_t2(n=n, m=m)
        N = Nanotube.compute_N(n=n, m=m)

        p = None
        q = None

        for i in xrange(0, t1 + n + 1):
            for j in xrange(t2, m + 1):
                R = t1 * j - t2 * i
                if R == 1:
                    M = m * i - n * j
                    if M > 0 and M <= N:
                        p = i
                        q = j

        if magnitude:
            if bond is None:
                bond = ccbond
            return bond * np.sqrt(3 * (p**2 + q**2 + p * q))
        else:
            return (p, q)

    @property
    def M(self):
        """The number of :math:`\\mathbf{T}` in :math:`N\\mathbf{R}`"""
        return self._M

    @classmethod
    def compute_M(cls, n=int, m=int):
        """Compute :math:`M = mp - nq`

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.

        Returns
        -------
        int
            :math:`M = mp - nq`

        """
        p, q = Nanotube.compute_R(n=n, m=m)
        return m * p - n * q

    @property
    def Ch(self):
        """Nanotube circumference :math:`\\mathbf{C}_{h} = (n, m)`"""
        return self._Ch

    @classmethod
    def compute_Ch(cls, n=int, m=int, bond=None):
        """Compute the nanotube circumference.

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.
        bond : float, optional
            distance between nearest neighbor atoms.
            Must be in units of **Angstroms**.

        Returns
        -------
        float
            nanotube circumference in Angstroms

        """
        if bond is None:
            bond = ccbond

        return bond * np.sqrt(3 * (n**2 + m**2 + n * m))

    @property
    def dt(self):
        """Nanotube diameter :math:`d_{t} = \\frac{|\\mathbf{C}_{h}|}{\\pi}`"""
        return self._dt

    @classmethod
    def compute_tube_diameter(cls, n=int, m=int, bond=None):
        """Alias for :meth:`Nanotube.compute_dt`"""
        return Nanotube.compute_dt(n, m, bond)

    @classmethod
    def compute_dt(cls, n=int, m=int, bond=None):
        """Compute nanotube diameter :math:`d_{t}`

        .. math::

           d_{t} = \\frac{|\\mathbf{C}_{h}|}{\\pi}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} =
            (n, m)`.
        bond : float, optional
            distance between nearest neighbor atoms.
            Must be in units of **Angstroms**.

        Returns
        -------
        float
            nanotube diameter in Angstroms

        """
        Ch = Nanotube.compute_Ch(n, m, bond)
        return Ch / pi

    @property
    def rt(self):
        """Nanotube radius :math:`r_{t} = \\frac{|\\mathbf{C}_{h}|}{2\\pi}`"""
        return self._rt

    @classmethod
    def compute_tube_radius(cls, n=int, m=int, bond=None):
        """Alias for :meth:`Nanotube.compute_rt`"""
        return Nanotube.compute_rt(n, m, bond)

    @classmethod
    def compute_rt(cls, n=int, m=int, bond=None):
        """Compute nanotube radius :math:`r_{t}`

        .. math::

           r_{t} = \\frac{|\\mathbf{C}_{h}|}{2\\pi}


        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2}
            = (n, m)`.
        bond : float, optional
            distance between nearest neighbor atoms.
            Must be in units of **Angstroms**.

        Returns
        -------
        float
            nanotube radius in Angstroms

        """
        Ch = Nanotube.compute_Ch(n, m, bond)
        return Ch / (2 * pi)

    @property
    def d(self):
        """:math:`d=\\gcd{(n, m)}`"""
        return self._d

    @classmethod
    def compute_d(cls, n=int, m=int):
        """Compute :math:`d=\\gcd{(n, m)}`

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} +
            m\\mathbf{a}_{2} = (n, m)`.

        Returns
        -------
        int
            greatest common divisor of :math:`n` and :math:`m`

        """
        return gcd(n, m)

    @property
    def dR(self):
        """:math:`d_R=\\gcd{(2n + m, 2m + n)}`"""
        return self._dR

    @classmethod
    def compute_dR(cls, n=int, m=int):
        """Compute :math:`d_R=\\gcd{(2n + m, 2m + n)}`

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} +
            m\\mathbf{a}_{2} = (n, m)`.

        Returns
        -------
        int
            greatest common divisor of :math:`2n+m` and :math:`2m+n`

        """
        return gcd(2 * m + n, 2 * n + m)

    @property
    def chiral_angle(self):
        """Chiral angle :math:`\\theta_{c}`.

        .. math::

           \\theta_{c} = \\atan{\\frac{\\sqrt{3} m}{2n + m}}

        """
        return self._chiral_angle

    @classmethod
    def compute_chiral_angle(cls, n=int, m=int):
        """Compute chiral angle :math:`\\theta_{c}`

        .. math::

           \\theta_{c} = \\atan{\\frac{\\sqrt{3} m}{2n + m}}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} +
            m\\mathbf{a}_{2} = (n, m)`.

        """
        #return np.arccos((2*n + m) / (2 * np.sqrt(n**2 + m**2 + n*m)))
        return np.degrees(np.arctan(np.sqrt(3) * m / (2 * n + m)))

    @property
    def T(self):
        """Unit cell length.

        :math:`|\\mathbf{T}| = \\frac{\\sqrt{3} |\\mathbf{C}_{h}|}{d_{R}}`

        """
        return self._T

    @classmethod
    def compute_T(cls, n=None, m=None, Ch=None, dR=None, bond=None):
        """Compute unit cell length :math:`|\\mathbf{T}|`

        .. math::

           |\\mathbf{T}| = \\frac{\\sqrt{3} |\\mathbf{C}_{h}|}{d_{R}}

        Parameters
        ----------
        n, m : int
            Chiral indices defining the nanotube chiral vector
            :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} +
            m\\mathbf{a}_{2} = (n, m)`.
        Ch : float, optional
            nanotube circumference in Angstroms
        dR : int, optional
            greatest common divisor of :math:`2n + m` :math:`2m + n`
        bond : float, optional
            distance between nearest neighbor atoms.
            Must be in units of **Angstroms**.

        Returns
        -------
        float
            length of unit cell in Angstroms

        Raises
        ------
        TypeError
            if the parameters not valid.

        """
        if bond is None:
            bond = ccbond
        if n is not None and m is not None:
            Ch = Nanotube.compute_Ch(n=n, m=m, bond=bond)
            dR = Nanotube.compute_dR(n=n, m=m)
            return np.sqrt(3) * Ch / dR
        elif Ch is not None and dR is not None:
            return np.sqrt(3) * Ch / dR
        else:
            raise TypeError("Invalid parameters")

    @property
    def Ltube(self):
        return self._Ltube

    @Ltube.setter
    def Ltube(self, value):
        self._Ltube = value

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
        self._Natoms_per_cell = self.compute_Natoms_per_cell()
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
        self._update_Natoms_per_tube()
        self._update_tube_mass()
        #self._update_L()

    def _update_L(self):
        self._L = self._Nz_unit_cells * self._T / 10

    def _update_tube_mass(self):
        self._tube_mass = self.cgs_mass_C * \
            self._Natoms_per_cell * self._Nz_unit_cells

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
            (9 * np.sqrt(3) * (a_CC * cm_per_angstrom)**3 *
                (np.sqrt(n**2 + m**2 + n*m) +
                    pi * d_vdw / (np.sqrt(3) * a_CC))**2)
        return bundle_density

    def register_observer(self, observer):
        self.observers.append(observer)

    def remove_observer(self, observer):
        self.observers.remove(observer)

    def notify_observers(self):
        for observer in self.observers[:]:
            observer.update_app_view()
