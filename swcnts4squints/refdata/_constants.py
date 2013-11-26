"""
====================================================
Referece data (:mod:`tasr.tools.refdata._constants`)
====================================================

Collection of physical constants and conversion factors
built upon :mod:`scipy.constants` module.

.. currentmodule:: tasr.tools.refdata._constants

"""
from __future__ import division, print_function, absolute_import
from scipy.constants import codata
from scipy.constants.constants import nano

h = codata.value('Planck constant in eV s')
c = codata.value('speed of light in vacuum')
hc = h * c / nano

ccbond = 1.421  # angstroms

mass_units = {
    'kg': codata.value('atomic mass constant'),
    'eV': codata.value('atomic mass unit-electron volt relationship')
}

__all__ = ['h', 'c', 'hc', 'ccbond', 'mass_units']
