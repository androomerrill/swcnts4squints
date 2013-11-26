# -*- coding: utf-8 -*-
"""
======================================================================
Sub-package for generating nano-structures (:mod:`tasr.tools.nanogen`)
======================================================================

.. currentmodule:: tasr.tools.nanogen

TASR package for generating nano-structures.

Contents
========

Abstract nanostructure objects
------------------------------

.. autosummary::
   :toctree: generated/

   Nanotube

Structure generator classes
---------------------------

.. autosummary::
   :toctree: generated/

   GrapheneGenerator
   BiLayerGrapheneGenerator
   NanotubeGenerator
   TubeGen

.. note::
   The :py:class:`~tasr.tools.nanogen.NanotubeGenerator` class
   does not yet generate nanotube *bundles*. Only single tubes.

.. seealso:: CLI module :py:mod:`tasr.scripts.nanogen`

Classes for changing structures
-------------------------------

.. autosummary::
   :toctree: generated/

   DefectGenerator
   VacancyGenerator

Sub-packages
------------

.. autosummary::
   :toctree: generated/

   structure_io

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

from .graphene import *
from .nanotube import *

__all__ = [s for s in dir() if not s.startswith('_')]
