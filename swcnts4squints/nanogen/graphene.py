# -*- coding: utf-8 -*-
"""
=================================================================
Graphene structure generator (:mod:`tasr.tools.nanogen.graphene`)
=================================================================

.. currentmodule:: tasr.tools.nanogen.graphene

.. autosummary::
   :toctree: generated/

   GrapheneGenerator
   BiLayerGrapheneGenerator

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext'

from math import ceil, sqrt
import itertools
import sys

import numpy as np

from ..arrayfuncs import rotation_matrix
from ..chemistry import Atom, Atoms
from ..refdata import ccbond
from .structure_io import XYZWriter

__all__ = ['GrapheneGenerator', 'BiLayerGrapheneGenerator']


class GrapheneGenerator(object):

    """Class for generating `n`-layer graphene nanostructures.

    .. versionadded:: 0.3.8

    Parameters
    ----------
    width : float
        Width of graphene sheet in **nanometers**
    length : float
        Length of graphene sheet in **nanometers**
    edge : {'AC', 'armchair', 'ZZ', 'zigzag'}, optional
        .. versionchanged:: 0.3.10
           Now recognizes keyword values ``armchair`` and ``zigzag``.
        **A**\ rm\ **C**\ hair or **Z**\ ig\ **Z**\ ag edge along
        the ``length`` of the sheet sheet.
    element1, element2 : {str, int}, optional
        .. versionadded:: 0.3.14
           Element symbol or atomic number of basis atoms 1 and 2.
    bond : float, optional
        bond length between nearest-neighbor atoms in **Angstroms**.
    nlayers : int, optional
        Number of graphene layers.
    layer_spacing : float, optional
        Distance between layers in **Angstroms**.
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers
    autogen : bool, optional
        automatically generate unit cell and full structure

        .. versionchanged:: 0.3.20

           if ``True``, also calls the ``generate_structure`` method.
    verbose : bool, optional
        verbose output

    Notes
    -----
    For now, the graphene structure is generated using a
    conventional unit cell, not the primitive unit cell.

    .. todo::

       Add notes on unit cell calculation.

    Examples
    --------

    Start an interactive python or ipython session, then import the
    GrapheneGenerator class.

    >>> from tasr.tools.nanogen import GrapheneGenerator

    Now generate a ``1 nm x 20 nm`` armchair edge graphene nano-ribbon.

    >>> ACG = GrapheneGenerator(width=1, length=20, edge='armchair')

    Save structure data in `xyz` format:

    >>> ACG.save_data(fname='1nmx20nm_AC_edge.xyz')

    The rendered structure look like:

    .. image:: /images/1nmx20nm_AC_edge.png

    Now let's generate a ``1 nm x 20 nm`` zigzag edge graphene nano-ribbon.

    >>> ZZG = GrapheneGenerator(width=1, length=20, edge='zigzag')
    >>> ZZG.save_data(fname='1nmx20nm_ZZ_edge.xyz')

    The rendered structure looks like:

    .. image:: /images/1nmx20nm_ZZ_edge.png

    Now generate ``5 nm`` by ``25 nm``, ``armchair`` edge,
    5 layer, ``AB``-stacked graphene.

    >>> ACG_5layers = GrapheneGenerator(width=5, length=25,
    ...                                 edge='armchair', nlayers=5)
    >>> ACG_5layers.save_data(fname='5nmx25nm_5layer_AC_graphene.xyz')

    The rendered structure looks like:

    .. image:: /images/5nmx25nm_5layer_AC_graphene.png

    Now generate single layer, ``10 nm x 10 nm`` sheet of BN Graphene.

    >>> BN_graphene = GrapheneGenerator(width=10, length=10, edge='AC',
    ...                                 element1='B', element2='N')
    >>> BN_graphene.save_data(fname='10nmx10nm_1_layer_BN_graphene.xyz')

    The rendered structure looks like:

    .. image:: /images/10nmx10nm_single_layer_BN_graphene.png

    Now, just because we can, generate a ``5 nm x 5 nm`` sheet of
    Uranium-Einsteinium Graphene.

    >>> UEs_graphene = GrapheneGenerator(width=5, length=5, edge='zigzag',
    ...                                  element1='U', element2='Es')
    >>> UEs_graphene.save_data(fname='5nmx5nm_1_layer_UEs_graphene.xyz')

    The rendered structure looks like:

    .. image:: /images/5nmx5nm_single_layer_UEs_graphene.png

    """

    def __init__(self, width=float, length=float, edge='armchair',
                 element1='C', element2='C', bond=ccbond, nlayers=1,
                 layer_spacing=3.35, stacking_order='AB',
                 autogen=True, verbose=False):

        self.element1 = element1
        self.element2 = element2

        self.Lx = width
        self.Ly = length
        self.edge = edge
        self.bond = bond
        self.verbose = verbose

        self.lx = 0.
        self.ly = 0.

        self.Nx = 0
        self.Ny = 0

        self.nlayers = nlayers
        self.layer_spacing = layer_spacing
        self.stacking_order = stacking_order

        self.layer_shift = np.zeros(3)

        if nlayers > 1 and stacking_order == 'AB':
            if edge in ('AC', 'armchair'):
                self.layer_shift[1] = self.bond
            elif edge in ('ZZ', 'zigzag'):
                self.layer_shift[0] = self.bond
            else:
                print('unrecognized edge parameter: {}'.format(edge))
                sys.exit(1)

        self.atom1 = Atom(element1)
        self.atom2 = Atom(element2)
        self.atom3 = Atom(element1)
        self.atom4 = Atom(element2)

        self.Natoms = 0
        self.atoms = Atoms(atoms=[self.atom1,
                                  self.atom2,
                                  self.atom3,
                                  self.atom4])
        self.structure_atoms = None

        if autogen:
            self.generate_unit_cell()
            self.generate_structure()

    def generate_unit_cell(self):
        """Generate the unit cell.

        Called automatically if ``autogen`` is True.

        """

        if self.edge in ('AC', 'armchair'):
            # Set up the unit cell with the armchair edge aligned
            # along the `y`-axis.
            self.lx = sqrt(3) * self.bond
            self.ly = 3 * self.bond

            # Set up 4 atom basis
            # Leave atom 1 at the origin

            # Move atom 2 to 2nd basis position
            self.atom2.x = -sqrt(3) / 2 * self.bond
            self.atom2.y = self.bond / 2

            # Move atom 3 along a1 primitive vector
            self.atom3.x = -sqrt(3) / 2 * self.bond
            self.atom3.y = 3 / 2 * self.bond

            # Move atom 4 from atom 2 along a2 primitive vector
            self.atom4.y = 2 * self.bond

        elif self.edge in ('ZZ', 'zigzag'):
            # Set up the unit cell with the zigzag edge aligned
            # along the `y`-axis.
            self.lx = 3 * self.bond
            self.ly = sqrt(3) * self.bond

            # Set up 4 atom basis
            # Leave atom 1 at the origin

            # Move atom 2 to the right
            self.atom2.x = self.bond

            # Move atom 3 to the left and up
            self.atom3.x = 3 / 2 * self.bond
            self.atom3.y = sqrt(3) / 2 * self.bond

            # Move atom 4 to the right and up
            self.atom4.x = -self.bond / 2
            self.atom4.y = sqrt(3) / 2 * self.bond

        else:
            print('unrecognized edge parameter: {}'.format(self.edge))
            sys.exit(1)

        self.Nx = int(ceil(10 * self.Lx / self.lx))
        self.Ny = int(ceil(10 * self.Ly / self.ly))

    def generate_structure(self):
        """Generate the full structure coordinates.

        .. versionchanged:: 0.3.20

           Now called automatically if ``autogen`` is True.

        """

        self.structure_atoms = []

        for nlayer in xrange(self.nlayers):
            layer_atoms = Atoms()
            for nx in xrange(self.Nx):
                for ny in xrange(self.Ny):
                    dr = np.array([nx * self.lx,
                                   ny * self.ly,
                                   nlayer * self.layer_spacing])

                    for atom in self.atoms.atomlist:
                        layer_atom = Atom(atom.symbol)
                        layer_atom.r = atom.r + dr
                        layer_atoms.append(layer_atom)
                        self.Natoms += 1

            # translate layer to put its center of mass at the origin
            layer_atoms.center_CM(r_indices=[0, 1])
            if (nlayer % 2) != 0:
                layer_atoms.translate(self.layer_shift)

            self.structure_atoms.append(layer_atoms)

        if self.verbose:
            print(self.structure_atoms)

    def save_data(self, fname=None, structure_format='xyz', rotation_angle=-90,
                  rot_axis='x', deg2rad=True, center_CM=True):
        """Save structure data.

        Parameters
        ----------
        fname : str, optional
            file name string
        structure_format : str, optional
            chemical file format of saved structure data.
        rotation_angle : {None, float}, optional
            Angle of rotation
        rot_axis : {'x', 'y', 'z'}, optional
            Rotation axis
        deg2rad : bool, optional
            Convert ``rotation_angle`` from degrees to radians.
        center_CM : bool, optional
            .. versionadded:: 0.3.21

            Center center-of-mass on on origin.

        """
        structure_atoms = list(itertools.chain(*self.structure_atoms))
        structure_atoms = Atoms(structure_atoms)
        if center_CM and self.nlayers > 1:
            structure_atoms.center_CM(r_indices=[2])
        if rotation_angle is not None:
            R_matrix = rotation_matrix(rotation_angle,
                                       rot_axis=rot_axis,
                                       deg2rad=deg2rad)
            structure_atoms.rotate(R_matrix)
            #for layer_atoms in self.structure_atoms:
            #    layer_atoms.rotate(R_matrix)
        if fname is None:
            dimensions = '{}nmx{}nm'.format(self.Lx, self.Ly)
            nlayer = '{}layer'.format(self.nlayers)
            edge = 'AC' if self.edge in ('AC', 'armchair') else 'ZZ'
            atombond = '{}{}'.format(self.atom1.symbol, self.atom2.symbol)
            fname_wordlist = (dimensions, nlayer, edge, atombond,
                              '.'.join(('graphene', structure_format)))
            fname = '_'.join(fname_wordlist)
        if structure_format == 'xyz':
            XYZWriter.write(fname=fname, atoms=structure_atoms)


class BiLayerGrapheneGenerator(GrapheneGenerator):

    """Class for generating bi-layer graphene structures.

    .. versionadded:: 0.3.9

    Parameters
    ----------
    width : float
        Width of graphene sheet in **nanometers**
    length : float
        Length of graphene sheet in **nanometers**
    edge : {'AC', 'armchair', 'ZZ', 'zigzag'}, optional
        .. versionchanged:: 0.3.10
           Now recognizes keyword values ``armchair`` and ``zigzag``.
        **A**\ rm\ **C**\ hair or **Z**\ ig\ **Z**\ ag edge along
        the ``length`` of the sheet sheet.
    element1, element2 : {str, int}, optional
        .. versionadded:: 0.3.14
           Element symbol or atomic number of basis atoms 1 and 2.
    bond : float, optional
        bond length between nearest-neighbor atoms in **Angstroms**.
    layer_spacing : float, optional
        Distance between layers in **Angstroms**.
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers
    rotation_angle : {None, float}, optional
        .. versionchanged:: 0.3.10
           Default changed from radians to degrees.
        Rotation angle of second layer specified in degrees.
        If specified in radians, then you must set ``deg2rad=False``
    deg2rad : bool, optional
        .. versionchanged:: 0.3.10
           Default changed from ``False`` to ``True``
        The rotation angle is specified in degrees and needs to be converted
        to radians.
    autogen : bool, optional
        automatically generate unit cell and full structure

        .. versionchanged:: 0.3.20

           if ``True``, also calls the ``generate_structure`` method.
    verbose : bool, optional
        verbose output

    Examples
    --------

    Import the BiLayerGrapheneGenerator class

    >>> from tasr.tools.nanogen import BiLayerGrapheneGenerator

    Generate ``1 nm`` wide by ``10 nm`` long ``AB`` stacked
    bilayer-graphene with a ``ZZ`` edge:

    >>> bi_graphene = BiLayerGrapheneGenerator(width=1, length=10, edge='ZZ')

    Save structure data in `xyz` format:

    >>> bi_graphene.save_data(fname='1nmx10nm_bilayer.xyz')

    The rendered structure looks like (after rotating 90 degrees so that
    it better fits the page):

    .. image:: /images/1nmx10nm_bilayer.png

    Now generate bilayer-graphene with top layer rotated by 45 degrees.

    >>> rotated_bilayer = BiLayerGrapheneGenerator(width=10, length=10,
    ...                                            edge='armchair',
    ...                                            rotation_angle=45)
    >>> rotated_bilayer.save_data(fname='bilayer_rotation=45deg.xyz')

    The rendered structure looks like:

    .. image:: /images/rotated_bilayer.png

    Now generate BN bilayer-graphene with top layer rotated 45 degrees.

    >>> rotated_BN_bilayer = BiLayerGrapheneGenerator(width=10, length=10,
    ...                                               edge='zigzag',
    ...                                               element1='B',
    ...                                               element2='N',
    ...                                               rotation_angle=45)
    >>> rotated_BN_bilayer.save_data(fname='BN_bilayer_rotated_45deg.xyz')

    The rendered structure looks like:

    .. image:: /images/BN_bilayer_rotated_45deg.png

    """

    def __init__(self, width=float, length=float, edge='armchair',
                 element1='C', element2='C', bond=ccbond,
                 layer_spacing=3.35, stacking_order='AB',
                 rotation_angle=None, deg2rad=True, autogen=True,
                 verbose=False):

        super(BiLayerGrapheneGenerator, self).__init__(
            width=width, length=length, edge=edge,
            element1=element1, element2=element2, bond=bond,
            nlayers=2, layer_spacing=layer_spacing, autogen=False,
            verbose=verbose)

        self.rotation_matrix = None
        if rotation_angle is not None:
            self.rotation_matrix = rotation_matrix(rotation_angle,
                                                   rot_axis='z',
                                                   deg2rad=deg2rad)

        if autogen:
            super(BiLayerGrapheneGenerator, self).generate_unit_cell()
            self.generate_structure()

    def generate_structure(self):
        """Generate the full structure coordinates

        .. versionchanged:: 0.3.20

           Now called automatically if ``autogen`` is True.

        """
        super(BiLayerGrapheneGenerator, self).generate_structure()

        if self.rotation_matrix is not None:
            for n in xrange(self.nlayers):
                if (n % 2) != 0:
                    self.structure_atoms[n].rotate(self.rotation_matrix)
