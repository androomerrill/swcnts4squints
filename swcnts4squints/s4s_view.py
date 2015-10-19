# -*- coding: utf-8 -*-
"""
=============================================================
swcnts4squints view (:mod:`swcnts4squints.s4s_view`)
=============================================================

.. currentmodule:: swcnts4squints.s4s_view

"""
from __future__ import division, print_function, absolute_import
from __future__ import unicode_literals
from builtins import str
__docformat__ = 'restructuredtext'

from PyQt4.QtCore import pyqtSlot
from PyQt4.QtGui import QMainWindow

from .ui_s4s import Ui_S4S

__all__ = ['S4SView']


class S4SView(QMainWindow, Ui_S4S):
    """swcnts4squints MVC view class.

    Parameters
    ----------
    controller : :py:class:`swcnts4squints.s4s_controller` instance
    model : :py:class:`swcnts4squints.s4s_model` instance

    """
    def __init__(self, controller=None, model=None):
        self.controller = controller
        self.model = model
        model.register_observer(self)
        super(S4SView, self).__init__()
        self.setupUi(self)

    def init(self):
        self.show()

    @pyqtSlot(int)
    def on_n_spinBox_valueChanged(self, value):
        self.model.n = value

    @pyqtSlot(int)
    def on_m_spinBox_valueChanged(self, value):
        self.model.m = value

    @pyqtSlot()
    def on_bond_doubleSpinBox_editingFinished(self):
        self.model.bond = self.bond_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_bond_doubleSpinBox_valueChanged(self, value):
        self.model.bond = value

    @pyqtSlot()
    def on_nzcells_doubleSpinBox_editingFinished(self):
        self.model.nzcells = self.nzcells_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_nzcells_doubleSpinBox_valueChanged(self, value):
        self.model.nzcells = value

    @pyqtSlot(int)
    def on_fix_buttonGroup_buttonClicked(self):
        if self.fix_Ncells_radioButton.isChecked():
            self.model.fix_nzcells = True
            self.model.fix_tube_length = False
            self.nzcells_doubleSpinBox.setReadOnly(True)
            self.tube_length_lineEdit.setReadOnly(False)
        else:
            self.model.fix_nzcells = False
            self.model.fix_tube_length = True
            self.nzcells_doubleSpinBox.setReadOnly(False)
            self.tube_length_lineEdit.setReadOnly(True)

    @pyqtSlot()
    def on_tube_length_lineEdit_editingFinished(self):
        self.model.tube_length = float(self.tube_length_lineEdit.text())

    #@pyqtSlot()
    #def on_menuFile_Notes_on_Numbers_triggered(self):
    #    #notes = QFile(QString(":/max.density.calc.pdf"))
    #    if sys.platform == 'win32':
    #        try:
    #            import os
    #            os.startfile(notes.fileName())
    #        except WindowsError as e:
    #            print(e)
    #    elif sys.platform == 'darwin':
    #        try:
    #            import subprocess
    #            subprocess.call(["open", notes.fileName()])
    #        except Exception as e:
    #            print(e)

    def update_app_view(self):
        self.n_spinBox.setValue(self.model.n)
        self.m_spinBox.setValue(self.model.m)

        self.Ch_lineEdit.setText('{:.3f}'.format(self.model.Ch))
        self.dt_lineEdit.setText('{:.3f}'.format(self.model.dt))
        self.T_lineEdit.setText('{:.3f}'.format(self.model.T))
        self.chiral_angle_lineEdit.setText(
            '{:.2f}'.format(self.model.chiral_angle))
        self.N_lineEdit.setText(str(self.model.N))
        self.Natoms_lineEdit.setText(str(self.model.Natoms))
        self.bond_doubleSpinBox.setValue(self.model.bond)

        self.nzcells_doubleSpinBox.setValue(self.model.nzcells)
        self.tube_length_lineEdit.setText(
            '{:.3f}'.format(self.model.tube_length))
        self.tube_mass_lineEdit.setText(
            '{:.3e}'.format(self.model.tube_mass))
        self.Natoms_per_tube_lineEdit.setText(
            str(self.model.Natoms_per_tube))
        #self.Ntubes_mantissa_spinBox.setValue(
        #    str(self.model.Ntubes))
        self.bundle_mass_lineEdit.setText(
            str(self.model.bundle_mass))
        self.bundle_density_lineEdit.setText(
            '{:.3f}'.format(self.model.bundle_density))

        self.d_lineEdit.setText(str(self.model.d))
        self.dR_lineEdit.setText(str(self.model.dR))
        self.t1_lineEdit.setText(str(self.model.t1))
        self.t2_lineEdit.setText(str(self.model.t2))
