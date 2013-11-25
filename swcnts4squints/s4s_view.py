from __future__ import absolute_import, print_function, division

from PyQt4.QtCore import pyqtSlot
from PyQt4.QtGui import QMainWindow

from .ui_s4s import Ui_S4S


class S4SView(QMainWindow, Ui_S4S):

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
    def on_ccbond_doubleSpinBox_editingFinished(self):
        self.model.ccbond = self.ccbond_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_ccbond_doubleSpinBox_valueChanged(self, value):
        self.model.ccbond = value

    @pyqtSlot()
    def on_Nz_unit_cells_doubleSpinBox_editingFinished(self):
        self.model.Nz_unit_cells = self.Nz_unit_cells_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_Nz_unit_cells_doubleSpinBox_valueChanged(self, value):
        self.model.Nz_unit_cells = value

    #@pyqtSlot(str)
    #def on_L_lineEdit_textEdited(self, text):
    #    self.model.L = float(text)

    @pyqtSlot()
    def on_L_lineEdit_editingFinished(self):
        self.model.L = float(self.L_lineEdit.text())

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
        self.ccbond_doubleSpinBox.setValue(self.model.ccbond)
        self.Nz_unit_cells_doubleSpinBox.setValue(self.model.Nz_unit_cells)
        self.L_lineEdit.setText('{:.3f}'.format(self.model.L))
        self.Ch_lineEdit.setText('{:.3f}'.format(self.model.Ch))
        self.chiral_angle_lineEdit.setText('{:.2f}'.format(self.model.chiral_angle))
        self.dt_lineEdit.setText('{:.3f}'.format(self.model.dt))
        self.T_lineEdit.setText('{:.3f}'.format(self.model.T))
        self.d_lineEdit.setText(str(self.model.d))
        self.dR_lineEdit.setText(str(self.model.dR))
        self.Natoms_per_cell_lineEdit.setText(str(self.model.Natoms_per_cell))
        self.Natoms_per_tube_lineEdit.setText(str(self.model.Natoms_per_tube))
        self.t1_lineEdit.setText(str(self.model.t1))
        self.t2_lineEdit.setText(str(self.model.t2))
        self.p_lineEdit.setText(str(self.model.p))
        self.q_lineEdit.setText(str(self.model.q))
        self.M_lineEdit.setText(str(self.model.M))
        self.mass_lineEdit.setText('{:.4e}'.format(self.model.tube_mass))
        self.density_lineEdit.setText('{:.3f}'.format(self.model.bundle_density))
