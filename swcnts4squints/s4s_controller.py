# -*- coding: utf-8 -*-
"""
=================================================================
swcnts4squints controller (:mod:`swcnts4squints.s4s_controller`)
=================================================================

.. currentmodule:: swcnts4squints.s4s_controller

"""
from __future__ import absolute_import, print_function, division
from __future__ import unicode_literals
from builtins import object

from PyQt4.QtGui import QApplication

from .s4s_view import S4SView

__all__ = ['S4SController']


class S4SController(object):
    """swcnts4squints MVC controller class.

    Parameters
    ----------
    args : `<python:sys.argv>`
    model : :py:class:`swcnts4squints.s4s_model` instance

    """
    def __init__(self, args, model=None):
        app = QApplication(args)
        self.model = model

        self.view = S4SView(self, self.model)
        self.model.init()
        self.view.show()
        app.exec_()

    def refresh_view(self):
        """refresh view."""
        self.view.update_app_view()
