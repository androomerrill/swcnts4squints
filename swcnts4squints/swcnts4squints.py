#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
=============================================================
swcnts4squints CLI (:mod:`swcnts4squints.swcnts4squints`)
=============================================================

.. currentmodule:: swcnts4squints.swcnts4squints

"""
from __future__ import absolute_import, print_function, division

import sys

from .s4s_controller import S4SController
from .s4s_model import S4SModel

__all__ = ['S4S']


class S4S(object):
    """Base class for instantiating the S4S MVC."""
    def __init__(self, args):
        model = S4SModel()
        S4SController(args, model)


def main():
    S4S(sys.argv)

if __name__ == "__main__":
    sys.exit(main())
