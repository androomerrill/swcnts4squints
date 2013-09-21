#!/usr/bin/env python

import sys

from s4s_controller import S4SController
from s4s_model import S4SModel

class S4S(object):

    def __init__(self, args):
        model = S4SModel()
        S4SController(args, model)

def main():
    S4S(sys.argv)

if __name__ == "__main__":
    sys.exit(main())
