#!/usr/bin/env python
"""
find-spcgrp: Determine the space group of a crystal structure
"""

from __future__ import print_function
import numpy as np
from pyspglib import spglib as sp
from ase import Atoms
from ase.io import vasp
import sys
import string as st
from ase.visualize import view

class msegui:

  def __init__(self):

      pass

  def run(self):

      atoms=Atoms()
      view(atoms)

if __name__ == "__main__":
  # stand-alone mode
  import argparse
  parser = argparse.ArgumentParser(
    description='MSE GUI program')
  parser.add_argument('--reportlevel', type=int, default=4,
                      help='Report level index (0-5)')
  # parser.add_argument('infile', help='input file', default=None)

  args = parser.parse_args()

  gui = msegui()
  stat = gui.run()
