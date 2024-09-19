#!/usr/bin/env python
"""
vasp2msf: Convert vasp file to msf format
"""

import sys
from xml.etree import ElementTree as ET
import os
from qmutil import vasp, msf

class vasp2msf:

  def __init__(self, reportlevel=1, coord="fraction"):

    self.reportlevel = reportlevel
    self.coord = coord

  def run(self, infname, outfname):

    try:
      infile = open(infname, 'r')
      print '%vasp2msf: Reading the input file "' + infname + '"'
    except IOError:
      print '%vasp2msf: Failed to open the input file "' + infname + '"'
      return False

    self.mol = vasp.read(infile, self.reportlevel)
    infile.close()
    self.mol.finalize()
    self.mol.id = 0
    if self.reportlevel >= 3:
      print "%vasp2msf: Resultant molecule"
      self.mol.report()

    try:
      outfile = open(outfname, 'w')
      print '%vasp2msf: Writing output file "' + outfname + '"'
    except IOError:
      print '%vasp2msf: Failed to open an output file "' + outfname + '"'
      return False

    msf.write(self.mol, outfile, self.coord, self.reportlevel)
    outfile.close()

    # done
    print "%vasp2msf: END"
    return True

if __name__ == "__main__":
  # stand-alone mode
  if len(sys.argv) < 3:
    print "%vasp2msf: Usage: vasp2msf.py infile outfile"
  else:
    app = vasp2msf(reportlevel=4)
    stat = app.run(sys.argv[1], sys.argv[2])
