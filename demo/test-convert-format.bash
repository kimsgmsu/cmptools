#!/bin/bash

convert-format.py test1.vasp tmp1.msf

convert-format.py tmp1.msf new1.vasp

convert-format.py GaAs.vasp tmp2.msf

convert-format.py tmp2.msf new2.vasp --coord=cart

