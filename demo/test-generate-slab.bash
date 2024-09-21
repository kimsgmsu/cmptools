#!/bin/bash

generate-slab.py Al4.fcc.vasp Al-slab.vasp --miller_index="[1,1,1]" --min_slab_size=4.0 --min_vacuum_size=1.0 --in_unit_planes=True --primitive=False --max_normal_search=10 --reorient_lattice=True


