Gd2C view settings
==================

* Boundary: Objects >> Boundary
  * ranges: (-1:3, -1:3, -1:3)

* Cutoff planes: Objects >> boundary >> cutoff planes
  * (1 1 1) at 4d     top
  * (-1 -1 -1) at -1d bottom
  * (-1 1 0) at 1.2d  back
  * (1 -1 0) at 0.6d  front

  * (1 0 -1) at 2d
    
  * (1 1 -2) at 4d    right
  * (-1 -1 2) at 4d   left

* Objects >> Orientation

  * Project along the normal to (hkl)

  * View direction

    * projection vector = (1,-1, 0)
    * upward vector = (1,1,1)

* Slice plane: Edit >> Lattice Planes

  * (1,-1, 0) at 0d back
  * (1, 1, -2) at 4d bottom
  * (-1 -1 -1) at -2d right

* Set value range: Objects >> Properties >> Sections
  * Max = 0.006, Min = -0.006

* Bonds: Edit >> Bonds

  o Gd-C 3.0


* Isosurface: Objects >> Properties >> Isosurface
  * level= 0.0045
