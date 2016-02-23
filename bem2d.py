#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
from input_parameters import PARAMETERS as P
from terminal_output import print_output as po

po().prog_title('1.0.100115a')

if P['SW_FSI']:
    if P['SW_SPRING']:
        # Run the FSI solver with leading edge spring model
        execfile("spring_bem2d.py")
    else:
        # Run the Fluid Structure Interaction Solver with the BEM Solver
        execfile("FSI_bem2d.py")
else:
    # Run the Boundary Element Solver
    execfile("rigid_bem2d.py")