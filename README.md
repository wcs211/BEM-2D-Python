# BEM-2D-Python
2D Python boundary element method solver

This is a boundary element solver library devloped and maintained by the Biofluids Research Group, Department of Mechanical Engineering and Mechanics, P.C. Rossin College of Engineering and Applied Science, Lehigh University.

Other related projects:

* [kmoored/BEM-2D-Matlab](https://github.com/kmoored/BEM-2D-Matlab)
* [wcs211/BEM-3D-Python](https://github.com/wcs211/BEM-3D-Python)

## Get Started
* Edit the simulation input file [input_parameters.py](https://github.com/wcs211/BEM-2D-Python/blob/master/input_parameters.py)
* Verrify a 'movies' directory exists in the execuition path for timestep plots to be saved to.
* Execuite the file 'bem2d.py' at the terminal prompt:
    $ python bem2d.py

## Features

* Modular code structure allows easier implementation of new features
* Multiple body interactions
* Implicit and Explicit Kutta condition enforcement

## Future Features
The following features have planned implementation in the code:

* Vortex Particle wake representation
* Lumped wake representation
* Equations of motion solver
* Boundary layer solver for skin friction estimation
* Quadtree collision detection (fencing scheme)
* Fast Multipole solver
* Parallel processing
* GPGPU processing
