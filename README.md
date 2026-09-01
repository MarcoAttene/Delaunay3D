# Delaunay3D - A generic data structure to represent and compute 3D Delaunay tetrahedrizations
This code implements a templated C++ class to represent Delaunay tetrahedral meshes and to compute them starting from points or from sets of indexed tetrahedra. The class is templated on the underlying type used to represent 3D points, and the code provides a basic type that can be readily used.
The underlying data structure is described in "**Constrained Delaunay Tetrahedrization: A robust and practical approach**" by L. Diazzi, D. Panozzo, A. Vaxman and <a href="http://saturno.ge.imati.cnr.it/ima/personal-old/attene/PersonalPage/attene.html">M. Attene</a> (ACM Trans Graphics Vol 42, N. 6, Procs of SIGGRAPH Asia 2023). 
You may download a copy here: http://arxiv.org/abs/2309.09805

## Usage
Clone this repository with:
```
git clone https://github.com/MarcoAttene/Delaunay3D
```

Once done, you may build an example executable as follows:
```
cmake -B build -S .
```

This will produce an appropriate building configuration for your system.
On Windows MSVC, this will produce a cdt.sln file.
On Linux/MacOS, this will produce a Makefile. 
Use it as usual to compile cdt. Alternatively, you can use the command line:
```
cmake --build build --config Release
```

When compiled, the example code generates an executable called ``delaunay3d``.
Launch it with no command line parameters to have a list of supported options.

To use within your program, just include "delaunay.h" and configure your build
with the same options included in the CMakeLists provided.

## License
This program is distributed under the terms of either the GNU GPL or the GNU LGPL license.
The code can be compiled in two ways, depending on how CMake is invoked.
If you build using ``CMake -DLGPL=ON ..``, you may choose between GPL and LGPL at your option.
If you build using ``CMake -DLGPL=OFF ..`` or just ``CMake ..``, the code makes use of modified 
parts of a third-party code which requires you to accept the terms of the GPL license.
See ``src/delaunay.h`` for details.

In either case, the program is distributed in the hope that it will be      
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of   
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.