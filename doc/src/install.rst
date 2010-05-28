Installation
============

Normally, you will need to satisfy some prerequisities before you can compile and install Hermes3D.
These are:

* cmake
* make
* C++ compiler
* Fortran compiler
* Judy
* BLAS
* UMFPACK

Having these packages installed, you can do:

.. code-block:: bash

   $ cmake . # as yourself
   $ make
   $ make install # as root, if your install path is system-wide

In order to install Hermes3D into different location, you need to change ``CMAKE_INSTALL_PREFIX``
variable.


Customization
-------------

Hermes3D is quite modular and can be build with several options. Configuration can be done via
``CMake.vars`` file that has to be placed in the build directory. Look into ``CMake.vars.example``
file for the file format.

- Element types

  * ``WITH_TETRA`` -- enable/dislable teterahedral elements.
  * ``WITH_HEX`` -- enable/dislable hexahedral elements.

- Library type:

  * ``REAL`` -- build real version of the library (scalar type is ``double``)
  * ``COMPLEX`` -- build complex version of the library (scalar type is ``std::complex<double>``)

- Modules

  * ``WITH_UMFPACK`` -- build with support for UMFPACK solver.
  * ``WITH_PETSC`` -- build with support for PETSc solver.
  * ``WITH_PARDISO`` -- build with support for PARDISO solver.
  * ``WITH_MUMPS`` -- build with support for MUMPS solver.
  * ``WITH_TRILINOS`` -- build with Trilinos support (AztecOO, Epetra, ML, Ifpack, NOX, Amesos).
  * ``WITH_METIS`` -- build with METIS support.
  * ``WITH_OPENMP`` -- build with OpenMP support.
  * ``WITH_MPI`` -- build with MPI support (currently no effect).

  If you have problems with CMake not finding your packages, you might want to check 
  ``cmake/FindXYZ.cmake`` files for further details and configuration options.

- Misc

  * ``DEBUG`` -- build debugging version of the library.
  * ``DEBUG_ORDER`` -- use the maximal integration order for integral evaluation.
  * ``WITH_TESTS`` -- build the tests to check that Hermes3D is doing what it is supposed to.
  * ``DEV_TESTS`` -- build developers tests. It is not recommended for normal users, these tests
    take very long time to finish (approx. weeks)
  * ``ADDITIONAL_LIBS`` -- the list of additional libraries that you need to link the binary files
    against in case you have some more requirements to fulfill. For example, if your PETSc is
    compiled with X11 support you need to link against X11 libs. Then you list all the libraries
    here.
  * ``OUTPUT_DIR`` -- set this to a directory that will contain the output files like solutions,
    convergence graphs, etc. Setting this to ``"."`` will write these files into the current
    directory.

Another way how to configure the build process is to use ``ccmake``. Then the most of the options
listed above can be set via textual UI. The paths where your packages sit have to be specified in
``CMake.vars`` the same way as mentioned above.

Notes
-----

* To build documentation, you will need to install the following packages:

   - Doxygen
   - breathe
   - sphinx (0.6.1 works)
   - TeX
   - dvipng

  ``make doc`` builds the documentation in html format; there is also ``make doc-tex`` which builds
  the TeX files with the documentation.

* When building the complex version of Hermes3D with PETSc support, you will need PETSc build with
  C++ support (i.e. ``--with-clanguage=C++`` when building PETSc)

* **Warning:** if you try to build the release version of Hermes3D on 64 bit machine using gcc 4.4.x
  you will get an error. This seems to be a problem in gcc, not in Hermes3D. To workaround this
  issue, build the debug version, if you cannot use different compiler or different version of gcc.   
 
* **Trilinos** can be build using our howto_.


.. _howto: http://hpfem.org/main/howto/howto.html