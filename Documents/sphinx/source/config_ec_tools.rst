Configuring the EC tools on the Mac for building UFO
====================================================

After `configuring your Mac </config_mac>`_, the next step is to build the EC
tools which will then be used to build UFO and its components.
The primary tool that you will need is "ecbuild" which
provides a comprehensive set of `CMake <https://cmake.org/>`_ extensions
for software compiliation, test and installation.

The packages for this step live on `GitHub <https://github.com/>`_ in respositories
belonging to user "UCAR".
The general process is to clone the GitHub repository, create a build directory
(i.e., "out-of-source" style of build), run cmake followed by "make", "ctest" and
"make install".
This sequence will download source code, compile the code, run unit tests and 
finally install the newly created binaries and libraries in /usr/local.

.. index:: Build and install ecbuild

The first step is to create the "ecbuild" system which will subsequently be used
to build all other packages.
For the following examples, it will be assumed that the user will build all of their
packages in a subdirectory of their home directory called "projects",
(\$HOME/projects).
Additionally, it is assumed that for each package the user will run the build
process in a subdirectory of that package called "build".
In actuality, the ecbuild/CMake system is quite flexible in regard to where the
user can locate their builds and the only restriction is that the build location
live outside the package source directory.

.. code:: bash

  cd $HOME/projects
  git clone https://github.com/UCAR/ecbuild.git
  git checkout master    # make sure on branch "master"
  cd ecbuild
  mkdir build
  cd build
  cmake $HOME/projects/ecbuild   # Check CMakeFiles/CMake*.log for errors/warnings
  make
  make install

.. note::
  If you have a Mac with multiple cores, the make process can be sped up
  (via parallel processing) by running "make -jN", where N is the number
  of processes to use for running the job.

  The command "make VERBOSE=1" can be run if you want to see the commands that
  the make process is issuing.

.. index:: Build and install eckit and fckit

Once ecbuild is installed, the EC tool build and install process can be completed
by building "eckit" and "fckit".
These two packages are built in an identical fashion, therefore only eckit will be
shown as an example.
Also, it will be assumed that the Mac has 4 cores available for the make process.

.. code:: bash

  cd $HOME/projects
  git clone https://github.com/UCAR/eckit.git
  git checkout master    # make sure on branch "master"
  cd eckit
  mkdir build
  cd build
  ecbuild $HOME/projects/eckit   # Check ecbuild.log for errors/warnings
  make -j4
  ctest
  make install

