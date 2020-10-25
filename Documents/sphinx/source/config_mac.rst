Installing third party software on the Mac for building UFO
===========================================================

Many installer tools for OS X (MacBook, iMac, etc.) exist.
MacPorts and homebrew are two popular examples.
Homebrew has been succussfully used for installing tools on the Mac for the purpose
of compiling software (e.g., UFO) in the JEDI project.
For this reason, it is recommended to use homebrew for configuring your Mac.

By default, homebrew installs packages into /usr/local.
It is recommended to stick with the /usr/local install location since many creators
of homebrew packages are used to using /usr/local as their target install location.
It is recommended to change ownership of the subdirectories of /usr/local
to your user id.
Once this is done, you will be able to manage installations with homebrew without
the need for superuser privilages.

.. index:: Install homebrew

To install homebrew itself on your Mac, run the following Ruby command (bash syntax is
being shown in the example below).

.. code:: bash

  /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"


.. note::
  Homebrew is Ruby based and the installation scripts (called "formulas")
  are Ruby scripts.

Once installed, the primary command for the user is "brew".
See the homebrew `website <https:/brew.sh/>`_ for details, or enter
the command "man brew" for even more details.

.. index:: Install gcc 7

The first step is to install the gcc version 7 compilers package.
This is accomplished by running:

.. code:: bash

  brew install gcc

This results in the gcc tools being installed in /usr/local/bin, but they are
suffixed with "-7" so that they do not conflict with existing installed compilers.
However, the existing compilers (likely the clang versions in /usr/bin) are not
compatible with some of the libraries used by UFO.
To work around this, create the following links in /usr/local/bin:

.. code:: bash

  cd /usr/local/bin
  ln -s c++-7 c++
  ln -s gcc-7 cc
  ln -s cpp-7 cpp
  ln -s g++-7 g++
  ln -s gcc-7 gcc

And make sure that "/usr/local/bin" is first in your PATH environment variable.

.. index:: Install Eigen and Boost

UFO uses the third party libraries `Eigen <http://eigen.tuxfamily.org/>`_ 
and `Boost <http://www.boost.org/>`_.
Eigen provides C++ templates for linear algebra, and Boost provides a wide range
of C++ libraries that have been peer reviewed.
UFO is utilizing Boost libraries that assist with unit testing.
To load Eigen and Boost, run the following:

.. code:: bash

  brew install eigen
  brew install --cc=gcc-7 boost

.. note::
  The --cc=gcc-7 option for the boost install is utilized so that the boost
  libraries will be installed using the c++11 code standard.
  This is necessary to get the proper function name mangling in the boost link
  libraries for the downstream UFO link steps.

.. index:: Install CMake

The UFO build system is based on `CMake <https://cmake.org/>`_.
To install CMake, run:

.. code:: bash

  brew install cmake

.. index:: Install HDF5, netCDF, MPI

UFO used `HDF5 <https://www.hdfgroup.org/>`_,
`netCDF <https://www.unidata.ucar.edu/software/netcdf/>`_ and
MPI (e.g., `mpich <https://www.mpich.org/>`_,
`openMPI <https://www.open-mpi.org/>`_) libraries.
To install these, run:

.. code:: bash

  brew install hdf5
  brew install netcdf
  brew install mpich

.. index:: Install Python

When testing the process of compiling UFO on the Mac, it was found that the
`Anaconda <https://www.anaconda.com/>`_ installation of python was not compatible.
It may be useful to fix the Anaconda incompatibility, but as of this time a
solution has not been found.
Homebrew can install python and it is recommended to install version 3 which
can be accomplished by running:

.. code:: bash

  brew install python3

.. index:: Install optional packages

Additional packages, which are optional, that you might find handy are as follows:

.. code:: bash

  brew install tkdiff    # a nice GUI for diff
  brew install graphviz  # for creating diagrams, flow charts, etc.
  brew install macvim    # a nice GUI for vim ("vi improved")
  brew install jupyter   # notebooks
  brew install iterm     # nice terminal emulator (this command loads iTerm2)


