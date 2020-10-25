.. Building UFO in OS X documentation master file, created by
   sphinx-quickstart on Fri Jan 19 15:58:54 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Building UFO in OS X
====================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   config_mac
   config_ec_tools

Purpose
=======

This document is meant to show the user how to prepare their OS X operating
system (iMac, MacBook, etc.) for compiling UFO on their system, and then how
to run the build of the UFO code.

.. index:: Background

Background
==========

The Unified Forward Operator (UFO) is the software that implements the "observation operators"
piece of the data assimilation cost function.

The data assimilation cost function is:

.. math::
  J(x_0) = \left\{\sum_{i=0}^{N} [H(x_i) - y_i]^T R^{-1} [H(x_i) -y_i]\right\} + (x_o - x_b)^T B^{-1} (x_o - X_b)

where the observation operatores (UFO) are defined by :math:`H(x)`.

.. index:: Prepare the Mac with third party software

Prepare the Mac with third party software
=========================================

In preparation for building UFO, the Mac needs to be configured properly with
third party software and the EC build tools.
The tool of choice for installing third party software packages on the Mac
is `Homebrew <https://brew.sh/>`_.
Homebrew will be used to install `GNU compilers <https://gcc.gnu.org/>`_,
and the `Eigen <http://eigen.tuxfamily.org/>`_ 
and `Boost <http://www.boost.org/>`_ libraries.

This is accomplished in two steps:

#. :doc:`Install third party software on the Mac </config_mac>`
#. :doc:`Build the EC tools and libraries </config_ec_tools>`.

.. index:: Build UFO components

Build UFO components
====================

Upon completion of the Mac configuration, you should be in good shape to
build UFO.
UFO has three components:

#. CRTM (Community Radiation Transfer Model)
#. OOPS (Object Oriented Prgramming System)
#. IODA (Interface for Observational Data Assimilation)

Each component can be built with the same sequence of commands with only
a variation on the Git branch that is used.
The sequence for CRTM will be shown, and the branches for all three
components will be noted.

To build CRTM:

.. code:: bash

  cd $HOME/projects
  git clone https://github.com/UCAR/crtm.git
  git checkout master
  cd crtm
  mkdir build
  cd build
  ecbuild $HOME/projects/crtm    # check ecbuild.log for errors/warnings
  make -j4
  ctest

.. note::
  There is no "make install" step necessary for this process.
  This is so because UFO only needs the libraries created by building
  CRTM, OOPS, and IODA.
  When building UFO, the ecbuild command will be told where to look
  for the component libraries (see below).

Repeat the above command sequence to build OOPS (substitute the string
"oops" for "crtm") and IODA (substitute the string "ioda" for "crtm").
Make sure that the proper Git branch is used for each component and noted
below.

.. topic:: Git branches for UFO components

  For OOPS: git checkout feature/ufo

  For IODA: git checkout develop

.. index:: Build UFO

Build UFO 
=========

Once CRTM, OOPS and IODA are successfully built (tests pass) UFO can
then be compiled.
This is accomplished with a similar command sequence to that
which was used for the components.

.. code:: bash

  cd $HOME/projects
  git clone https://github.com/UCAR/ufo.git
  git checkout develop
  cd ufo
  mkdir build
  cd build
  ecbuild -Dcrtm_SOURCE_DIR=$HOME/projects/crtm -DIODA_PATH=$HOME/projects/ioda/build -DCRTM_PATH=$HOME/projects/crtm/build -DOOPS_PATH=$HOME/projects/oops/build ~/projects/ufo
  make -j4
  ctest

.. note::
  The long string of -D options on the ecbuild command are there to tell
  the UFO build process where the libraries for CRTM, OOPS and IODA are
  located (the "build" directories for each component).
  An additional option, -Dcrtm_SOURCE_DIR, is needed to tell the UFO build
  process where test data from CRTM is located.


.. index:: jedi-bundle

jedi-bundle
===========

You may have noticed some awkwardness in the UFO build process.
Specifically, the need to remember the proper Git branches for the
checkout commands, and the need to specify the locations of the
component build directories greatly increase the possibility of
making mistakes along the way.
While it is instructive to go through the above process, a separate
UCAR GitHub repository (called "jedi-bundle") has been created
to mitigate the opportunity for errors during the build process.
Jedi-bundle takes care of checking out the proper git branches, building
all of the components (along with UFO) and passing on the locations
of the components to the UFO build process.
Jedi-bundle can be downloaded and built by running the following:

.. code:: bash

  cd $HOME/projects
  git clone https://github.com/UCAR/jedi-bundle.git
  cd jedi-bundle
  mkdir build
  cd build
  ecbuild $HOME/projects/jedi-bundle
  make -j4
  ctest

Much easier!


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
