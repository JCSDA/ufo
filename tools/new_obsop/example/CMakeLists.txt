# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( example_files
    ObsExample.h
    ObsExample.cc
    ObsExampleTLAD.h
    ObsExampleTLAD.cc
    ObsExample.interface.F90
    ObsExample.interface.h
    ObsExampleParameters.h
    ObsExampleTLAD.interface.F90
    ObsExampleTLAD.interface.h
    ufo_example_mod.F90
    ufo_example_tlad_mod.F90
)

PREPEND( _p_example_files     "example"     ${example_files} )

set ( example_src_files
      ${_p_example_files}
      PARENT_SCOPE
)
