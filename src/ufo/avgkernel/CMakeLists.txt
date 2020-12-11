# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( avgkernel_files
    ObsAvgKernel.h
    ObsAvgKernel.cc
    ObsAvgKernel.interface.F90
    ObsAvgKernel.interface.h
    ufo_avgkernel_mod.F90
    ufo_satcolumn_mod.F90
    ObsAvgKernelTLAD.h
    ObsAvgKernelTLAD.cc
    ObsAvgKernelTLAD.interface.F90
    ObsAvgKernelTLAD.interface.h
    ufo_avgkernel_tlad_mod.F90
)

PREPEND( _p_avgkernel_files     "avgkernel"     ${avgkernel_files} )

set ( avgkernel_src_files
      ${_p_avgkernel_files}
      PARENT_SCOPE
)
