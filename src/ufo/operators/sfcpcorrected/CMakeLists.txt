# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( sfcpcorrected_files
    ObsSfcPCorrected.h
    ObsSfcPCorrected.cc
    ObsSfcPCorrectedParameters.cc
    ObsSfcPCorrectedParameters.h
    ObsSfcPCorrected.interface.F90
    ObsSfcPCorrected.interface.h
    ufo_sfcpcorrected_mod.F90
)

PREPEND( _p_sfcpcorrected_files     "operators/sfcpcorrected"     ${sfcpcorrected_files} )

set ( sfcpcorrected_src_files
      ${_p_sfcpcorrected_files}
      PARENT_SCOPE
)
