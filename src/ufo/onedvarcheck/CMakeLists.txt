# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( onedvarcheck_files
    RTTOVOneDVarCheck.h
    RTTOVOneDVarCheck.cc
    RTTOVOneDVarCheck.interface.F90
    RTTOVOneDVarCheck.interface.h
    ufo_rttovonedvarcheck_mod.F90
    ufo_rttovonedvarcheck_utils_mod.F90
    ufo_rttovonedvarcheck_init_mod.F90
    ufo_rttovonedvarcheck_process_mod.F90
    ufo_rttovonedvarcheck_minimize_newton_mod.F90
    ufo_rttovonedvarcheck_minimize_ml_mod.F90
    ufo_rttovonedvarcheck_forward_model_mod.F90
)

PREPEND( _p_onedvarcheck_files        "onedvarcheck"                     ${onedvarcheck_files} )

set ( onedvarcheck_src_files
      ${_p_onedvarcheck_files}
      PARENT_SCOPE
)

