# (C) Copyright 2017-2020 Met Office.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( rttovonedvarcheck_files
    RTTOVOneDVarCheck.h
    RTTOVOneDVarCheck.cc
    RTTOVOneDVarCheck.interface.F90
    RTTOVOneDVarCheck.interface.h
    ufo_rttovonedvarcheck_bmatrix_mod.f90
    ufo_rttovonedvarcheck_constants_mod.f90
    ufo_rttovonedvarcheck_minimize_jacobian_mod.f90
    ufo_rttovonedvarcheck_minimize_ml_mod.f90
    ufo_rttovonedvarcheck_minimize_newton_mod.f90
    ufo_rttovonedvarcheck_minimize_utils_mod.f90
    ufo_rttovonedvarcheck_mod.f90
    ufo_rttovonedvarcheck_pcemis_mod.f90
    ufo_rttovonedvarcheck_profindex_mod.f90
    ufo_rttovonedvarcheck_rmatrix_mod.f90
    ufo_rttovonedvarcheck_rsubmatrix_mod.f90
    ufo_rttovonedvarcheck_ob_mod.f90
    ufo_rttovonedvarcheck_obs_mod.f90
    ufo_rttovonedvarcheck_utils_mod.f90
)

PREPEND( _p_rttovonedvarcheck_files        "rttovonedvarcheck"        ${rttovonedvarcheck_files} )

set ( rttovonedvarcheck_src_files
      ${_p_rttovonedvarcheck_files}
      PARENT_SCOPE
)

