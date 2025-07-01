# (C) Copyright 2025 Met Office.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( refractivityonedvarcheck_files
    RefractivityOneDVarCheck.h
    RefractivityOneDVarCheck.cc
    RefractivityOneDVarCheck.interface.F90
    RefractivityOneDVarCheck.interface.h
    RefractivityOneDVarCheckParameters.h
    ufo_refractivityonedvarcheck_do1dvar_mod.f90
    ufo_refractivityonedvarcheck_mod.f90
    ufo_refractivityonedvarcheck_rootsolv_mod.f90
    ufo_refractivityonedvarcheck_utils_mod.f90
)

PREPEND( _p_refractivityonedvarcheck_files        "refractivityonedvarcheck"        ${refractivityonedvarcheck_files} )

set ( refractivityonedvarcheck_src_files
      ${_p_refractivityonedvarcheck_files}
      PARENT_SCOPE
)

