# (C) Copyright 2017-2020 Met Office.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( gnssroonedvarcheck_files
    GNSSROOneDVarCheck.h
    GNSSROOneDVarCheck.cc
    GNSSROOneDVarCheck.interface.F90
    GNSSROOneDVarCheck.interface.h
    GNSSROOneDVarCheckParameters.h
    ufo_gnssroonedvarcheck_do1dvar_mod.f90
    ufo_gnssroonedvarcheck_eval_derivs_mod.f90
    ufo_gnssroonedvarcheck_get_bmatrix_mod.f90
    ufo_gnssroonedvarcheck_humidcheck_mod.f90
    ufo_gnssroonedvarcheck_mod.f90
    ufo_gnssroonedvarcheck_pen_mod.f90
    ufo_gnssroonedvarcheck_rootsolv_mod.f90
    ufo_gnssroonedvarcheck_setom1_mod.f90
    ufo_gnssroonedvarcheck_utils_mod.f90
)

PREPEND( _p_gnssroonedvarcheck_files        "gnssroonedvarcheck"        ${gnssroonedvarcheck_files} )

set ( gnssroonedvarcheck_src_files
      ${_p_gnssroonedvarcheck_files}
      PARENT_SCOPE
)

