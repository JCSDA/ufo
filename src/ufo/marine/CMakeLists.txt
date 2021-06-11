# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND marine_src_files
    FortranMarine.h
)
add_subdirectory( adt )
add_subdirectory( coolskin )
add_subdirectory( seaicefraction )
add_subdirectory( seaicethickness )
add_subdirectory( chleuzintegr )
add_subdirectory( utils )
if( ${gsw_FOUND} )
  add_subdirectory( insitutemperature )
  add_subdirectory( marinevertinterp )
  add_subdirectory( obsop )
endif( ${gsw_FOUND} )

PREPEND( _p_utils_files               "marine/utils"               ${utils_files} )
PREPEND( _p_adt_files                 "marine/adt"                 ${adt_files} )
PREPEND( _p_coolskin_files            "marine/coolskin"            ${coolskin_files} )
PREPEND( _p_seaicefraction_files      "marine/seaicefraction"      ${seaicefraction_files} )
PREPEND( _p_seaicethickness_files     "marine/seaicethickness"     ${seaicethickness_files} )
PREPEND( _p_chleuzintegr_files        "marine/chleuzintegr"        ${chleuzintegr_files} )
if( ${gsw_FOUND} )
  PREPEND( _p_insitutemperature_files   "marine/insitutemperature"   ${insitutemperature_files} )
  PREPEND( _p_marinevertinterp_files    "marine/marinevertinterp"    ${marinevertinterp_files} )
  PREPEND( _p_obsop_files               "marine/obsop"               ${obsop_files} )
endif( ${gsw_FOUND} )

set ( marine_src_files
    ${_p_utils_files}
    ${_p_adt_files}
    ${_p_coolskin_files}
    ${_p_marinevertinterp_files}
    ${_p_obsop_files}
    ${_p_seaicefraction_files}
    ${_p_seaicethickness_files}
    ${_p_insitutemperature_files}
    ${_p_chleuzintegr_files}
    PARENT_SCOPE
)
