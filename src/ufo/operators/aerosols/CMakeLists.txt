# (C) Copyright 2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


add_subdirectory( AODExt )
add_subdirectory( AODMetOffice )
if( ${geos-aero_FOUND} )
   add_subdirectory( AODLuts )
endif( ${geos-aero_FOUND} )


PREPEND( _p_aodext_files  "operators/aerosols/AODExt"  ${aodext_src_files} )
PREPEND( _p_aodmetoffice_files  "operators/aerosols/AODMetOffice"  ${aodmetoffice_src_files} )
if( ${geos-aero_FOUND} )
  PREPEND( _p_aodluts_files  "operators/aerosols/AODLuts"  ${aodluts_src_files} )
endif( ${geos-aero_FOUND} )

set ( aerosols_src_files
    ${_p_aodext_files}
    ${_p_aodmetoffice_files}
    ${_p_aodluts_files}
    PARENT_SCOPE
)

