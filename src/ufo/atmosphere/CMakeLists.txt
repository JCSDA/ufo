# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND atmosphere_src_files
    FortranAtmosphere.h
)
add_subdirectory( radiosonde )
add_subdirectory( radiance )
add_subdirectory( aircraft )
add_subdirectory( gnssro )

PREPEND( _p_radiosonde_files  "atmosphere/radiosonde"  ${radiosonde_files} )
PREPEND( _p_radiance_files    "atmosphere/radiance"    ${radiance_files} )
PREPEND( _p_aircraft_files    "atmosphere/aircraft"    ${aircraft_files} )
PREPEND( _p_gnssro_files      "atmosphere/gnssro"      ${gnssro_files}  ) 

set ( atmosphere_src_files
    ${_p_radiosonde_files}
    ${_p_radiance_files}
    ${_p_aircraft_files}
    ${_p_gnssro_files} 
    PARENT_SCOPE
)
