add_subdirectory( utils )
add_subdirectory( radiosonde )
add_subdirectory( radiance )
add_subdirectory( aircraft )

PREPEND( _p_utils_files       "atmosphere/utils"       ${utils_files} )
PREPEND( _p_radiosonde_files  "atmosphere/radiosonde"  ${radiosonde_files} )
PREPEND( _p_radiance_files    "atmosphere/radiance"    ${radiance_files} )
PREPEND( _p_aircraft_files    "atmosphere/aircraft"    ${aircraft_files} )

set ( atmosphere_src_files
    ${_p_utils_files}
    ${_p_radiosonde_files}
    ${_p_radiance_files}
    ${_p_aircraft_files}
    PARENT_SCOPE
)
