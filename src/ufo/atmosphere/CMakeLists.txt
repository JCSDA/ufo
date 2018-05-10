add_subdirectory( utils )
add_subdirectory( radiosonde )
add_subdirectory( radiance )

PREPEND( _p_utils_files       "atmosphere/utils"       ${utils_files} )
PREPEND( _p_radiosonde_files  "atmosphere/radiosonde"  ${radiosonde_files} )
PREPEND( _p_radiance_files    "atmosphere/radiance"    ${radiance_files} )

set ( atmosphere_src_files
    ${_p_utils_files}
    ${_p_radiosonde_files}
    ${_p_radiance_files}
    PARENT_SCOPE
)
