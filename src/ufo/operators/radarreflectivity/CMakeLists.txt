##########################
add_subdirectory( directZDA )
add_subdirectory( genericreflectivity )

PREPEND( _p_directZDA_files     "operators/radarreflectivity/directZDA"     ${directZDA_src_files} )
PREPEND( _p_genericreflectivity_files     "operators/radarreflectivity/genericreflectivity"     ${genericreflectivity_src_files} )

set ( radarreflectivity_src_files
        ${_p_directZDA_files}
        ${_p_genericreflectivity_files}
    PARENT_SCOPE  
)
