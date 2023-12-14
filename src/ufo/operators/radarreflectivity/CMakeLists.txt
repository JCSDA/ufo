##########################
add_subdirectory( directZDA )

PREPEND( _p_directZDA_files     "operators/radarreflectivity/directZDA"     ${directZDA_src_files} )

set ( radarreflectivity_src_files
        ${_p_directZDA_files}
    PARENT_SCOPE  
)
