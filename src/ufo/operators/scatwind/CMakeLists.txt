##########################
add_subdirectory( NeutralMetOffice )

PREPEND( _p_metoffice_files  "operators/scatwind/NeutralMetOffice" ${neutralmetoffice_src_files} )

set ( scatwind_src_files
        ${_p_metoffice_files}
    PARENT_SCOPE  
)
