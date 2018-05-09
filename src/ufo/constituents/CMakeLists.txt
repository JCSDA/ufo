add_subdirectory( aod )

PREPEND( _p_aod_files  "constituents/aod"  ${aod_files} )

set ( constituents_src_files
    ${_p_aod_files}
    PARENT_SCOPE
)
