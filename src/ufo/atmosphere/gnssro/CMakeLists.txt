##########################
add_subdirectory( RefGSI )
add_subdirectory( BndGSI )
add_subdirectory( BndROPP1D )

PREPEND( _p_refgsi_files     "RefGSI"    ${refgsi_src_files} )
PREPEND( _p_bndgsi_files     "BndGSI"    ${bndgsi_src_files} )
PREPEND( _p_bndropp1d_files  "BndROPP1D" ${bndropp1d_src_files} )

set ( gnssro_files
        utils/gnssro_mod_constants.F90
        utils/gnssro_mod_transform.F90
        utils/gnssro_mod_grids.F90
        utils/lag_interp.F90
        ${_p_refgsi_files}
        ${_p_bndgsi_files}
        ${_p_bndropp1d_files}
    PARENT_SCOPE  
)
