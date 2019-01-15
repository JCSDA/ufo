##########################
add_subdirectory( RefGSI )
add_subdirectory( BndGSI )
add_subdirectory( BndROPP1D )
add_subdirectory( utils )
add_subdirectory( QC )

PREPEND( _p_refgsi_files     "RefGSI"    ${refgsi_src_files} )
PREPEND( _p_bndgsi_files     "BndGSI"    ${bndgsi_src_files} )
PREPEND( _p_bndropp1d_files  "BndROPP1D" ${bndropp1d_src_files} )
PREPEND( _p_utils_files      "utils"     ${utils_src_files} )
PREPEND( _p_qc_files         "QC"        ${qc_src_files} )

set ( gnssro_files
        ${_p_refgsi_files}
        ${_p_bndgsi_files}
        ${_p_bndropp1d_files}
        ${_p_utils_files}
        ${_p_qc_files}
    PARENT_SCOPE  
)
