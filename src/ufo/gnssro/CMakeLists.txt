##########################
add_subdirectory( RefGSI )
add_subdirectory( BndGSI )
add_subdirectory( BndROPP1D )
add_subdirectory( BndROPP2D )
add_subdirectory( utils )
add_subdirectory( QC )

PREPEND( _p_refgsi_files     "gnssro/RefGSI"    ${refgsi_src_files} )
PREPEND( _p_bndgsi_files     "gnssro/BndGSI"    ${bndgsi_src_files} )
PREPEND( _p_bndropp1d_files  "gnssro/BndROPP1D" ${bndropp1d_src_files} )
PREPEND( _p_bndropp2d_files  "gnssro/BndROPP2D" ${bndropp2d_src_files} )
PREPEND( _p_utils_files      "gnssro/utils"     ${utils_src_files} )
PREPEND( _p_qc_files         "gnssro/QC"        ${qc_src_files} )

set ( gnssro_src_files
        ${_p_refgsi_files}
        ${_p_bndgsi_files}
        ${_p_bndropp1d_files}
        ${_p_bndropp2d_files}
        ${_p_utils_files}
        ${_p_qc_files}
    PARENT_SCOPE  
)
