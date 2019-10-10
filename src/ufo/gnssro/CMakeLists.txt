##########################
add_subdirectory( Ref )
add_subdirectory( BndNBAM )
add_subdirectory( BndROPP1D )
add_subdirectory( BndROPP2D )
add_subdirectory( utils )
add_subdirectory( QC )

PREPEND( _p_ref_files        "gnssro/Ref"       ${ref_src_files} )
PREPEND( _p_bndnbam_files    "gnssro/BndNBAM"   ${bndnbam_src_files} )
PREPEND( _p_bndropp1d_files  "gnssro/BndROPP1D" ${bndropp1d_src_files} )
PREPEND( _p_bndropp2d_files  "gnssro/BndROPP2D" ${bndropp2d_src_files} )
PREPEND( _p_utils_files      "gnssro/utils"     ${utils_src_files} )
PREPEND( _p_qc_files         "gnssro/QC"        ${qc_src_files} )

set ( gnssro_src_files
        ${_p_ref_files}
        ${_p_bndnbam_files}
        ${_p_bndropp1d_files}
        ${_p_bndropp2d_files}
        ${_p_utils_files}
        ${_p_qc_files}
    PARENT_SCOPE  
)
