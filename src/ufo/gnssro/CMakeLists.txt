##########################
add_subdirectory( RefNBAM )
add_subdirectory( RefMetOffice )
add_subdirectory( BndNBAM )
add_subdirectory( BndROPP1D )
add_subdirectory( BndROPP2D )
add_subdirectory( utils )
add_subdirectory( QC )
add_subdirectory( BendMetOffice )

PREPEND( _p_refnbam_files      "gnssro/RefNBAM"       ${refnbam_src_files} )
PREPEND( _p_refmetoffice_files "gnssro/RefMetOffice"  ${refmetoffice_src_files} )
PREPEND( _p_bndnbam_files      "gnssro/BndNBAM"       ${bndnbam_src_files} )
PREPEND( _p_bndropp1d_files    "gnssro/BndROPP1D"     ${bndropp1d_src_files} )
PREPEND( _p_bndropp2d_files    "gnssro/BndROPP2D"     ${bndropp2d_src_files} )
PREPEND( _p_utils_files        "gnssro/utils"         ${utils_src_files} )
PREPEND( _p_qc_files           "gnssro/QC"            ${qc_src_files} )
PREPEND( _p_metoffice_files    "gnssro/BendMetOffice" ${bendmetoffice_src_files} )

set ( gnssro_src_files
        ${_p_refnbam_files}
        ${_p_refmetoffice_files}
        ${_p_bndnbam_files}
        ${_p_bndropp1d_files}
        ${_p_bndropp2d_files}
        ${_p_utils_files}
        ${_p_qc_files}
        ${_p_metoffice_files}
    PARENT_SCOPE  
)
