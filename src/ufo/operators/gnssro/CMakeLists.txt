##########################
add_subdirectory( RefNCEP )
add_subdirectory( RefMetOffice )
add_subdirectory( BndNBAM )
add_subdirectory( BndROPP1D )
add_subdirectory( BndROPP2D )
add_subdirectory( utils )
add_subdirectory( QC )
add_subdirectory( BendMetOffice )

PREPEND( _p_refncep_files      "operators/gnssro/RefNCEP"       ${refncep_src_files} )
PREPEND( _p_refmetoffice_files "operators/gnssro/RefMetOffice"  ${refmetoffice_src_files} )
PREPEND( _p_bndnbam_files      "operators/gnssro/BndNBAM"       ${bndnbam_src_files} )
PREPEND( _p_bndropp1d_files    "operators/gnssro/BndROPP1D"     ${bndropp1d_src_files} )
PREPEND( _p_bndropp2d_files    "operators/gnssro/BndROPP2D"     ${bndropp2d_src_files} )
PREPEND( _p_utils_files        "operators/gnssro/utils"         ${utils_src_files} )
PREPEND( _p_qc_files           "operators/gnssro/QC"            ${qc_src_files} )
PREPEND( _p_metoffice_files    "operators/gnssro/BendMetOffice" ${bendmetoffice_src_files} )

set ( gnssro_src_files
        ${_p_refncep_files}
        ${_p_refmetoffice_files}
        ${_p_bndnbam_files}
        ${_p_bndropp1d_files}
        ${_p_bndropp2d_files}
        ${_p_utils_files}
        ${_p_qc_files}
        ${_p_metoffice_files}
    PARENT_SCOPE  
)
