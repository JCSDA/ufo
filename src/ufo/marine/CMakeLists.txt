add_subdirectory( gsw )
add_subdirectory( obsop )
add_subdirectory( seaicethickness )
add_subdirectory( seaicefraction )
add_subdirectory( stericheight )
add_subdirectory( insitutemperature )

PREPEND( _p_gsw_files              "marine/gsw"              ${gsw_files} )
PREPEND( _p_obsop_files            "marine/obsop"            ${obsop_files} )
PREPEND( _p_seaicethickness_files  "marine/seaicethickness"  ${seaicethickness_files} )
PREPEND( _p_seaicefraction_files   "marine/seaicefraction"   ${seaicefraction_files} )
PREPEND( _p_stericheight_files     "marine/stericheight"     ${stericheight_files} )
PREPEND( _p_insitutemperature_files     "marine/insitutemperature"     ${insitutemperature_files} )

set ( marine_src_files
    ${_p_gsw_files}
    ${_p_obsop_files}
    ${_p_seaicethickness_files}
    ${_p_seaicefraction_files}
    ${_p_stericheight_files}
    ${_P_insitutemperature_files}    
    PARENT_SCOPE
)
