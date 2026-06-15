# (C) Copyright 2017-2025 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

##########################
add_subdirectory( ppro )

PREPEND( _p_ppro_files     "operators/radarpolarimetric/ppro"     ${ppro_src_files} )

set ( radarpolarimetric_src_files
	${_p_ppro_files}
    PARENT_SCOPE  
)

