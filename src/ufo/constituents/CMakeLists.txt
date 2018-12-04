# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND constituents_src_files
)

add_subdirectory( aod )

PREPEND( _p_aod_files  "constituents/aod"  ${aod_files} )

set ( constituents_src_files
    ${_p_aod_files}
    PARENT_SCOPE
)
