# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND generic_src_files
)
add_subdirectory( identity )

PREPEND( _p_identity_files  "generic/identity"  ${identity_files} )

set ( generic_src_files
    ${_p_identity_files}
    PARENT_SCOPE
)
