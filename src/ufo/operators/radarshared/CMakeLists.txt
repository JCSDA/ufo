# (C) Copyright 2017-2025 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Shared radar operator utilities

set ( _radarshared_files
    MicrophysicsOptions.h
)

PREPEND( _p_radarshared_files  "operators/radarshared"  ${_radarshared_files} )

set ( radarshared_src_files
    ${_p_radarshared_files}
    PARENT_SCOPE
)

