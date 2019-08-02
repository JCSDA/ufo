# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

add_subdirectory( Correction )

PREPEND( _p_surface_files       "surface"       ${surface_files} )
PREPEND( _p_correction_files         "surface/Correction"        ${correction_src_files} )

set ( surface_src_files
      ${_p_surface_files}
      ${_p_correction_files}
      PARENT_SCOPE
)
