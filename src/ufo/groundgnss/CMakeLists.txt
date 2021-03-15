# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

add_subdirectory( ZenithTotalDelayMetOffice )
add_subdirectory( RefROPP1D )

PREPEND( _p_metoffice_files  "groundgnss/ZenithTotalDelayMetOffice" ${groundgnssmetoffice_src_files}  )
PREPEND( _p_refropp1d_files  "groundgnss/RefROPP1D" ${refropp1d_src_files} )

set ( groundgnss_src_files
      ${_p_refropp1d_files}
      ${_p_metoffice_files}
      PARENT_SCOPE
)
