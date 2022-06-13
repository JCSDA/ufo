# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

add_subdirectory( ZenithTotalDelayMetOffice )
add_subdirectory( ZenithTotalDelayROPP )

PREPEND( _p_metoffice_files  "operators/groundgnss/ZenithTotalDelayMetOffice" ${groundgnssmetoffice_src_files}  )
PREPEND( _p_groundgnssropp_files  "operators/groundgnss/ZenithTotalDelayROPP" ${groundgnssropp_src_files} )

set ( groundgnss_src_files
      ${_p_metoffice_files}
      ${_p_groundgnssropp_files}
      PARENT_SCOPE
)
