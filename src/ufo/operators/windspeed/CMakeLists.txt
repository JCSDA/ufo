# Crown Copyright 2021 Met Office

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( windspeed_files
    ObsWindSpeed.h
    ObsWindSpeed.cc
    ObsWindSpeedParameters.h

)

PREPEND( _p_windspeed_files     "operators/windspeed"     ${windspeed_files} )

set ( windspeed_src_files
      ${_p_windspeed_files}
      PARENT_SCOPE
)
