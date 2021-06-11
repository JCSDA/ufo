# (C) Copyright 2021 UK Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( categoricaloper_files
    ObsCategorical.h
    ObsCategorical.cc
    ObsCategoricalParameters.h
)
PREPEND( _p_categoricaloper_files     "categoricaloper"     ${categoricaloper_files} )

set ( categoricaloper_src_files
      ${_p_categoricaloper_files}
      PARENT_SCOPE
)
