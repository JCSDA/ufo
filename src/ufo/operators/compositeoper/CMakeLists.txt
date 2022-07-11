# (C) Copyright 2021 UK Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( compositeoper_files
    ObsComposite.h
    ObsComposite.cc
    ObsCompositeParameters.h
    ObsCompositeTLAD.h
    ObsCompositeTLAD.cc
)
PREPEND( _p_compositeoper_files     "operators/compositeoper"     ${compositeoper_files} )

set ( compositeoper_src_files
      ${_p_compositeoper_files}
      PARENT_SCOPE
)
