# (C) Copyright 2017-2026 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( columnretrieval_files
    ObsColumnRetrieval.h
    ObsColumnRetrieval.cc
    ObsColumnRetrievalParameters.h
    ObsColumnRetrievalParameters.cc
    ObsColumnRetrievalTLAD.h
    ObsColumnRetrievalTLAD.cc
    SatColumnAlgorithm.h
    SatColumnAlgorithm.cc
)

PREPEND( _p_columnretrieval_files     "operators/columnretrieval"     ${columnretrieval_files} )

set ( columnretrieval_src_files
      ${_p_columnretrieval_files}
      PARENT_SCOPE
)
