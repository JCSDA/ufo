# (C) Copyright 2017-2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( columnretrieval_files
    ObsColumnRetrieval.h
    ObsColumnRetrieval.cc
    ObsColumnRetrieval.interface.F90
    ObsColumnRetrieval.interface.h
    ufo_columnretrieval_mod.F90
    ufo_satcolumn_mod.F90
    ObsColumnRetrievalParameters.h
    ObsColumnRetrievalTLAD.h
    ObsColumnRetrievalTLAD.cc
    ObsColumnRetrievalTLAD.interface.F90
    ObsColumnRetrievalTLAD.interface.h
    ufo_columnretrieval_tlad_mod.F90
)

PREPEND( _p_columnretrieval_files     "operators/columnretrieval"     ${columnretrieval_files} )

set ( columnretrieval_src_files
      ${_p_columnretrieval_files}
      PARENT_SCOPE
)
