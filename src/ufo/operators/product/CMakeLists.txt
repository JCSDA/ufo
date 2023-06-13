# (C) Copyright 2023- UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( product_files
    ObsProduct.h
    ObsProduct.cc
    ObsProductParameters.h
    ObsProductTLAD.h
    ObsProductTLAD.cc
)

PREPEND( _p_product_files       "operators/product"       ${product_files} )

set ( product_src_files
      ${_p_product_files}
      PARENT_SCOPE
)

