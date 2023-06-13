/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVALPARAMETERS_H_
#define UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVALPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// Configuration options recognized by the ColumnRetrieval operator.
class ObsColumnRetrievalParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsColumnRetrievalParameters, ObsOperatorParametersBase)

 public:
  oops::Parameter <int> nlayers_retrieval
    {"nlayers_retrieval",
     "Number of layers in the retrieval kernel",
     1,
     this};

  oops::RequiredParameter <std::vector<std::string>> tracerVariables
    {"tracer variables",
     "Names of model tracer variables",
     this};

  oops::Parameter <bool> isApriori
    {"isApriori",
     "Add the a priori retrieval term to the AK smoothed model profile",
     false,
     this};

  oops::Parameter <bool> isAveragingKernel
    {"isAveragingKernel",
     "Add the AK to smooth the model profile",
     false,
     this};

  oops::Parameter <std::string> stretchVertices
    {"stretchVertices",
     "match Obs and Bkg vertices, options: top, bottom, topbottom and none (default)",
     "none",
     this};

  oops::Parameter <double> modelUnitsCoeff
    {"model units coeff",
     "Conversion between model units",
     1.0,
     this};

  oops::Parameter <bool> totalNoVertice
    {"totalNoVertice",
     "Total column calculation without AK and apriori (no vertices needed)",
     false,
     this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVALPARAMETERS_H_
