/*
 * (C) Copyright 2022-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVALPARAMETERS_H_
#define UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVALPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// How the obs and model pressure vertices are aligned before integration.
enum class StretchVertices {
  none, top, bottom, topbottom
};

struct StretchVerticesParameterTraitsHelper {
  typedef StretchVertices EnumType;
  static constexpr char enumTypeName[] = "StretchVertices";
  static constexpr util::NamedEnumerator<StretchVertices> namedValues[] = {
    { StretchVertices::none,      "none" },
    { StretchVertices::top,       "top" },
    { StretchVertices::bottom,    "bottom" },
    { StretchVertices::topbottom, "topbottom" }
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::StretchVertices> :
    public EnumParameterTraits<ufo::StretchVerticesParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// Configuration options recognized by the ColumnRetrieval operator.
class ObsColumnRetrievalParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsColumnRetrievalParameters, ObsOperatorParametersBase)

 public:
  /// Checks that the combination of options is consistent.
  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  oops::Parameter <int> nlayers_retrieval
    {"nlayers_retrieval",
     "Number of layers in the retrieval kernel",
     1,
     this,
     {oops::minConstraint(1)}};

  oops::RequiredParameter <std::vector<std::string>> tracerVariables
    {"tracer variables",
     "Name of the model tracer variable (exactly one; the operator simulates a single variable)",
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

  oops::Parameter <StretchVertices> stretchVertices
    {"stretchVertices",
     "match Obs and Bkg vertices, options: top, bottom, topbottom and none (default)",
     StretchVertices::none,
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

  oops::Parameter <double> tropopause
    {"tropopause pressure",
     "Force tropopause pressure value",
     0.0,
     this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_COLUMNRETRIEVAL_OBSCOLUMNRETRIEVALPARAMETERS_H_
