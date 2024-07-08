/*
 * (C) Crown copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_RADARSCANEDGEFLAG_H_
#define UFO_FILTERS_OBSFUNCTIONS_RADARSCANEDGEFLAG_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Options controlling the RadarScanEdgeFlag ObsFunction
class RadarScanEdgeFlagParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(RadarScanEdgeFlagParameters, Parameters)

 public:
  /// Conditions used to select locations at which the RadarScanEdgeFlag
  /// ObsFunction should be applied.
  /// If not specified, all locations will be selected.
  oops::Parameter<std::vector<WhereParameters>> where{"where", {}, this};

  /// Operator used to combine the results of successive `where` options at the same location.
  /// The available operators are `and` and `or`.
  oops::Parameter<WhereOperator> whereOperator{"where operator", WhereOperator::AND, this};

  /// Laplace filter threshold.
  /// If this value is zero or negative, do not apply the filter.
  oops::Parameter<double> thresholdLaplaceFilter{"Laplace filter threshold", 0.0, this};

  /// Cleaning filter threshold.
  /// If this value is zero or negative, do not apply the filter.
  oops::Parameter<int> thresholdCleanFilter{"cleaning filter threshold", 0, this};

  /// Double cleaning filter threshold.
  /// If this value is zero or negative, do not apply the filter.
  oops::Parameter<int> thresholdDoubleCleanFilter{"double cleaning filter threshold", 0, this};

  /// OPS compatibility mode.
  /// If true, the wrapping of the observation value count and count arrays are wrapped
  /// in different ways in the azimuthal direction.
  oops::Parameter<bool> opsCompatibilityMode{"OPS compatibility mode", false, this};
};

// -----------------------------------------------------------------------------

/// \brief Flag noisy edges of radar scans.
///
/// Two methods are used to flag the edges:
/// * Laplace filter,
/// * (Double) edge cleaning.
class RadarScanEdgeFlag : public ObsFunctionBase<int> {
 public:
  explicit RadarScanEdgeFlag(const eckit::LocalConfiguration &);
  ~RadarScanEdgeFlag();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<int> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  RadarScanEdgeFlagParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_RADARSCANEDGEFLAG_H_
