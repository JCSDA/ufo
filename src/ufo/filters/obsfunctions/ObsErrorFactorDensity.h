/*
 * (C) Copyright 2026 NOAA/NWS/NCEP/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORDENSITY_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORDENSITY_H_

#include <string>

#include "oops/util/Duration.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling ObsErrorFactorDensity ObsFunction
class ObsErrorFactorDensityParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorFactorDensityParameters, Parameters)

 public:
  /// Variable whose QC flags decide which observations are counted
  oops::RequiredParameter<std::string> variable{"variable", this};
  /// Name of the data group holding the QC flags (only observations with flag 0 are counted)
  oops::Parameter<std::string> testQCflag{"test QCflag", "QCflagsData", this};
  /// Horizontal size (km) of the counting volume
  oops::RequiredParameter<float> horizontalMesh{"horizontal mesh", this,
                                                {oops::exclusiveMinConstraint(0.0f)}};
  /// Vertical coordinate, e.g. MetaData/height
  oops::RequiredParameter<Variable> verticalCoordinate{"vertical coordinate", this};
  /// Vertical size of the counting volume, in the units of the vertical coordinate
  oops::RequiredParameter<float> verticalMesh{"vertical mesh", this,
                                              {oops::exclusiveMinConstraint(0.0f)}};
  /// Length of the time intervals, starting at the beginning of the window (must be positive).
  /// Like the window, each interval excludes its beginning and includes its end.
  oops::RequiredParameter<util::Duration> timeMesh{"time mesh", this};
  /// Density (observations per box of the mesh sizes) above which errors are inflated
  oops::RequiredParameter<float> densityThreshold{"density threshold", this,
                                                  {oops::exclusiveMinConstraint(0.0f)}};
};

// -----------------------------------------------------------------------------

/// \brief Error inflation factor for observations in dense regions.
///
/// The local density is the number of observations passing QC within a cylinder centred on each
/// observation, in the same time interval and across all MPI tasks. The cylinder has the horizontal
/// area of a square of the horizontal mesh size (radius mesh/sqrt(pi)) and the height of the
/// vertical mesh size, so its volume is that of a box of the mesh sizes. Each observation gets the
/// factor sqrt(max(1, density / threshold)), so that dense observations carry about the weight of
/// `threshold` observations per box. This follows the GSI down-weighting of dense aircraft winds
/// (read_prepbufr), which counts observations in the same volume along randomized Hilbert curves
/// and inflates the error by sqrt(density / threshold) where the density exceeds the threshold.
///
/// Observations with a missing location, time or vertical coordinate get a factor of 1 and are not
/// counted.
///
/// ### example configuration ###
///
///     - filter: Perform Action
///       filter variables:
///       - name: windEastward
///       - name: windNorthward
///       action:
///         name: inflate error
///         inflation variable:
///           name: ObsFunction/ObsErrorFactorDensity
///           options:
///             variable: windEastward
///             horizontal mesh: 100.
///             vertical coordinate: MetaData/height
///             vertical mesh: 800.
///             time mesh: PT1H
///             density threshold: 2.
///
class ObsErrorFactorDensity : public ObsFunctionBase<float> {
 public:
  explicit ObsErrorFactorDensity(const eckit::LocalConfiguration &);
  ~ObsErrorFactorDensity();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;

 private:
  ufo::Variables invars_;
  ObsErrorFactorDensityParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORDENSITY_H_
