/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORLATRAD_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORLATRAD_H_

#include <vector>

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Options controlling the observation error bound reduction in Tropical regions
///
class ObsErrorFactorLatRadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorFactorLatRadParameters, Parameters)

 public:
  /// Parameters for reducing observation error bounds within latitude band in Tropics
  /// params[0] defines the latitude bound for which the observation error function applies.
  /// params[1-3] determine the error function within the latitude bound given by params[0].
  /// The error function gives the maximum error bound reduction at equator and decreasing
  /// towards params[0].
  /// Error Function = params[1] * ( |latitude| * params[2] + params[3] )
  oops::RequiredParameter<std::vector<float>> latitudeParameters{"latitude_parameters", this};
};

///
/// \brief Function determines the observation error bound reduction within Tropics.
/// The function gives the maximum error bound reduction at equator and decreasing
/// towards higher latitudes.
///
class ObsErrorFactorLatRad : public ObsFunctionBase<float> {
 public:
  explicit ObsErrorFactorLatRad(const eckit::LocalConfiguration &);
  ~ObsErrorFactorLatRad();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ObsErrorFactorLatRadParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORLATRAD_H_
