/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSLOCALIZATION_OBSVERTLOCPARAMETERS_H_
#define UFO_OBSLOCALIZATION_OBSVERTLOCPARAMETERS_H_

#include <cmath>
#include <string>
#include <utility>

#include "oops/base/ObsLocalizationParametersBase.h"
#include "oops/util/missingValues.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"

namespace ufo {

/// \brief Options controlling vertical localization
class ObsVertLocParameters : public oops::ObsLocalizationParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsVertLocParameters, oops::ObsLocalizationParametersBase)

 public:
  /// Localization lengthscale (find all obs within the distance from reference point)
  oops::Parameter<double> lengthscale{"vertical lengthscale",
                        "lengthscale for localization beyond which localization is set to zero",
                        0.0, this};

  oops::Parameter<bool> logTransform{"apply log transformation",
                        "apply localization function to the logarithm of the distance",
                        false, this};

  oops::Parameter<std::string> iodaVerticalCoordinateGroup{"ioda vertical coordinate group",
                       "group in the ioda file that stores vertical coordinate (default MetaData)",
                       "MetaData", this};

  oops::Parameter<std::string> iodaVerticalCoordinate{"ioda vertical coordinate",
                       "field in the ioda file that stores vertical coordinate", this};

  oops::Parameter<std::string> localizationFunction{"localization function",
                       "localization function (e.g. Box Car, Gaspari Cohn, SOAR", this};

  oops::OptionalParameter<int> maxnobs{"max nobs",
                       "maximum number of obs to include in localization", this};

  oops::Parameter<double> SOARexpDecayH{"soar decay",
                       "soar decay",
                       util::missingValue<double>(), this};

  oops::Parameter<bool> assignConstantVcoordToObs{"assign constant vertical coordinate to obs",
                        "assign constant vertical coordinate to obs (e.g. for surface obs)",
                        false, this};

  oops::Parameter<float> constantVcoordValue{"constant vertical coordinate value",
                       "value of the constant vertical coordinate",
                       util::missingValue<float>(), this};

  /// returns distance between points \p p1 and \p p2
  double distance(const double & vCoord1, const double & vCoord2) const {
      return std::abs(vCoord1 - vCoord2);
    }
};

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSVERTLOCPARAMETERS_H_
