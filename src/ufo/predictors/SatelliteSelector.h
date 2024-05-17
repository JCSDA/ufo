/*
 * (C) Crown copyright 2021, MetOffice
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_SATELLITESELECTOR_H_
#define UFO_PREDICTORS_SATELLITESELECTOR_H_

#include <memory>
#include <string>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/predictors/PredictorBase.h"

namespace eckit {
  class LocalConfiguration;
}

namespace oops {
  class ObsVariables;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

/// Configuration parameters of the SatelliteSelector wrapper for a predictor.
class SatelliteSelectorParameters : public PredictorParametersBase {
  OOPS_CONCRETE_PARAMETERS(SatelliteSelectorParameters, PredictorParametersBase);

 public:
  /// A wrapper to allow for satellite selection when applying any of the
  /// predictors.  Where the satid is not in the list the predictor is set
  /// to zero.
  ///
  /// The satellite id of the satellite which the predictor should apply to.
  oops::RequiredParameter<int> satelliteId{"satellite id", this};

  /// The configuration for a specific predictor
  oops::RequiredParameter<eckit::LocalConfiguration> predictor{"owned predictor", this};

  /// Name for the metadata item which will be compared to allow for:
  /// MetaData/satelliteIdentifier - the default
  /// MetaData/satelliteIdentifier
  /// This will hopefully be removed in the future.
  oops::Parameter<std::string> metadataName{"metadata name", "satelliteIdentifier", this};
};

// -----------------------------------------------------------------------------

class SatelliteSelector : public PredictorBase {
 public:
  /// The type of parameters accepted by the constructor of this predictor.
  /// This typedef is used by the PredictorFactory.
  typedef SatelliteSelectorParameters Parameters_;

  SatelliteSelector(const Parameters_ &, const oops::ObsVariables &);

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               const ObsBias &,
               ioda::ObsVector &) const override;

 private:
  /// The local predictor specified from yaml
  std::unique_ptr<PredictorBase> predictor_;
  const int satid_;
  const std::string metadata_name_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_SATELLITESELECTOR_H_
