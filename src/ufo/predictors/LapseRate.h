/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_LAPSERATE_H_
#define UFO_PREDICTORS_LAPSERATE_H_

#include <map>

#include "ufo/predictors/PredictorBase.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

// -----------------------------------------------------------------------------

class LapseRate : public PredictorBase {
 public:
  LapseRate(const eckit::Configuration &, const oops::Variables &);

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               ioda::ObsVector &) const override;

 private:
  std::map<int, float> tlapmean_;  // <channel, tlaps>
  int order_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_LAPSERATE_H_
