/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_SINEOFLATITUDE_H_
#define UFO_PREDICTORS_SINEOFLATITUDE_H_

#include <string>
#include <vector>

#include "ufo/predictors/PredictorBase.h"

namespace eckit {
  class Configuration;
  class Comm;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

// -----------------------------------------------------------------------------

class SineOfLatitude : public PredictorBase {
 public:
  SineOfLatitude(const eckit::Configuration &,
                 const std::vector<int> &,
                 const std::string &,
                 const eckit::mpi::Comm &);
  ~SineOfLatitude() {}

  void write(const eckit::Configuration &,
             ObsBiasIO< Record > &) override {}

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               ioda::ObsVector &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_SINEOFLATITUDE_H_
