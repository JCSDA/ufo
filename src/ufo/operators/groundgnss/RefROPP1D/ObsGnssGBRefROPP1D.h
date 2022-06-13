/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_GROUNDGNSS_REFROPP1D_OBSGNSSGBREFROPP1D_H_
#define UFO_OPERATORS_GROUNDGNSS_REFROPP1D_OBSGNSSGBREFROPP1D_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/groundgnss/RefROPP1D/ObsGnssGBRefROPP1D.interface.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------

/// GnssGBRefROPP1D observation operator
class ObsGnssGBRefROPP1D : public ObsOperatorBase,
                        private util::ObjectCounter<ObsGnssGBRefROPP1D> {
 public:
  static const std::string classname() {return "ufo::ObsGnssGBRefROPP1D";}

  ObsGnssGBRefROPP1D(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsGnssGBRefROPP1D();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssGBRefROPP1D_;}
  const int & toFortran() const {return keyOperGnssGBRefROPP1D_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssGBRefROPP1D_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_GROUNDGNSS_REFROPP1D_OBSGNSSGBREFROPP1D_H_
