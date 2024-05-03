/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_GNSSRO_BNDROPP1D_OBSGNSSROBNDROPP1DTLAD_H_
#define UFO_OPERATORS_GNSSRO_BNDROPP1D_OBSGNSSROBNDROPP1DTLAD_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/gnssro/BndROPP1D/ObsGnssroBndROPP1D.h"
#include "ufo/operators/gnssro/BndROPP1D/ObsGnssroBndROPP1DTLAD.interface.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// GnssroBndROPP1D observation operator
class ObsGnssroBndROPP1DTLAD : public LinearObsOperatorBase,
                          private util::ObjectCounter<ObsGnssroBndROPP1DTLAD> {
 public:
  typedef GnssroBndROPP1DParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsGnssroBndROPP1DTLAD";}

  ObsGnssroBndROPP1DTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsGnssroBndROPP1DTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssroBndROPP1D_;}
  const int & toFortran() const {return keyOperGnssroBndROPP1D_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssroBndROPP1D_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_GNSSRO_BNDROPP1D_OBSGNSSROBNDROPP1DTLAD_H_
