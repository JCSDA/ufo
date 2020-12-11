/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_AVGKERNEL_OBSAVGKERNELTLAD_H_
#define UFO_AVGKERNEL_OBSAVGKERNELTLAD_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/avgkernel/ObsAvgKernelTLAD.interface.h"
#include "ufo/LinearObsOperatorBase.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;

// -----------------------------------------------------------------------------
/// AvgKernel TL/AD observation operator class
class ObsAvgKernelTLAD : public LinearObsOperatorBase,
                       private util::ObjectCounter<ObsAvgKernelTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsAvgKernelTLAD";}

  ObsAvgKernelTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAvgKernelTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperAvgKernel_;}
  const int & toFortran() const {return keyOperAvgKernel_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperAvgKernel_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_AVGKERNEL_OBSAVGKERNELTLAD_H_
