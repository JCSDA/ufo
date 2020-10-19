/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_AVGKERNEL_OBSAVGKERNEL_H_
#define UFO_AVGKERNEL_OBSAVGKERNEL_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/avgkernel/ObsAvgKernel.interface.h"
#include "ufo/ObsOperatorBase.h"

/// Forward declarations
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
/// AvgKernel observation operator class
class ObsAvgKernel : public ObsOperatorBase,
                   private util::ObjectCounter<ObsAvgKernel> {
 public:
  static const std::string classname() {return "ufo::ObsAvgKernel";}

  ObsAvgKernel(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAvgKernel();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_AVGKERNEL_OBSAVGKERNEL_H_
