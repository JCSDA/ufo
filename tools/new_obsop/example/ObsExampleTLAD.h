/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLETLAD_H_
#define TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLETLAD_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/example/ObsExampleParameters.h"
#include "ufo/example/ObsExampleTLAD.interface.h"
#include "ufo/LinearObsOperatorBase.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

// -----------------------------------------------------------------------------
/// Example TL/AD observation operator class
class ObsExampleTLAD : public LinearObsOperatorBase,
                       private util::ObjectCounter<ObsExampleTLAD> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the LinearObsOperatorFactory.
  typedef ObsExampleParameters Parameters_;

  static const std::string classname() {return "ufo::ObsExampleTLAD";}

  ObsExampleTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsExampleTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLETLAD_H_
