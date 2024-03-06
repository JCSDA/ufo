/*
 * (C) Copyright 2021- UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLE_H_
#define TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLE_H_

#include <ostream>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"

#include "ufo/example/ObsExampleParameters.h"
#include "ufo/Fortran.h"
#include "ufo/ObsOperatorBase.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// Example observation operator class
class ObsExample : public ObsOperatorBase {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsExampleParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t; 

  ObsExample(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsExample();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

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
#endif  // TOOLS_NEW_OBSOP_EXAMPLE_OBSEXAMPLE_H_
