/*
 * (C) Copyright 2022 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_OASIM_OBSRADIANCEOASIM_H_
#define UFO_OPERATORS_OASIM_OBSRADIANCEOASIM_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/oasim/ObsRadianceOASIM.interface.h"
#include "ufo/operators/oasim/ObsRadianceOASIMParameters.h"

/// Forward declarations

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// OASIM observation operator class
class ObsRadianceOASIM : public ObsOperatorBase,
                   private util::ObjectCounter<ObsRadianceOASIM> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsRadianceOASIM";}

  typedef ObsRadianceOASIMParameters Parameters_;

  ObsRadianceOASIM(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsRadianceOASIM();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;

  ObsRadianceOASIMParameters parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_OASIM_OBSRADIANCEOASIM_H_
