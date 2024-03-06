/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_MARINE_COOLSKIN_OBSCOOLSKIN_H_
#define UFO_OPERATORS_MARINE_COOLSKIN_OBSCOOLSKIN_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/marine/coolskin/ObsCoolSkin.interface.h"
#include "ufo/operators/marine/coolskin/ObsCoolSkinParameters.h"

/// Forward declarations

namespace ioda {
class ObsSpace;
class ObsVector;
}  // namespace ioda

namespace ufo {
class GeoVaLs;
class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// coolskin observation operator class
class ObsCoolSkin : public ObsOperatorBase,
                    private util::ObjectCounter<ObsCoolSkin> {
 public:
  static const std::string classname() { return "ufo::ObsCoolSkin"; }
  typedef ObsCoolSkinParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  ObsCoolSkin(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsCoolSkin();

  // Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

  // Other
  const oops::Variables &requiredVars() const override { return *varin_; }

  int &toFortran() { return keyOper_; }
  const int &toFortran() const { return keyOper_; }

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  const ioda::ObsSpace & odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_MARINE_COOLSKIN_OBSCOOLSKIN_H_
