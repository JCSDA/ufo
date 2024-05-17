/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_OBSPROFILEAVERAGETLAD_H_
#define UFO_PROFILE_OBSPROFILEAVERAGETLAD_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"

#include "ufo/profile/ObsProfileAverageData.h"
#include "ufo/profile/ObsProfileAverageParameters.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;
  class ObsBiasIncrements;

/// \brief TL/AD code for the ProfileAverage observation operator.
class ObsProfileAverageTLAD : public LinearObsOperatorBase,
  private util::ObjectCounter<ObsProfileAverageTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsProfileAverageTLAD";}
  typedef ObsProfileAverageParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  ObsProfileAverageTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsProfileAverageTLAD();

  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override { return requiredVars_; }

  oops::ObsVariables simulatedVars() const override { return data_.simulatedVars(); }

 private:
  void print(std::ostream &) const override;

 private:
  /// ObsSpace.
  const ioda::ObsSpace& odb_;

  /// Data handler for the ProfileAverage operator and TL/AD code.
  ObsProfileAverageData data_;

  /// Required variables.
  oops::Variables requiredVars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_PROFILE_OBSPROFILEAVERAGETLAD_H_
