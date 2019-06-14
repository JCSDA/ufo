/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_RTTOV_OBSRADIANCERTTOV_H_
#define UFO_RTTOV_OBSRADIANCERTTOV_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/rttov/ObsRadianceRTTOV.interface.h"

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

// -----------------------------------------------------------------------------
/// RadianceRTTOV observation operator class
class ObsRadianceRTTOV : public ObsOperatorBase,
                   private util::ObjectCounter<ObsRadianceRTTOV> {
 public:
  static const std::string classname() {return "ufo::ObsRadianceRTTOV";}

  ObsRadianceRTTOV(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsRadianceRTTOV();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &) const override;

// Other
  const oops::Variables & variables() const override {return *varin_;}
  const std::string & obstype() const override {return obsname_;}

  int & toFortran() {return keyOperRadianceRTTOV_;}
  const int & toFortran() const {return keyOperRadianceRTTOV_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperRadianceRTTOV_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
  std::vector<int> channels_;
  std::string obsname_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_RTTOV_OBSRADIANCERTTOV_H_
