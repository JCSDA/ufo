/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_RTTOV_OBSRADIANCERTTOV_H_
#define UFO_RTTOV_OBSRADIANCERTTOV_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

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
  class ObsBias;

// -----------------------------------------------------------------------------
/// RadianceRTTOV observation operator class
class ObsRadianceRTTOV : public ObsOperatorBase,
                   private util::ObjectCounter<ObsRadianceRTTOV> {
 public:
  static const std::string classname() {return "ufo::ObsRadianceRTTOV";}

  ObsRadianceRTTOV(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsRadianceRTTOV();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}
  const oops::Variables & observed() const {return *varout_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
  boost::scoped_ptr<const oops::Variables> varout_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_RTTOV_OBSRADIANCERTTOV_H_
