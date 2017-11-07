/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSWSPEED_H_
#define UFO_OBSWSPEED_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "oops/interface/ObsOperatorBase.h"
#include "ObsSpace.h"
#include "UfoTrait.h"
#include "util/ObjectCounter.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsVector;

// -----------------------------------------------------------------------------
/// Wind speed observation for UFO.

class ObsWSpeed : public oops::ObsOperatorBase<UfoTrait>,
                  private util::ObjectCounter<ObsWSpeed> {
 public:
  static const std::string classname() {return "ufo::ObsWSpeed";}

  ObsWSpeed(const ObsSpace &, const eckit::Configuration &);
  virtual ~ObsWSpeed();

// Obs Operator
  void obsEquiv(const GeoVaLs &, ObsVector &, const ObsBias &) const;

// Other
  boost::shared_ptr<const Variables> variables() const {return varin_;}

  int & toFortran() {return keyOperWspeed_;}
  const int & toFortran() const {return keyOperWspeed_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperWspeed_;
  boost::shared_ptr<const Variables> varin_;
};
// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSWSPEED_H_
