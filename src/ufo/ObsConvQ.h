/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSCONVQ_H_
#define UFO_OBSCONVQ_H_

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

class ObsConvQ : public oops::ObsOperatorBase<UfoTrait>,
                  private util::ObjectCounter<ObsConvQ> {
 public:
  static const std::string classname() {return "ufo::ObsConvQ";}

  ObsConvQ(const ObsSpace &, const eckit::Configuration &);
  virtual ~ObsConvQ();

// Obs Operator
  void obsEquiv(const GeoVaLs &, ObsVector &, const ObsBias &) const;

// Other
  boost::shared_ptr<const Variables> variables() const {return varin_;}

  int & toFortran() {return keyOperConvQ_;}
  const int & toFortran() const {return keyOperConvQ_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperConvQ_;
  boost::shared_ptr<const Variables> varin_;
};
// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSCONVQ_H_
