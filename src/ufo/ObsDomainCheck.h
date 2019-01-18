/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSDOMAINCHECK_H_
#define UFO_OBSDOMAINCHECK_H_

#include <ostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace ioda {
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

/// Domain check: generic check that obs are within domain

// Domain is defined by metadata criteria regardless of obs value.
// If obs value is required, use ObsBoundsCheck.

// The same effect can be achieved with opposite criteria through BlackList,
// the choice is a matter of convenience or which seems more natural.

class ObsDomainCheck : public util::Printable,
                       private util::ObjectCounter<ObsDomainCheck> {
 public:
  static const std::string classname() {return "ufo::ObsDomainCheck";}

  ObsDomainCheck(ioda::ObsSpace &, const eckit::Configuration &);
  ~ObsDomainCheck();

  void priorFilter(const GeoVaLs &) const;
  void postFilter(const ioda::ObsVector &) const {}

  const oops::Variables & requiredGeoVaLs() const {return geovars_;}

 private:
  void print(std::ostream &) const;

  ioda::ObsSpace & obsdb_;
  const eckit::LocalConfiguration config_;
  const oops::Variables geovars_;
};

}  // namespace ufo

#endif  // UFO_OBSDOMAINCHECK_H_
