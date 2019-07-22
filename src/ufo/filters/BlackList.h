/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_BLACKLIST_H_
#define UFO_FILTERS_BLACKLIST_H_

#include <ostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "eckit/config/LocalConfiguration.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace ioda {
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

/// BlackList: generic black listing of observations

// Filters observations out regardless of obs value
// The same effect can be achieved with opposite criteria through the "Domain Check",
// the choice is a matter of convenience or which seems more natural.

class BlackList : public util::Printable,
                  private util::ObjectCounter<BlackList> {
 public:
  static const std::string classname() {return "ufo::BlackList";}

  BlackList(ioda::ObsSpace &, const eckit::Configuration &,
            boost::shared_ptr<ioda::ObsDataVector<int> >,
            boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~BlackList();

  void preProcess() const {}
  void priorFilter(const GeoVaLs &) const;
  void postFilter(const ioda::ObsVector &) const {}

  const oops::Variables & requiredGeoVaLs() const {return geovars_;}

 private:
  void print(std::ostream &) const;

  ioda::ObsSpace & obsdb_;
  const eckit::LocalConfiguration config_;
  oops::Variables geovars_;
  ioda::ObsDataVector<int> & flags_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_BLACKLIST_H_
