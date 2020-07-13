/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ANALYTICINIT_H_
#define UFO_ANALYTICINIT_H_

#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/ObjectCounter.h"

namespace ufo {
  class GeoVaLs;
  class Locations;

/// AnalyticInit: filling GeoVaLs with analytic formula
class AnalyticInit : private util::ObjectCounter<AnalyticInit> {
 public:
  static const std::string classname() {return "ufo::AnalyticInit";}

  explicit AnalyticInit(const eckit::Configuration &);
  void fillGeoVaLs(const Locations &, GeoVaLs &) const;

 private:
  const eckit::LocalConfiguration config_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_ANALYTICINIT_H_
