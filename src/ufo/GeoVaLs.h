/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_GEOVALS_H_
#define UFO_GEOVALS_H_

#include <ostream>
#include <string>

#include "Fortran.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace ufo {
  class ObsSpace;

/// GeoVaLs: geophysical values at locations

class GeoVaLs : public util::Printable,
                private util::ObjectCounter<GeoVaLs> {
 public:
  static const std::string classname() {return "ufo::GeoVaLs";}

  GeoVaLs(const ObsSpace &, const oops::Variables &,
          const util::DateTime &, const util::DateTime &);
  GeoVaLs(const eckit::Configuration &);

  explicit GeoVaLs(): keyGVL_(0) {}
  explicit GeoVaLs(int & fgvl): keyGVL_(fgvl) {}

  ~GeoVaLs();

  void zero();
  void random();
  double dot_product_with(const GeoVaLs & other) const;
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

  int & toFortran() {return keyGVL_;}
  const int & toFortran() const {return keyGVL_;}

 private:
  void print(std::ostream &) const;
  F90goms keyGVL_;
};

}  // namespace ufo

#endif  // UFO_GEOVALS_H_
