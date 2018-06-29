/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_GEOVALS_H_
#define UFO_GEOVALS_H_

#include <ostream>
#include <string>

#include "Fortran.h"
#include "FortranGeoVals.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class Locations;
}

namespace ufo {

/// GeoVaLs: geophysical values at locations

class GeoVaLs : public util::Printable,
                private util::ObjectCounter<GeoVaLs> {
 public:
  static const std::string classname() {return "ufo::GeoVaLs";}

  GeoVaLs(const ioda::Locations &, const oops::Variables &);
  GeoVaLs(const eckit::Configuration &, const oops::Variables &);
  GeoVaLs(const GeoVaLs &);

//  explicit GeoVaLs(): keyGVL_(0) {}
//  explicit GeoVaLs(int & fgvl): keyGVL_(fgvl) {}

  ~GeoVaLs();

  void abs();
  void zero();
  void random();
  double norm() const;
  GeoVaLs & operator = (const GeoVaLs &);    
  GeoVaLs & operator*=(const double);
  GeoVaLs & operator += (const GeoVaLs &);
  GeoVaLs & operator -= (const GeoVaLs &);
  GeoVaLs & operator /= (const GeoVaLs &);
  double dot_product_with(const GeoVaLs & other) const;
  void read(const eckit::Configuration &);
  void analytic_init(const ioda::Locations &, const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

  int & toFortran() {return keyGVL_;}
  const int & toFortran() const {return keyGVL_;}

 private:
  void print(std::ostream &) const;

  F90goms keyGVL_;
  oops::Variables vars_;
};

}  // namespace ufo

#endif  // UFO_GEOVALS_H_
