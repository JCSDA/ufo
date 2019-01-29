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
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/GeoVaLs.interface.h"

namespace eckit {
  class Configuration;
}

namespace ufo {
  class Locations;

// -----------------------------------------------------------------------------

/// GeoVaLs: geophysical values at locations

class GeoVaLs : public util::Printable,
                private util::ObjectCounter<GeoVaLs> {
 public:
  static const std::string classname() {return "ufo::GeoVaLs";}

  GeoVaLs(const Locations &, const oops::Variables &);
  GeoVaLs(const eckit::Configuration &, const oops::Variables &);
  GeoVaLs(const GeoVaLs &);

  ~GeoVaLs();

  GeoVaLs & operator = (const GeoVaLs &);
  GeoVaLs & operator*=(const double);
  GeoVaLs & operator += (const GeoVaLs &);
  GeoVaLs & operator -= (const GeoVaLs &);
  GeoVaLs & operator /= (const GeoVaLs &);
  double dot_product_with(const GeoVaLs &) const;

  void abs();
  void zero();
  void random();
  double norm() const;

  bool has(const std::string & var) const {return vars_.has(var);}
  void get(std::vector<float> &, const std::string &, const int lev = 1) const;

  void read(const eckit::Configuration &);
  void analytic_init(const Locations &, const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

  int & toFortran() {return keyGVL_;}
  const int & toFortran() const {return keyGVL_;}

 private:
  void print(std::ostream &) const;

  F90goms keyGVL_;
  oops::Variables vars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_GEOVALS_H_
