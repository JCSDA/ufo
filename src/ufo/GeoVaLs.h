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

#include "eckit/mpi/Comm.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/GeoVaLs.interface.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class Locations;

// -----------------------------------------------------------------------------

/// GeoVaLs: geophysical values at locations

class GeoVaLs : public util::Printable,
                private util::ObjectCounter<GeoVaLs> {
 public:
  static const std::string classname() {return "ufo::GeoVaLs";}

  explicit GeoVaLs(const eckit::mpi::Comm &);
  GeoVaLs(const Locations &, const oops::Variables &);
  GeoVaLs(const eckit::Configuration &, const ioda::ObsSpace &,
          const oops::Variables &);
  GeoVaLs(const GeoVaLs &, const int &);
  GeoVaLs(const GeoVaLs &);

  ~GeoVaLs();

  GeoVaLs & operator = (const GeoVaLs &);
  GeoVaLs & operator *= (const double);
  GeoVaLs & operator *= (const std::vector<float> &);
  GeoVaLs & operator += (const GeoVaLs &);
  GeoVaLs & operator -= (const GeoVaLs &);
  GeoVaLs & operator *= (const GeoVaLs &);
  double dot_product_with(const GeoVaLs &) const;
  void split(GeoVaLs &, GeoVaLs &) const;
  void merge(const GeoVaLs &, const GeoVaLs &);

  void zero();
  void reorderzdir(const std::string &, const std::string &);
  void random();
  double rms() const;
  double normalizedrms(const GeoVaLs &) const;

  bool has(const std::string & var) const {return vars_.has(var);}
  const oops::Variables & getVars() const {return vars_;}

  size_t nlevs(const std::string & var) const;
  void get(std::vector<float> &, const std::string &) const;
  void get(std::vector<float> &, const std::string &, const int) const;
  void get(std::vector<double> &, const std::string &, const int) const;
  void put(const std::vector<double> &, const std::string &, const int) const;

  void read(const eckit::Configuration &, const ioda::ObsSpace &);
  void write(const eckit::Configuration &) const;
  size_t nlocs() const;

  int & toFortran() {return keyGVL_;}
  const int & toFortran() const {return keyGVL_;}

 private:
  void print(std::ostream &) const;

  F90goms keyGVL_;
  oops::Variables vars_;
  const eckit::mpi::Comm & comm_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_GEOVALS_H_
