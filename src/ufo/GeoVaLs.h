/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_GEOVALS_H_
#define UFO_GEOVALS_H_

#include <memory>
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

namespace ioda {
  class ObsSpace;
  class Distribution;
}

namespace ufo {
  class Locations;

// -----------------------------------------------------------------------------

/// GeoVaLs: geophysical values at locations

class GeoVaLs : public util::Printable,
                private util::ObjectCounter<GeoVaLs> {
 public:
  static const std::string classname() {return "ufo::GeoVaLs";}

  GeoVaLs(std::shared_ptr<const ioda::Distribution>, const oops::Variables &);
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

  /// \brief Allocate GeoVaLs for \p vars variables with \p nlev number of levels
  /// \details Fails if at least one of the \p vars doesn't exist in GeoVaLs.
  ///          Only allocates variables that haven't been allocated before.
  ///          Fails if one of \p vars is already allocated with number of levels
  ///          different than \p nlev; doesn't reallocate variables that are already
  ///          allocated with \p nlev.
  void allocate(const int & nlev, const oops::Variables & vars);

  void zero();
  void reorderzdir(const std::string &, const std::string &);
  void random();
  double rms() const;
  double normalizedrms(const GeoVaLs &) const;

  bool has(const std::string & var) const {return vars_.has(var);}
  const oops::Variables & getVars() const {return vars_;}

  size_t nlevs(const std::string & var) const;
  void get(std::vector<float> &, const std::string &, const int) const;
  void get(std::vector<double> &, const std::string &, const int) const;
  /// Get 2D GeoVaLs for variable \p var (fails for 3D GeoVaLs)
  void get(std::vector<double> &, const std::string & var) const;
  /// Get 2D GeoVaLs for variable \p var (fails for 3D GeoVaLs), and convert to float
  void get(std::vector<float> &, const std::string & var) const;
  /// Get 2D GeoVaLs for variable \p var (fails for 3D GeoVaLs), and convert to int
  void get(std::vector<int> &, const std::string & var) const;
  /// Get GeoVaLs at a specified location
  void getAtLocation(std::vector<double> &, const std::string &, const int) const;
  /// Get GeoVaLs at a specified location and convert to float
  void getAtLocation(std::vector<float> &, const std::string &, const int) const;
  /// Get GeoVaLs at a specified location and convert to int
  void getAtLocation(std::vector<int> &, const std::string &, const int) const;
  /// Put GeoVaLs for double variable \p var at level \p lev.
  void put(const std::vector<double> & vals, const std::string & var, const int lev) const;
  /// Put GeoVaLs for float variable \p var at level \p lev.
  void put(const std::vector<float> & vals, const std::string & var, const int lev) const;
  /// Put GeoVaLs for int variable \p var at level \p lev.
  void put(const std::vector<int> & vals, const std::string & var, const int lev) const;
  /// Put GeoVaLs for double variable \p var at location \p loc.
  void putAtLocation(const std::vector<double> & vals, const std::string & var,
                     const int loc) const;
  /// Put GeoVaLs for float variable \p var at location \p loc.
  void putAtLocation(const std::vector<float> & vals, const std::string & var, const int loc) const;
  /// Put GeoVaLs for int variable \p var at location \p loc.
  void putAtLocation(const std::vector<int> & vals, const std::string & var, const int loc) const;
  void read(const eckit::Configuration &, const ioda::ObsSpace &);
  void write(const eckit::Configuration &) const;
  size_t nlocs() const;

  int & toFortran() {return keyGVL_;}
  const int & toFortran() const {return keyGVL_;}

 private:
  void print(std::ostream &) const;

  F90goms keyGVL_;
  oops::Variables vars_;
  std::shared_ptr<const ioda::Distribution> dist_;   /// observations MPI distribution
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_GEOVALS_H_
