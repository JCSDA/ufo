/*
 * (C) Copyright 2018  UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSDIAGNOSTICS_H_
#define UFO_OBSDIAGNOSTICS_H_

#include <ostream>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/util/Printable.h"
#include "ufo/GeoVaLs.h"

// Forward declarations
namespace oops {
  template <typename OBS> class Locations;
  class Variables;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  struct ObsTraits;

// -----------------------------------------------------------------------------

class ObsDiagnostics : public util::Printable,
                       private boost::noncopyable {
 public:
  typedef GeoVaLs::Parameters_ Parameters_;
  typedef oops::Locations<ObsTraits> Locations_;

  ObsDiagnostics(const ioda::ObsSpace &, const Locations_ &,
                 const oops::Variables &);
  ObsDiagnostics(const Parameters_ &, const ioda::ObsSpace &,
                 const oops::Variables &);
  ~ObsDiagnostics() {}

  /// \brief Allocate diagnostics for variables \p vars with \p nlev number of levels
  /// \details Fails if at least one of the \p vars doesn't exist in the ObsDiagnostics object.
  ///          Only allocates variables that haven't been allocated before.
  ///          Fails if one of \p vars is already allocated with a number of levels
  ///          different than \p nlev; doesn't reallocate variables that are already
  ///          allocated with \p nlev.
  void allocate(const int nlev, const oops::Variables & vars);

  void save(const std::vector<double> &, const std::string &, const int);

// Interfaces
  int & toFortran() {return gdiags_.toFortran();}
  const int & toFortran() const {return gdiags_.toFortran();}

  bool has(const std::string & var) const {return gdiags_.has(var);}
  size_t nlevs(const std::string &) const;
  template <typename T>
  void get(std::vector<T> & vals, const std::string & var, const int lev) const {
    gdiags_.getAtLevel(vals, var, lev);
  }
  template <typename T>
  void get(std::vector<T> & vals, const std::string & var) const {
    gdiags_.get(vals, var);
  }

  void write(const Parameters_ & params) const {
    gdiags_.write(params);}
 private:
  void print(std::ostream &) const;
  const ioda::ObsSpace & obsdb_;

  GeoVaLs gdiags_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSDIAGNOSTICS_H_
