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

#include <boost/noncopyable.hpp>

#include "oops/util/Printable.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

// -----------------------------------------------------------------------------

class ObsDiagnostics : public util::Printable,
                       private boost::noncopyable {
 public:
  ObsDiagnostics(const ioda::ObsSpace &, const oops::Variables &);
  ~ObsDiagnostics() {}

// I/O
  void save(const std::string &) const;

 private:
  void print(std::ostream &) const;
  const ioda::ObsSpace & obsdb_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSDIAGNOSTICS_H_
