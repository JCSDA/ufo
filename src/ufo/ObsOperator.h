/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSOPERATOR_H_
#define UFO_OBSOPERATOR_H_

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace util {
  class DateTime;
}

namespace ioda {
  class Locations;
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsOperatorBase;

// -----------------------------------------------------------------------------

class ObsOperator : public util::Printable,
                    private boost::noncopyable {
 public:
  ObsOperator(const ioda::ObsSpace &, const eckit::Configuration &);
  ~ObsOperator();

/// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

/// Operator input required from Model
  const oops::Variables & variables() const;

// TODO(YT): move locations method down to actual obs operators
//           and remove obsdb_ from here
  ioda::Locations * locations(const util::DateTime &, const util::DateTime &) const;

 private:
  void print(std::ostream &) const;
  const ioda::ObsSpace & obsdb_;
  boost::scoped_ptr<ObsOperatorBase> oper_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSOPERATOR_H_
