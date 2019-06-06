/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSOPERATOR_H_
#define UFO_OBSOPERATOR_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

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
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
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

/// Operator locations
  Locations * locations(const util::DateTime &, const util::DateTime &) const;

/// Operator type
  const std::string & obstype() const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsOperatorBase> oper_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSOPERATOR_H_
