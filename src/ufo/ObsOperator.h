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

#include "ufo/ObsOperatorBase.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
  class ObsBias;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------

class ObsOperator : public util::Printable,
                    private boost::noncopyable {
 public:
  typedef ObsOperatorParametersWrapper Parameters_;

  ObsOperator(ioda::ObsSpace &, const Parameters_ &);

/// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, const ObsBias &, ioda::ObsVector &,
                   ObsDiagnostics &) const;

/// Operator input required from Model
  const oops::Variables & requiredVars() const;

/// Operator locations
  std::unique_ptr<Locations> locations() const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsOperatorBase> oper_;
  ioda::ObsSpace & odb_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSOPERATOR_H_
