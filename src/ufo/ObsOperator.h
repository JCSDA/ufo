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

#include "ioda/ObsDataVector.h"
#include "oops/util/Printable.h"

#include "ufo/ObsOperatorBase.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  template <typename OBS> class Locations;
  class Variables;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsDiagnostics;
  struct ObsTraits;

// -----------------------------------------------------------------------------

class ObsOperator : public util::Printable,
                    private boost::noncopyable {
 public:
  typedef oops::Locations<ObsTraits> Locations_;
  typedef ObsOperatorParametersWrapper Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;
  ObsOperator(ioda::ObsSpace &, const eckit::Configuration &);

/// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, const ObsBias &,
                   const QCFlags_t &, ioda::ObsVector &, ObsDiagnostics &) const;
/// Operator input required from Model
  const oops::Variables & requiredVars() const;

/// Model variable interpolation paths
  Locations_ locations() const;

  void computeReducedVars(const oops::Variables & vars, GeoVaLs & geovals) const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsOperatorBase> oper_;
  ioda::ObsSpace & odb_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSOPERATOR_H_
