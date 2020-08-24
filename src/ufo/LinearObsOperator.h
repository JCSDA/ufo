/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_LINEAROBSOPERATOR_H_
#define UFO_LINEAROBSOPERATOR_H_

#include <Eigen/Core>

#include <memory>
#include <vector>

#include <boost/noncopyable.hpp>

#include "ioda/ObsDataVector.h"

#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;
  class LinearObsOperatorBase;

// -----------------------------------------------------------------------------

class LinearObsOperator : public util::Printable,
                          private boost::noncopyable {
 public:
  LinearObsOperator(ioda::ObsSpace &, const eckit::Configuration &);
  ~LinearObsOperator();

/// Obs Operator
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

/// Operator input required from Model
  const oops::Variables & requiredVars() const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<LinearObsOperatorBase> oper_;
  ioda::ObsSpace & odb_;
  std::unique_ptr<ioda::ObsDataVector<double>> biaspreds_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_LINEAROBSOPERATOR_H_
