/*
 *
 * (C) Copyright 2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_SATTCWV_SATTCWV_H_
#define UFO_SATTCWV_SATTCWV_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

  // -----------------------------------------------------------------------------
  /*! \brief Precipitable water observation operator
  *
  * \details This is based on Roger Saunders' Fortran original "sattcwv" 
  * observation operator and also derives from Heather Lawrence's MetOffice AOD
  * forward operator code. 
  *
  * The TCWV (total column water vapour) is calculated as:
  * \f[ TCWV = \sum_{k=1}^{Nlevels}(1/g) q(k) |\Delta P(k)| \f],
  *
  * where: \f$g\f$ is the gravity constant and \f$q(k)\f$ is the specific humidity
  * per atmospheric level. All inputs are taken from the geovals. Note that the 
  * geoval input of pressure needs to be on staggered levels relative to specific
  * humidity so that \f$\|\Delta P|\f$ is at the same height as specific humidity.
  *
  * \date Sept. 2021: Created by J. Hocking (Met Office)
  */
  // -----------------------------------------------------------------------------
class SatTCWV : public ObsOperatorBase,
                  private util::ObjectCounter<SatTCWV> {
 public:
  static const std::string classname() {return "ufo::SatTCWV";}

  SatTCWV(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~SatTCWV();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<const oops::Variables> varin_;  // list of the required geovals
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_SATTCWV_SATTCWV_H_
