/*
 *
 * Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_AEROSOLS_AODMETOFFICE_OBSAODMETOFFICE_H_
#define UFO_OPERATORS_AEROSOLS_AODMETOFFICE_OBSAODMETOFFICE_H_

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/operators/aerosols/AODMetOffice/ObsAodMetOfficeParameters.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/ObsOperatorParametersBase.h"


/// Forward declarations
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


class ObsAodMetOffice : public ObsOperatorBase,
                   private util::ObjectCounter<ObsAodMetOffice> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsAodMetOfficeParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() {return "ufo::ObsAodMetOffice";}

  // -----------------------------------------------------------------------------
  /*! \brief A dust-only AOD forward model operator at one wavelength
  *
  * \details Calculates the dust Aerosol Optical Depth at one wavelength from the dust mass concentration, for a variable number
  * of dust size bins. The Met Office forecast model has 2 bins operationally but it can also be run with 6 bins.
  *
  * The AOD is calculated as:
  * \f[ AOD = \sum_{d=1}^{NBins}\sum_{k=1}^{Nlevels}(1/g) k_{ext}(d) r(k) |\Delta P(k)| \f],
  *
  * where: \f$g\f$ is the gravity constant, \f$k_{ext}(d)\f$ is the extinction coefficient per dust bin \f$d\f$, \f$r(k)\f$ is the dust
  * mass concentration per atmospheric level \f$k\f$, and the summation is performed over all bins (\f$NBins\f$) and all levels
  * (\f$Nlevels\f$). The number of bins and the extinction coefficients per bin are read in from yaml files and other inputs are
  * taken from the geovals. Note that the geoval input of pressure needs to be on staggered levels relative to mass concentration
  * so that \f$\|\Delta P|\f$ is at the same height as mass concentration.
  *
  * \date Aug. 2021: Created by H. Lawrence (Met Office)
  *
  */
  // -----------------------------------------------------------------------------
  ObsAodMetOffice(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAodMetOffice();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t&) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

 private:
  void print(std::ostream &) const override;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;  // list of the required geovals
  std::size_t NDustBins_;  // Number of dust bins
  std::vector<double> AodKExt_;  // Extinction coefficients per bin, independent of humidity
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_AEROSOLS_AODMETOFFICE_OBSAODMETOFFICE_H_
