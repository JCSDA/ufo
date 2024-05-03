/*
 *
 * Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_AEROSOLS_AODMETOFFICE_OBSAODMETOFFICETLAD_H_
#define UFO_OPERATORS_AEROSOLS_AODMETOFFICE_OBSAODMETOFFICETLAD_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/aerosols/AODMetOffice/ObsAodMetOfficeParameters.h"

// Forward declarations
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
class ObsAodMetOfficeTLAD : public LinearObsOperatorBase,
                       private util::ObjectCounter<ObsAodMetOfficeTLAD> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ioda::ObsDataVector<int> QCFlags_t;

  typedef ObsAodMetOfficeParameters Parameters_;

  static const std::string classname() {return "ufo::ObsAodMetOfficeTLAD";}

  ObsAodMetOfficeTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAodMetOfficeTLAD();

  // Obs Operators
  // -----------------------------------------------------------------------------
  /*! \brief Compute Jacobian matrix d(AOD)/d(mass concentration).
  *
  * \details This matrix is of size: number dust bins x number levels x number profiles.
  * This method must be called before calling the TL/AD methods.
  *
  * \date Oct. 2021: Created by H. Lawrence (Met Office)
  */
  // -----------------------------------------------------------------------------
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;

  // -----------------------------------------------------------------------------
  /*! \brief Given an increment to the model state (dust mass concentration), calculate
  *  an increment to the observation (AOD).
  *
  * \date Oct. 2021: Created by H. Lawrence (Met Office)
  */
  // -----------------------------------------------------------------------------
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;

  // -----------------------------------------------------------------------------
  /*! \brief Given an increment to the observation (AOD), calculate the equivalent
  *   increment to the model state (dust mass concentration).
  *
  * \date Oct. 2021: Created by H. Lawrence (Met Office)
  */
  // -----------------------------------------------------------------------------
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;


  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

 private:
  void print(std::ostream &) const override;
  oops::Variables varin_;
  std::vector<std::vector<std::vector<double>>> kMatrix_;
  bool trajInit_;
  std::size_t NDustBins_;  // Number of dust bins
  std::vector<double> AodKExt_;  // Extinction coefficients per bin, independent of humidity
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_AEROSOLS_AODMETOFFICE_OBSAODMETOFFICETLAD_H_
