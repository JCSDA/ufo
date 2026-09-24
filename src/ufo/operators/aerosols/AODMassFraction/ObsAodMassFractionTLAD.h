/*
 *
 * Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_AEROSOLS_AODMASSFRACTION_OBSAODMASSFRACTIONTLAD_H_
#define UFO_OPERATORS_AEROSOLS_AODMASSFRACTION_OBSAODMASSFRACTIONTLAD_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ObsAodMassFractionParameters.h"

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
  class ObsAodMassFraction;


// -----------------------------------------------------------------------------
class ObsAodMassFractionTLAD : public LinearObsOperatorBase,
                       private util::ObjectCounter<ObsAodMassFractionTLAD> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  using QCFlags_t = LinearObsOperatorBase::QCFlags_t;

  typedef ObsAodMassFractionParameters Parameters_;

  static const std::string classname() {return "ufo::ObsAodMassFractionTLAD";}

  ObsAodMassFractionTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAodMassFractionTLAD();


  // Obs Operators
  // -----------------------------------------------------------------------------
  /*! \brief Compute Jacobian matrix d(AOD)/d(mass concentration) and optionally
  *   d(AOD)/d(number fraction).
  *
  * \details This matrix is of size: number species x number levels x number profiles.
  *  This method must be called before calling the TL/AD methods.
  *
  * \date May 2026: Created by C. Charlton-Perez (Met Office)
  */
  // -----------------------------------------------------------------------------
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;

  // -----------------------------------------------------------------------------
  /*! \brief Given an increment to the model state (species mass concentration
  *  and optionally number concentration);
  *  calculate an increment to the observation (AOD).
  *
  * \date May 2026: Created by C. Charlton-Perez (Met Office)
  */
  // -----------------------------------------------------------------------------
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;

  // -----------------------------------------------------------------------------
  /*! \brief Given an increment to the observation (AOD), calculate the equivalent
  *   increment to the model state (mass concentration, number fraction).
  *
  * \date May 2026: Created by C. Charlton-Perez (Met Office)
  */
  // -----------------------------------------------------------------------------
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

 private:
  void print(std::ostream &) const override;
  oops::Variables varin_;  // list of all required geovals
  std::vector<std::vector<std::vector<std::vector<double>>>> kMatrix_;
  bool trajInit_;
  std::size_t nspecies_;  // Number of species/modes
  const Parameters_ params_;
  std::vector<std::string> varatm_;  // list of required atmospheric geovals
  std::vector<std::string> speciesList_;  // list of aerosol species / modes
  std::vector<std::string> varaerosolmass_;  // list of mass fraction variables
                                             // for aerosol species
  std::vector<std::string> varaerosolnumber_;  // list of number fraction variables
                                               // for aerosol species
  std::vector<double> coarseParams_;  // UKCA dust coarse mode extinction coefficient
                                      // fit parameters
  std::vector<double> accumulationParams_;  // UKCA dust accumulation mode extinction
                                            // coefficient fit parameters
  double getUKCAExtinctionDerivatives(const double r,
                              const std::string & mode) const;
  double getUKCADiameterDerivatives(const double num, const double mass, const std::string & mode,
                              const std::string & var) const;
  std::unique_ptr<ObsAodMassFraction> obsOperator_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_AEROSOLS_AODMASSFRACTION_OBSAODMASSFRACTIONTLAD_H_
