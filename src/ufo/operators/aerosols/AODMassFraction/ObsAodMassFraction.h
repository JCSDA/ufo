/*
 *
 * Crown Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_AEROSOLS_AODMASSFRACTION_OBSAODMASSFRACTION_H_
#define UFO_OPERATORS_AEROSOLS_AODMASSFRACTION_OBSAODMASSFRACTION_H_

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/operators/aerosols/AODMassFraction/ObsAodMassFractionParameters.h"

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

class ObsAodMassFraction : public ObsOperatorBase,
                   private util::ObjectCounter<ObsAodMassFraction> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsAodMassFractionParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() {return "ufo::ObsAodMassFraction";}

  // -----------------------------------------------------------------------------
  /*! \brief Forward model operator to calculate aerosol AOD at one wavelength.
  *
  * \details This operator calculates the Aerosol Optical Depth at one wavelength from the mass
  * fractions of a list of aerosol species, defined for different aerosol models. Currently this is 
  * implemented for the Met Office UKCA/GLOMAP dust model only but other models could be added in 
  * future.
  *
  * The calculation is performed as follows:
  *
  * \f[ AOD = \sum_{d=1}^{NSpecies}\sum_{k=1}^{Nlevels}(1/g) k_{ext}(d)(k) r(k) |\Delta P(k)| \f],
  *
  * where: \f$g\f$ is the gravity constant, \f$k_{ext}(d)(k)\f$ and \f$r(k)\f$ are respectively the 
  * extinction coefficient and mass fraction per aerosol species \f$d\f$ and atmospheric level \f$k\f$,
  * and the summation is performed over all aerosol species (\f$NSpecies\f$) and all levels 
  * (\f$Nlevels\f$).
  *
  * The extinction coefficient can be calculated using one of a number of methods, specified in the yaml
  * configuration file. Currently only the "MetOfficeLUTFit" method is implemented, which uses best 
  * fit functions, fit to the Met Office RADAER LUTs. A pure LUT approach could be implemented in future.
  *
  * There is also the option to include the number fraction of aerosol particles in the calculation, by
  * defining these variables in ObsAodMassFraction.cc, as done for the UKCA_Dust model. These variables are required 
  * for modal models, such as the Met Office UKCA dust model, in order to calculate the modal diameter 
  * and then extinction coefficient. However if this list of number fraction variables is not specified 
  * for a given model in ObsAodMassFraction.cc then the number fraction variables list is left empty by default and 
  * not used.
  *
  * \date Dec. 2025: Created by H. Lawrence and P. Siwek (Met Office)
  *
  */
  // -----------------------------------------------------------------------------
  ObsAodMassFraction(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAodMassFraction();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t&) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

 private:
  void print(std::ostream &) const override;
  const ioda::ObsSpace& odb_;
  const Parameters_ params_;
  oops::Variables varin_;  // list of all required geovals
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
  double getUKCADustDiameter(const double num, const double mass, const std::string & mode) const;
  double getUKCADustKExt(const double diameter, const std::string & mode) const;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_AEROSOLS_AODMASSFRACTION_OBSAODMASSFRACTION_H_
