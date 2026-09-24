/*
 *
 * Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>
#include <vector>

#include "ufo/operators/aerosols/AODMassFraction/ObsAodMassFraction.h"
#include "ufo/operators/aerosols/AODMassFraction/ObsAodMassFractionTLAD.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAodMassFractionTLAD> makerAodMassFractionTL_("AodMassFraction");
// -----------------------------------------------------------------------------

ObsAodMassFractionTLAD::ObsAodMassFractionTLAD(const ioda::ObsSpace &odb, const Parameters_ &params)
    : LinearObsOperatorBase(odb), varin_(), kMatrix_(), trajInit_(false), nspecies_(0),
      params_(params),
      obsOperator_(std::make_unique<ObsAodMassFraction>(odb, params)) {

  // 1. Initialize aerosol species/modes and corresponding variable names for mass and
  // number fractions dust variables name in geovals.
  // By default, the list of variables is empty for number fractions for models
  // that do not use them.
  varaerosolnumber_.clear();

  // Define model-specific variables and parameters
  switch (params_.aerosolModelName.value()) {
    case aerosolModel::UKCA_Dust:

      oops::Log::trace() << "ObsAodMassFraction: using UKCA dust aerosol model." << std::endl;
      // aerosol fields for UKCA dust model
      speciesList_ = {"coarse", "accumulation"};
      varaerosolmass_ = {"mass_fraction_of_dust_coarse_aerosol_particles_in_air",
                     "mass_fraction_of_dust_accumulation_aerosol_particles_in_air"};
      varaerosolnumber_ = {"number_fraction_of_coarse_aerosol_particles_in_air",
                     "number_fraction_of_accumulation_aerosol_particles_in_air"};
      break;
    default:
      throw eckit::UserError("AODMassFraction aerosol model not recognised", Here());
  }

  // speciesList_, varaerosolmass_, and varaerosolnumber_ are now set based on the aerosol model
  nspecies_ = speciesList_.size();  // number of species/modes to include in AOD calc

  // check the number of dust bins/modes/species is supported
  if (nspecies_ < 1) {
    // raise error as we need to have some dust bins
    throw eckit::UserError("nspecies must be 1 or higher", Here());
  }

  // Define which model fields will be needed
  // pressure fields:
  varin_.push_back("air_pressure_levels");
  varin_.push_back("air_pressure_at_surface");

  for (const auto &var : varaerosolmass_) {
    varin_.push_back(var);
  }
  for (const auto &var : varaerosolnumber_) {
    varin_.push_back(var);
  }

  // Get parameters needed for extinction coefficient calculation method
  switch (params_.extinctionMethodName.value()) {
  case extinctionMethod::MetOfficeLUTFit:
    // Best fit parameters for UKCA dust extinction coefficient calculation
    if (params_.coarseModeParams.value()) {
      coarseParams_ = *params_.coarseModeParams.value();
    } else {
      throw eckit::UserError("Coarse mode extinction coefficient fit parameters "
                             "must be provided for MetOfficeLUTFit",
                             Here());
    }
    // Check that the correct number of fit parameters have been provided.
    if (coarseParams_.size() != 4) {
      throw eckit::UserError("Coarse mode extinction coefficient calculation requires "
                             "4 fit parameters",
                             Here());
    }
    if (params_.accumulationModeParams.value()) {
      accumulationParams_ = *params_.accumulationModeParams.value();
    } else {
      throw eckit::UserError("Accumulation mode extinction coefficient fit parameters "
                             "must be provided for MetOfficeLUTFit",
                             Here());
    }
    if (accumulationParams_.size() != 7) {
      throw eckit::UserError("Accumulation mode extinction coefficient calculation requires "
                             "7 fit parameters",
                             Here());
    }
    break;
  default:
    throw eckit::UserError("Extinction coefficient calculation method not recognised", Here());
  }

  // calculate dAOD for each species, level and profile

  oops::Log::trace() << "ObsAodMassFractionTLAD constructor done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodMassFractionTLAD::~ObsAodMassFractionTLAD() {
  trajInit_ = false;
  oops::Log::trace() << "ObsAodMassFractionTLAD destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodMassFractionTLAD::setTrajectory(const GeoVaLs &geovals, ObsDiagnostics &,
                                           const QCFlags_t &qc_flags) {
  oops::Log::trace() << "ObsAodMassFractionTLAD::setTrajectory start" << std::endl;
  // Get number of obs locations & number of model pressure levels:
  std::size_t nprofiles = geovals.nlocs();
  // number of full (rho) levels
  std::size_t nlevels = geovals.nlevs(oops::Variable{"air_pressure_levels"});

  // Check whether the model has number fraction fields for the aerosol species/modes.
  // If not, then the TL/AD will only consider the mass fraction contributions.
  const bool hasNumberFraction = !varaerosolnumber_.empty();

  // (Re)initialise the Jacobian matrix
  // Note: kMatrix_ is a 4-D vector of size: number species x number levels x number profiles x 2
  // The last dimension is for the mass fraction [0] and number fraction [1]; if number fraction
  // is not used in the model, this column will remain set to 0.
  if (!kMatrix_.empty()) kMatrix_.clear();

  kMatrix_.assign(
      nspecies_,
      std::vector<std::vector<std::vector<double>>>(
          nlevels - 1,
          std::vector<std::vector<double>>(
              nprofiles,
              std::vector<double>(2, 0.0))));

  // Get 2-D surface pressure
  std::vector<double> ps(nprofiles);  // surface pressure (Pa)
  geovals.get(ps, oops::Variable{"air_pressure_at_surface"});

  // Get 3-D air pressure on rho levels (Pa), one level at a time
  std::vector<double> plev(nlevels);
  geovals.getAtLocation(plev, oops::Variable{"air_pressure_levels"}, 0);

  // check model fields are ordered from top down, fail if not
  if (plev.front() > plev.back()) {
    throw eckit::BadValue("model fields must be ordered from top down", Here());
  }

  // Function to calculate the factor: vertical difference in pressure
  // divided by gravitational constant
  auto beta = [](double deltaP) {
        return (1.0 / Constants::grav) * deltaP;
  };

  for (std::size_t jloc = 0; jloc < nprofiles; ++jloc) {
    // get the air pressure column geovals at this location
    geovals.getAtLocation(plev, oops::Variable{"air_pressure_levels"}, jloc);

    for (std::size_t jmode = 0; jmode < nspecies_; ++jmode) {
      std::vector<double> massfrac(nlevels - 1, 0.0);  // initialise to zero
      std::vector<double> numfrac(nlevels - 1, 0.0);  // initialise to zero
      double gamma0 = 0.0;
      double gamma1 = 0.0;

      const oops::Variable massFraction{varaerosolmass_[jmode]};
      geovals.getAtLocation(massfrac, massFraction, jloc);
      if (hasNumberFraction) {
        const oops::Variable numberFraction{varaerosolnumber_[jmode]};
        geovals.getAtLocation(numfrac, numberFraction, jloc);
      }

      for (std::size_t k = 0; k < (nlevels - 1); ++k) {
        if (massfrac[k] == 0) {
            massfrac[k] = 1.0e-12;  // set to small value to avoid division by zero
        }
        if (hasNumberFraction && numfrac[k] == 0) {
            numfrac[k] = 1.0e-22;  // set to small value to avoid division by zero
        }
      }

      for (std::size_t k = 0; k < (nlevels - 2); ++k) {
        const double diameter =
            obsOperator_->getUKCADustDiameter(numfrac[k], massfrac[k], speciesList_[jmode]);

        gamma0 = obsOperator_->getUKCADustKExt(diameter, speciesList_[jmode]) +
                 massfrac[k] * getUKCAExtinctionDerivatives(
                                     diameter,
                                     speciesList_[jmode]) *
                               getUKCADiameterDerivatives(numfrac[k], massfrac[k],
                                                          speciesList_[jmode], "mass");

        gamma1 = massfrac[k] * getUKCAExtinctionDerivatives(
                                     diameter,
                                     speciesList_[jmode]) *
                               getUKCADiameterDerivatives(numfrac[k], massfrac[k],
                                                          speciesList_[jmode], "number");

        kMatrix_[jmode][k][jloc][0] = beta(plev[k + 1] - plev[k]) * gamma0;
        kMatrix_[jmode][k][jloc][1] = beta(plev[k + 1] - plev[k]) * gamma1;
      }

      const double diameter = obsOperator_->getUKCADustDiameter(
                              numfrac[nlevels - 2], massfrac[nlevels - 2], speciesList_[jmode]);

      gamma0 = obsOperator_->getUKCADustKExt(diameter, speciesList_[jmode]) +
               massfrac[nlevels - 2] * getUKCAExtinctionDerivatives(
                                             diameter,
                                             speciesList_[jmode]) *
                                       getUKCADiameterDerivatives(numfrac[nlevels - 2],
                                          massfrac[nlevels - 2],
                                          speciesList_[jmode], "mass");
      gamma1 = massfrac[nlevels - 2] * getUKCAExtinctionDerivatives(
                                             diameter,
                                             speciesList_[jmode]) *
                                       getUKCADiameterDerivatives(numfrac[nlevels - 2],
                                                                  massfrac[nlevels - 2],
                                                                  speciesList_[jmode], "number");

      kMatrix_[jmode][nlevels - 2][jloc][0] = beta(ps[jloc] - plev[nlevels - 2]) * gamma0;
      kMatrix_[jmode][nlevels - 2][jloc][1] = beta(ps[jloc] - plev[nlevels - 2]) * gamma1;
    }  // end loop over species/modes
  }  // end loop over profiles

  trajInit_ = true;

  oops::Log::trace() << "ObsAodMassFractionTLAD::setTrajectory done" << std::endl;
}

double ObsAodMassFractionTLAD::getUKCADiameterDerivatives(const double num_in, const double mass_in,
                                                      const std::string &mode,
                                                      const std::string &var) const {
  double num = num_in;
  double mass = mass_in;

  // Calculate derivatives df dust diameter with respect to both mass
  // and number fractions.
  const double sdev_a = 1.59;  // accumulation mode unitless modal width
  const double sdev_c = 2;     // coarse mode unitless modal width
  const double rho = 2650;     // assumed density of dust particles (kg/m^3)
  const double sdev_D = (mode == "accumulation") ? sdev_a : sdev_c;
  const double c_a = ((6.0 * Constants::k_B) / (M_PI * Constants::rd * rho)) *
                     std::exp(-4.5 * std::pow(std::log(sdev_D), 2));

  const double dD_dn = (-1.0 / 3.0) * std::pow(c_a, 1.0 / 3.0) * std::pow(mass, 1.0 / 3.0) *
                       std::pow(num, -4.0 / 3.0);

  const double dD_dm = (1.0 / 3.0) * std::pow(c_a, 1.0 / 3.0) * std::pow(mass, -2.0 / 3.0) *
                       std::pow(num, -1.0 / 3.0);

  return (var == "mass") ? dD_dm : dD_dn;
}
// -----------------------------------------------------------------------------

double ObsAodMassFractionTLAD::getUKCAExtinctionDerivatives(const double D,
                                                            const std::string &mode) const {
  // D must be positive
  if (D <= 0) {
    throw eckit::UserError("Diameter must be positive for extinction derivative calculation",
                           Here());
  }
  double dk_dD;

  if (mode == "coarse") {
    // Coarse mode
    const double c_1 = coarseParams_[0] * std::exp(coarseParams_[2] * D);
    const double c_2 = coarseParams_[1] * std::pow(D, coarseParams_[1] - 1.0) +
                       coarseParams_[2] * std::pow(D, coarseParams_[1]);

    dk_dD = c_1 * c_2;
  } else if (mode == "accumulation") {
    // Accumulation mode
    const double C = accumulationParams_[6] * std::pow(std::log(D), 6) +
                     accumulationParams_[5] * std::pow(std::log(D), 5) +
                     accumulationParams_[4] * std::pow(std::log(D), 4) +
                     accumulationParams_[3] * std::pow(std::log(D), 3) +
                     accumulationParams_[2] * std::pow(std::log(D), 2) +
                     accumulationParams_[1] * std::log(D) + accumulationParams_[0];
    const double dC = 6.0 * accumulationParams_[6] * std::pow(std::log(D), 5) +
                      5.0 * accumulationParams_[5] * std::pow(std::log(D), 4) +
                      4.0 * accumulationParams_[4] * std::pow(std::log(D), 3) +
                      3.0 * accumulationParams_[3] * std::pow(std::log(D), 2) +
                      2.0 * accumulationParams_[2] * std::log(D) + accumulationParams_[1];
    dk_dD = (1.0 / D) * dC * std::exp(C);
  } else {
    throw eckit::UserError("Aerosol mode not recognised for extinction derivative calculation",
                           Here());
  }

  return dk_dD;
}

// -----------------------------------------------------------------------------

void ObsAodMassFractionTLAD::simulateObsTL(const GeoVaLs &geovals, ioda::ObsVector &hofx) const {
  // Ensure trajectory has already been calculated
  ASSERT(trajInit_);
  ASSERT(geovals.nlocs() == hofx.nlocs());

  const bool hasNumberFraction = !varaerosolnumber_.empty();

  const std::size_t nprofiles = geovals.nlocs();
  const std::size_t nlevels = geovals.nlevs(oops::Variable{"air_pressure_levels"});

  hofx.zero();

  for (std::size_t jloc = 0; jloc < nprofiles; ++jloc) {
    for (std::size_t jmode = 0; jmode < nspecies_; ++jmode) {
      std::vector<double> massfrac_d(nlevels - 1, 0.0);
      std::vector<double> numfrac_d(nlevels - 1, 0.0);

      const oops::Variable massFraction{varaerosolmass_[jmode]};
      geovals.getAtLocation(massfrac_d, massFraction, jloc);

      if (hasNumberFraction) {
        const oops::Variable numberFraction{varaerosolnumber_[jmode]};
        geovals.getAtLocation(numfrac_d, numberFraction, jloc);
      }

      for (std::size_t k = 0; k < (nlevels - 1); ++k) {
        hofx[jloc] += kMatrix_[jmode][k][jloc][0] * massfrac_d[k];
        if (hasNumberFraction) {
          hofx[jloc] += kMatrix_[jmode][k][jloc][1] * numfrac_d[k];
        }
      }
    }
  }
  oops::Log::trace() << "ObsAodMassFractionTLAD::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodMassFractionTLAD::simulateObsAD(GeoVaLs &geovals, const ioda::ObsVector &hofx) const {
  // Ensure trajectory has already been calculated
  ASSERT(trajInit_);

  const bool hasNumberFraction = !varaerosolnumber_.empty();

  // Get number of obs locations & number of model pressure levels:
  std::size_t nprofiles = geovals.nlocs();
  // number of full (rho) levels
  std::size_t nlevels = geovals.nlevs(oops::Variable{"air_pressure_levels"});

  // Check hofx size
  ASSERT(geovals.nlocs() == hofx.nlocs());

  // Get the missing value indicator
  const double missing = util::missingValue<double>();

  // Loop through the obs, adding the increment to the model state
  for (std::size_t jloc = 0; jloc < nprofiles; jloc++) {
    if (hofx[jloc] != missing) {
      // loop over the aerosol species/modes (nspecies)
      for (std::size_t jmode = 0; jmode < nspecies_; jmode++) {
        std::vector<double> massfrac_d(nlevels - 1, 0.0);
        std::vector<double> numfrac_d(nlevels - 1, 0.0);

        const oops::Variable massFraction{varaerosolmass_[jmode]};
        geovals.getAtLocation(massfrac_d, massFraction, jloc);
        if (hasNumberFraction) {
          const oops::Variable numberFraction{varaerosolnumber_[jmode]};
          geovals.getAtLocation(numfrac_d, numberFraction, jloc);
        }

        // loop over the model layers
        for (std::size_t k = 0; k < (nlevels - 1); ++k) {
          massfrac_d[k] += kMatrix_[jmode][k][jloc][0] * hofx[jloc];
          if (hasNumberFraction) {
            numfrac_d[k]  += kMatrix_[jmode][k][jloc][1] * hofx[jloc];
          }
        }

        // Store the updated model state increments
        geovals.putAtLocation(massfrac_d, massFraction, jloc);
        if (hasNumberFraction) {
          const oops::Variable numberFraction{varaerosolnumber_[jmode]};
          geovals.putAtLocation(numfrac_d, numberFraction, jloc);
        }
      }
    }
  }
  oops::Log::trace() << "ObsAodMassFractionTLAD::simulateObsAD complete" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodMassFractionTLAD::print(std::ostream &os) const {
  os << "ObsAodMassFractionTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
