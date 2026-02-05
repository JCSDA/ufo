/*
 *
 * Crown Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <math.h>
#include <ostream>
#include <string>
#include <vector>

#include "ufo/operators/aerosols/AODMassFraction/ObsAodMassFraction.h"

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAodMassFraction> makerAodMassFraction_("AodMassFraction");
// -----------------------------------------------------------------------------

ObsAodMassFraction::ObsAodMassFraction(const ioda::ObsSpace & odb, const Parameters_ &params)
  : ObsOperatorBase(odb), odb_(odb), params_(params), varin_()
{
  oops::Log::trace() << "ObsAodMassFraction constructor starting" << std::endl;

  // 1. Define the aerosol model fields needed for AOD calculation,
  //    depending on the model.

  // By default the list of variables for number fractions is empty (for models that do not need it)
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

  // 2. Get parameters needed for extinction coefficient calculation method
  switch (params_.extinctionMethodName.value()) {
    case extinctionMethod::MetOfficeLUTFit:
    // Best fit parameters for UKCA dust extinction coefficient calculation
    if (params_.coarseModeParams.value()) {
      coarseParams_ = *params_.coarseModeParams.value();
    } else {
      throw eckit::UserError("Coarse mode extinction coefficient fit parameters "
        "must be provided for MetOfficeLUTFit", Here());
    }
    if (params_.accumulationModeParams.value()) {
      accumulationParams_ = *params_.accumulationModeParams.value();
    } else {
      throw eckit::UserError("Accumulation mode extinction coefficient fit parameters "
        "must be provided for MetOfficeLUTFit", Here());
    }
    // Check that the correct number of fit parameters have been provided.
    if (coarseParams_.size() != 4) {
        throw eckit::UserError("Coarse mode extinction coefficient calculation requires "
          "4 fit parameters", Here());
    }
    if (accumulationParams_.size() != 7) {
        throw eckit::UserError("Accumulation mode extinction coefficient calculation requires "
          "7 fit parameters", Here());
    }
      break;
    default:
      throw eckit::UserError("Extinction coefficient calculation method not recognised", Here());
  }

  // 3. Define the meteorological model fields needed for AOD calculation and add to varin_.
  varatm_ = {"air_pressure_levels", "air_pressure_at_surface"};
  for (const auto & var : varatm_) {
      varin_.push_back(var);
  }

  // 4. Add aerosol variables to varin_
  for (const auto & var : varaerosolmass_) {
          varin_.push_back(var);
      }
  if (!varaerosolnumber_.empty()) {
    for (const auto & var : varaerosolnumber_) {
            varin_.push_back(var);
        }
    }

  oops::Log::trace() << "ObsAodMassFraction created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodMassFraction::~ObsAodMassFraction() {
  oops::Log::trace() << "ObsAodMassFraction destructed" << std::endl;
}

// -----------------------------------------------------------------------------

double ObsAodMassFraction::getUKCADustKExt(const double diameter, const std::string & mode) const {
  // Returns the extinction coefficient for a given median diameter, for either accumulation or
  // coarse dust insoluble modes in the UKCA/GLOMAP dust model. Two best fit functions, one for
  // each mode, are used to derive the extinction coefficient, based on RADAER LUTs (Met
  // Office Radiation scheme).
  double AODKExt = 0.0;

  if (mode == "accumulation") {
    // Accumulation mode diameter best fit function.
    const double ln_x = std::log(diameter);
    const double exponent = accumulationParams_[6]*std::pow(ln_x, 6) +
                            accumulationParams_[5]*std::pow(ln_x, 5) +
                            accumulationParams_[4]*std::pow(ln_x, 4) +
                            accumulationParams_[3]*std::pow(ln_x, 3) +
                            accumulationParams_[2]*std::pow(ln_x, 2) +
                            accumulationParams_[1]*ln_x + accumulationParams_[0];
    AODKExt = std::exp(exponent);

  } else if (mode == "coarse") {
    // Coarse mode diameter best fit function.
    AODKExt = coarseParams_[0] * std::pow(diameter, coarseParams_[1]) *
              std::exp(coarseParams_[2] * diameter) + coarseParams_[3];
  }

  oops::Log::trace() << "ObsAodMassFraction:getUKCADustKExt ran for " << mode <<
                        " mode in ObsAodMassFraction::getUKCADustKExt: " << std::endl;
  return AODKExt;
}

double ObsAodMassFraction::getUKCADustDiameter(const double num, const double mass,
                               const std::string & mode) const {
  // Calculate UKCA dust modal diameter (per height) as a function of the mass
  // and number fractions.
  const double sdev_a = 1.59;  // accumulation mode unitless modal width
  const double sdev_c = 2;  // coarse mode unitless modal width
  const double rho = 2650;  // assumed density of dust particles (kg/m^3)

  const double c_a = (6.0 *Constants::k_B) / (M_PI * Constants::rd * rho*
                std::exp(4.5 * std::pow(std::log(sdev_a), 2)));
  const double c_c = (6.0 *Constants::k_B) / (M_PI * Constants::rd * rho*
                std::exp(4.5 * std::pow(std::log(sdev_c), 2)));
  const double c = (mode == "accumulation") ? c_a : c_c;

  // Return the diameter (m) if mass and number fractions are non-zero and positive.
  // Else return zero.
  return (num > 0 && mass > 0) ? std::pow(c*(mass / num), 1.0 / 3.0) : 0.0;
}

// -----------------------------------------------------------------------------

void ObsAodMassFraction::simulateObs(const GeoVaLs & geovals, ioda::ObsVector & hofx,
                         ObsDiagnostics &, const QCFlags_t & qc_flags) const {
  // Calculate aerosol AOD at one wavelength (e.g. 550nm)

  // Check hofx size and initialise hofx to zero:
  ASSERT(geovals.nlocs() == hofx.nlocs());
  hofx.zero();

  // Get number of obs locations, number of full (rho) model pressure levels,
  // number of aerosol species/modes
  const std::size_t nlocs = geovals.nlocs();
  const std::size_t nlevels = geovals.nlevs(oops::Variable{"air_pressure_levels"});
  std::size_t nspecies = speciesList_.size();  // number of species/modes to include in AOD calc

  // Get 2-D surface pressure
  std::vector<double> ps(nlocs);  // surface pressure (Pa)
  geovals.get(ps, oops::Variable{"air_pressure_at_surface"});

  // Function to calculate alpha, the factor to multiply mass fraction by in order
  // to get the AOD contribution
  auto alpha = [](double deltaP, double AODKExt) {
        return (1.0 / Constants::grav) * deltaP * AODKExt;
      };

  // Define variables needed for AOD calculation
  std::vector<double> plev(nlevels);  // air pressure on rho levels (Pa)
  std::vector<double> mass(nlevels - 1);  // mass fraction on half (theta) levels
  std::vector<double> num(nlevels - 1);  // number fraction on half (theta) levels
  const std::size_t lowest_level = nlevels - 2;  // index of lowest model layer (just above surface)

  // loop over profiles and species and calculate AOD
  for (std::size_t jloc = 0; jloc < nlocs; jloc++) {
    for (std::size_t jmode = 0; jmode < nspecies; jmode++) {
      // get the required 3-D geovals at this location
      geovals.getAtLocation(plev, oops::Variable{"air_pressure_levels"}, jloc);
      const oops::Variable massFraction{varaerosolmass_[jmode]};
      geovals.getAtLocation(mass, massFraction, jloc);
      if (!varaerosolnumber_.empty()) {
        const oops::Variable numberFraction{varaerosolnumber_[jmode]};
        geovals.getAtLocation(num, numberFraction, jloc);
      }

      // calculate AOD for lowest model layer, first getting the extinction coefficient
      double AODKExt = 0.0;
      if (params_.extinctionMethodName.value() == extinctionMethod::MetOfficeLUTFit) {
        const double diameter = getUKCADustDiameter(num[lowest_level], mass[lowest_level],
                                            speciesList_[jmode]);
        if (diameter > 0) {
          AODKExt = getUKCADustKExt(diameter, speciesList_[jmode]);
        }
      }
      hofx[jloc] += mass[lowest_level] * alpha(ps[jloc] - plev[lowest_level], AODKExt);

      // calculate AOD for the rest of the model layers
      for (std::size_t jlev = 0; jlev < (nlevels - 2); jlev++) {
        AODKExt = 0.0;
        if (params_.extinctionMethodName.value() == extinctionMethod::MetOfficeLUTFit) {
          auto diameter = getUKCADustDiameter(num[jlev], mass[jlev], speciesList_[jmode]);
          if (diameter > 0) {
            AODKExt = getUKCADustKExt(diameter, speciesList_[jmode]);
          }
        }
        hofx[jloc] += mass[jlev] * alpha(plev[jlev + 1] - plev[jlev], AODKExt);
      }
    }
  }

  oops::Log::trace() << "ObsAodMassFraction: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodMassFraction::print(std::ostream & os) const {
  os << "ObsAodMassFraction::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
