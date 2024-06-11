/*
 *
 * Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/aerosols/AODMetOffice/ObsAodMetOffice.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAodMetOffice> makerAodMetOffice_("AodMetOffice");
// -----------------------------------------------------------------------------

ObsAodMetOffice::ObsAodMetOffice(const ioda::ObsSpace & odb, const Parameters_ &params)
  : ObsOperatorBase(odb), odb_(odb), varin_()
{
  // 1. Read in number of dust bins and extinction coefficients from yaml file.
  NDustBins_ = params.NDustBins;  // number of dust bins
  AodKExt_ = params.AodKExt;  // Extinction coefficient for each dust bin

  // check the size of AodKExt_ is consistent with NDustBins_
  ASSERT(AodKExt_.size() == NDustBins_);

  // check the number of dust bins is supported
  if (NDustBins_ < 1) {
    // raise error as we need to have some dust bins
    throw eckit::UserError("NDustBins_ must be 1 or higher", Here());
  }

  // 2. Define model fields needed for AOD calculation
  // Note the dust model fields required depends on number of bins, defined in the yaml file.
  // 1 - 6 bins are supported.

  // pressure
  varin_.push_back("air_pressure_levels");
  varin_.push_back("surface_pressure");

  // Dust model fields
  std::string dust_var_name;  // dust variable name in geovals
  for (size_t d = 0; d < NDustBins_; d++) {
      dust_var_name = "mass_fraction_of_dust00" + std::to_string(d+1) + "_in_air";
      varin_.push_back(dust_var_name);
  }

  oops::Log::trace() << "ObsAodMetOffice created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodMetOffice::~ObsAodMetOffice() {
  oops::Log::trace() << "ObsAodMetOffice destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodMetOffice::simulateObs(const GeoVaLs & geovals, ioda::ObsVector & hofx,
                                  ObsDiagnostics &, const QCFlags_t & qc_flags) const {
  // Calculate dust AOD at one wavelength (e.g. 550nm) from dust mass concentration.

  // Check hofx size and initialise hofx to zero:
  ASSERT(geovals.nlocs() == hofx.nlocs());
  hofx.zero();

  // Get number of obs locations & number of model pressure levels:
  std::size_t nprofiles = geovals.nlocs();
  // number of full (rho) levels
  std::size_t nlevels   = geovals.nlevs(oops::Variable{"air_pressure_levels"});

  // Get 2-D surface pressure
  std::vector<double> ps(nprofiles);  // surface pressure (Pa)
  geovals.get(ps, oops::Variable{"surface_pressure"});

  // Get 3-D air pressure on rho levels (Pa), one level at a time
  std::vector<std::vector<double>> plev(nlevels, std::vector<double>(nprofiles));
  for (std::size_t i = 0; i < nlevels; ++i) {
        geovals.getAtLevel(plev[i], oops::Variable{"air_pressure_levels"}, i);
        }

  // check model fields are from top-down, fail if not
  if (plev.front() > plev.back()) {
    throw eckit::BadValue("model fields must be ordered from the top down", Here());
  }

  double alpha;
  std::vector<double> mass(nlevels - 1);  // mass concentration on half (theta) levels

  // loop over profiles
  for (size_t p = 0; p < nprofiles; p++) {
    // loop over dust bins
    for (size_t d = 0; d < NDustBins_; d++) {
      // get mass concentration for dust bin d
      oops::Variable dust_var{"mass_fraction_of_dust00" + std::to_string(d+1) + "_in_air"};
      geovals.getAtLocation(mass, dust_var, p);
      // start with lowest model layer, using surface pressure
      // NB this assumes the first mass layer is the lowest layer just above the surface
      alpha = (1.0 / Constants::grav) * (ps[p] - plev[nlevels - 2][p]) * AodKExt_[d];
      hofx[p] += mass[nlevels - 2] * alpha;
      // loop over the rest of the model layers
      for (size_t k = 0; k < (nlevels - 2); k++) {
        alpha = (1.0 / Constants::grav) * (plev[k+1][p] - plev[k][p]) * AodKExt_[d];
        hofx[p] += mass[k]*alpha;
      }
    }
  }

  oops::Log::trace() << "ObsAodMetOffice: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodMetOffice::print(std::ostream & os) const {
  os << "ObsAodMetOffice::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
