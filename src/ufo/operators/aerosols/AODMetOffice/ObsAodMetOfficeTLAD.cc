/*
 *
 * Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/aerosols/AODMetOffice/ObsAodMetOfficeTLAD.h"

#include <ostream>
#include <string>
#include <vector>

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
static LinearObsOperatorMaker<ObsAodMetOfficeTLAD> makerAodMetOfficeTL_("AodMetOffice");
// -----------------------------------------------------------------------------

ObsAodMetOfficeTLAD::ObsAodMetOfficeTLAD(const ioda::ObsSpace & odb, const Parameters_ &params):
  LinearObsOperatorBase(odb), varin_(), kMatrix_(), trajInit_(false)
{
  // 1. Read in number of dust bins and extinction coefficients from yaml file.
  NDustBins_ = params.NDustBins;  // number of dust bins
  AodKExt_ = params.AodKExt;  // Extinction coefficient for each dust bin

  // check the number of dust bins is supported
  if (NDustBins_ < 1) {
    // raise error as we need to have some dust bins
    throw eckit::UserError("NDustBins_ must be 1 or higher", Here());
  }

  // check the size of AodKExt_ is consistent with NDustBins_
  ASSERT(AodKExt_.size() == NDustBins_);

  // 2. Define which model fields will be needed
  // pressure fields:
  varin_.push_back("air_pressure_levels");
  varin_.push_back("surface_pressure");

  // dust mass concentration:
  std::string dust_var_name;  // dust variable name in geovals
  for (size_t d = 0; d < NDustBins_; d++) {
      dust_var_name = "mass_fraction_of_dust00" + std::to_string(d+1) + "_in_air";
      varin_.push_back(dust_var_name);
  }

  oops::Log::trace() << "ObsAodMetOfficeTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodMetOfficeTLAD::~ObsAodMetOfficeTLAD() {
  trajInit_ = false;
  oops::Log::trace() << "ObsAodMetOfficeTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodMetOfficeTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                        const QCFlags_t & qc_flags) {
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

    // check model fields are ordered from top down, fail if not
    if (plev.front() > plev.back()) {
      throw eckit::BadValue("model fields must be ordered from top down", Here());
    }

    // calculate dAOD/dmass for each bin, level and profile

    // (Re)initialise the Jacobian matrix
    if (!kMatrix_.empty()) kMatrix_.clear();

    kMatrix_.resize(NDustBins_);
    for (size_t d = 0; d < NDustBins_; d++)
    {
        kMatrix_[d].resize(nlevels - 1);
        for (size_t k = 0; k < (nlevels - 1); k++)
        {
           kMatrix_[d][k].resize(nprofiles);
        }
    }

    // loop over profiles
    for (size_t p = 0; p < nprofiles; p++) {
      // loop over dust bins and calculate gradient w.r.t mass concentration per bin:
      for (size_t d = 0; d < NDustBins_; d++) {
        // lowest model layer:
        kMatrix_[d][nlevels - 2][p] = (1.0 / Constants::grav) * (ps[p] - plev[nlevels - 2][p]) *
                AodKExt_[d];
        // loop over the rest of the model layers
        for (size_t k = 0; k < (nlevels - 2); k++) {
          kMatrix_[d][k][p] = (1.0 / Constants::grav) * (plev[k+1][p] - plev[k][p]) *
                  AodKExt_[d];
        }
      }
    }
    trajInit_ = true;
}

void ObsAodMetOfficeTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & hofx,
                                        const QCFlags_t & qc_flags) const {
    // Ensure trajectory has already been calculated
    ASSERT(trajInit_);

    // Get number of obs locations & number of model pressure levels:
    std::size_t nprofiles = geovals.nlocs();
    // number of full (rho) levels
    std::size_t nlevels   = geovals.nlevs(oops::Variable{"air_pressure_levels"});
    // Check hofx size and initialise hofx to zero
    ASSERT(geovals.nlocs() == hofx.nlocs());
    hofx.zero();

    // calculate tangent linear

    std::vector<double> mass_d(nlevels - 1);  // mass concentration increments
    // loop over profiles
    for (size_t p = 0; p < nprofiles; p++) {
      // loop over dust bins and calculate gradient w.r.t mass concentration per bin:
      for (size_t d = 0; d < NDustBins_; d++) {
        oops::Variable dust_var{"mass_fraction_of_dust00" + std::to_string(d+1) + "_in_air"};
        geovals.getAtLocation(mass_d, dust_var, p);
        // loop over model layers
        for (size_t k=0; k < (nlevels - 1); k++) {
          hofx[p] += kMatrix_[d][k][p] * mass_d[k];
        }
      }
    }

  oops::Log::trace() << "ObsAodMetOfficeTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodMetOfficeTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & hofx,
                                        const QCFlags_t & qc_flags) const {
  // Ensure trajectory has already been calculated
  ASSERT(trajInit_);

  // Get number of obs locations & number of model pressure levels:
  std::size_t nprofiles = geovals.nlocs();
  // number of full (rho) levels
  std::size_t nlevels   = geovals.nlevs(oops::Variable{"air_pressure_levels"});

  // Check hofx size
  ASSERT(geovals.nlocs() == hofx.nlocs());

  // Get the missing value indicator
  const double missing = util::missingValue<double>();

  // Loop through the obs, adding the increment to the model state

  std::vector<double> mass_d(nlevels - 1);  // mass concentration increments on half (theta) levels
  // loop over profiles
  for (size_t p = 0; p < nprofiles; p++) {
    if (hofx[p] != missing) {
      // loop over dust bins and calculate gradient w.r.t mass concentration per bin:
      for (size_t d = 0; d < NDustBins_; d++) {
        oops::Variable dust_var{"mass_fraction_of_dust00" + std::to_string(d+1) + "_in_air"};
        geovals.getAtLocation(mass_d, dust_var, p);
        // loop over the model layers
        for (size_t k=0; k < (nlevels - 1); k++) {
          mass_d[k] += kMatrix_[d][k][p] * hofx[p];
        }
        // Store the updated model state increments
        geovals.putAtLocation(mass_d, dust_var, p);
      }
    }
  }

  oops::Log::trace() << "ObsAodMetOfficeTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodMetOfficeTLAD::print(std::ostream & os) const {
  os << "ObsAodMetOfficeTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
