/*
 * (C) Copyright 2021- UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/lightning/ObsLightning.h"

#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ioda/ObsVector.h"
#include "oops/base/Locations.h"
#include "oops/base/Variables.h"
#include "oops/interface/SampledLocations.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsTraits.h"
#include "ufo/SampledLocations.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsLightning> makerLightning_("Lightning");
// -----------------------------------------------------------------------------

ObsLightning::ObsLightning(const ioda::ObsSpace & odb, const Parameters_ & params)
  : ObsOperatorBase(odb), keyOper_(0), odb_(odb), varin_(), nhoriz_(0)
{
  // Set local variables from YAML params.  The operator wants the square of number of gridpoints.
  nhoriz_ = params.num_gridpoints.value()*params.num_gridpoints.value();
  // l_fed_nonlinear_ = params.use_nonlinear.value(); //Keep sample logical paramter here

  ufo_lightning_setup_f90(keyOper_, nhoriz_, odb.obsvariables(), varin_);

  oops::Log::trace() << "ObsLightning created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsLightning::~ObsLightning() {
  ufo_lightning_delete_f90(keyOper_);
  oops::Log::trace() << "ObsLightning destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsLightning::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                               ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  ufo_lightning_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                           ovec.toFortran());
  oops::Log::trace() << "ObsLightning: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

// This method generates extended geographic locations (longitude and latitude)
// for each observation by calling a Fortran subroutine to calculate the 2D locations
// for a predefined area around each observation point.
  ObsLightning::Locations_ ObsLightning::locations() const {
  typedef oops::SampledLocations<ObsTraits> SampledLocations_;

  // Initialize vectors to hold extended longitudes and latitudes for each location
  std::vector<float> lons(odb_.nlocs()*nhoriz_);
  std::vector<float> lats(odb_.nlocs()*nhoriz_);
  std::vector<util::DateTime> times(odb_.nlocs()*nhoriz_);

  std::vector<util::DateTime> times_notduplicated(odb_.nlocs());
  odb_.get_db("MetaData", "dateTime", times_notduplicated);
  for (size_t jloc = 0; jloc < odb_.nlocs(); ++jloc) {
    for (size_t jhoriz = 0; jhoriz < nhoriz_; ++jhoriz) {
      times[jloc*nhoriz_ + jhoriz] = times_notduplicated[jloc];
    }
  }

  // Call to Fortran subroutine to fill 'lons' and 'lats' with calculated geographic locations
  ufo_lightning_extend_geovals_f90(keyOper_, odb_, lons.size(), lons[0], lats[0]);

  // Print dimention of lons
  // std::cout << "Lons size: " << lons.size() << std::endl;
  //
  // Find and print minimum vlaue of lons
  // auto min_lon = *std::min_element(lons.begin(), lons.end());
  // std::cout << "Min Lons value: " << min_lon << std::endl;
  //
  // Find and print maximum value of lons
  // auto max_lon = *std::max_element(lons.begin(), lons.end());
  // std::cout << "Max Lons value: " << max_lon << std::endl;
  //
  // Group the paths by location to associate each set of extended locations with
  // the original observation
  std::vector<util::Range<size_t>> pathsGroupedByLocation(odb_.nlocs());
  for (size_t jloc = 0; jloc < odb_.nlocs(); ++jloc) {
    pathsGroupedByLocation[jloc].begin = jloc * nhoriz_;
    pathsGroupedByLocation[jloc].end = (jloc + 1) * nhoriz_;
  }

  // Construct and return a SampledLocations object populated with the extended locations
  return SampledLocations_(
    std::make_unique<SampledLocations>(lons, lats, times, odb_.distribution(),
                                       std::move(pathsGroupedByLocation)));
}

// -----------------------------------------------------------------------------

void ObsLightning::computeReducedVars(const oops::Variables & reducedVars,
                                      GeoVaLs & /*geovals*/) const {
  // No method for reducing the set of nhoriz_ profiles sampling each location into a single
  // profile has been implemented so far, so when this obs operator is in use, neither it nor
  // any obs filters or bias predictors can request variables in the reduced format.
  // This procedure is called to bypass computeReducedVars
  if (reducedVars.size() != 0)
    throw eckit::NotImplemented("ObsLightning is unable to compute the reduced "
                                "representation of GeoVaLs", Here());
}

// -----------------------------------------------------------------------------

void ObsLightning::print(std::ostream & os) const {
  os << "ObsLightning::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
