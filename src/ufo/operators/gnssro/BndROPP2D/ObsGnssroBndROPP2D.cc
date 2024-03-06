/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/gnssro/BndROPP2D/ObsGnssroBndROPP2D.h"

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
static ObsOperatorMaker<ObsGnssroBndROPP2D> makerGnssroBndROPP2D_("GnssroBndROPP2D");
// -----------------------------------------------------------------------------

ObsGnssroBndROPP2D::ObsGnssroBndROPP2D(const ioda::ObsSpace & odb,
                                       const Parameters_ & params)
  : ObsOperatorBase(odb), keyOperGnssroBndROPP2D_(0), odb_(odb), varin_(),
    nhoriz_(params.options.value().nHoriz)
{
  const std::vector<std::string> vv{"air_temperature", "specific_humidity", "air_pressure",
                                    "geopotential_height", "surface_altitude"};
  varin_.reset(new oops::Variables(vv));

  ufo_gnssro_bndropp2d_setup_f90(keyOperGnssroBndROPP2D_,
                                 params.options.value().toConfiguration(), odb_.nlocs());
  oops::Log::trace() << "ObsGnssroBndROPP2D created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBndROPP2D::~ObsGnssroBndROPP2D() {
  ufo_gnssro_bndropp2d_delete_f90(keyOperGnssroBndROPP2D_);
  oops::Log::trace() << "ObsGnssroBndROPP2D destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2D::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                     ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  ufo_gnssro_bndropp2d_simobs_f90(keyOperGnssroBndROPP2D_, gom.toFortran(), odb_,
                                  ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
ObsGnssroBndROPP2D::Locations_
ObsGnssroBndROPP2D::locations() const {
  typedef oops::SampledLocations<ObsTraits> SampledLocations_;

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
  ufo_gnssro_2d_locs_init_f90(keyOperGnssroBndROPP2D_, odb_, lons.size(), lons[0], lats[0]);

  std::vector<util::Range<size_t>> pathsGroupedByLocation(odb_.nlocs());
  for (size_t jloc = 0; jloc < odb_.nlocs(); ++jloc) {
    pathsGroupedByLocation[jloc].begin = jloc * nhoriz_;
    pathsGroupedByLocation[jloc].end = (jloc + 1) * nhoriz_;
  }

  return SampledLocations_(
          std::make_unique<SampledLocations>(lons, lats, times, odb_.distribution(),
                                             std::move(pathsGroupedByLocation)));
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2D::computeReducedVars(const oops::Variables & reducedVars,
                                            GeoVaLs & /*geovals*/) const {
  // No method for reducing the set of nhoriz_ profiles sampling each location into a single profile
  // has been implemented so far, so when this obs operator is in use, neither it nor any obs
  // filters or bias predictors can request variables in the reduced format.
  if (reducedVars.size() != 0)
    throw eckit::NotImplemented("ObsGnssroBndROPP2D is unable to compute the reduced "
                                "representation of GeoVaLs", Here());
}

// -----------------------------------------------------------------------------

void ObsGnssroBndROPP2D::print(std::ostream & os) const {
  os << "ObsGnssroBndROPP2D::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
