/*
 * (C) Copyright 2025-2026 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/crtm/ObsExtCoeffProfCRTM.h"

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Locations.h"
#include "oops/interface/SampledLocations.h"
#include "oops/util/dateFunctions.h"
#include "oops/util/Range.h"
#include "oops/util/TimeWindow.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsTraits.h"
#include "ufo/operators/crtm/ObsExtCoeffProfCRTM.interface.h"
#include "ufo/SampledLocations.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsExtCoeffProfCRTM>
       makerExtCoeffProfCRTM_("ExtinctionCoefficientProfileCRTM");

// -----------------------------------------------------------------------------

ObsExtCoeffProfCRTM::ObsExtCoeffProfCRTM(const ioda::ObsSpace & odb,
                       const Parameters_ & parameters)
  : ObsOperatorBase(odb), keyOperExtCoeffProfCRTM_(0), odb_(odb), varin_(),
    parameters_(parameters)
{
  // parse channels from the config and create variable names
  const oops::ObsVariables & observed = odb.assimvariables();
  std::vector<int> channels_list = observed.channels();

  // get a single central of middle time from observation space
  const util::DateTime midPoint = odb.timeWindow().midpoint();
  std::string year, month, day, hour, minute, second;
  midPoint.toYYYYMMDDhhmmss(year, month, day, hour,  minute,  second);
  // Julian Day Number since noon Universal Time (UT) on January 1, 4713 BCE
  uint64_t midPointJulday = util::datefunctions::dateToJulian(std::stoi(year),
                                                              std::stoi(month),
                                                              std::stoi(day));

  // call Fortran setup routine
  ufo_extcoeffprofcrtm_setup_f90(keyOperExtCoeffProfCRTM_, parameters_.toConfiguration(),
                        channels_list.size(), channels_list[0], midPointJulday,
                        varin_, odb.comm());
  oops::Log::trace() << "ObsExtCoeffProfCRTM constructor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsExtCoeffProfCRTM::~ObsExtCoeffProfCRTM() {
  ufo_extcoeffprofcrtm_delete_f90(keyOperExtCoeffProfCRTM_);
  oops::Log::trace() << "ObsExtCoeffProfCRTM destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExtCoeffProfCRTM::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                             ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  // Max layer count across all profiles (records).
  int nobsLayer = 0;
  for (auto it = odb_.recidx_begin(); it != odb_.recidx_end(); ++it) {
    nobsLayer = std::max(nobsLayer, static_cast<int>(odb_.recidx_vector(it).size()));
  }
  ufo_extcoeffprofcrtm_simobs_f90(keyOperExtCoeffProfCRTM_, gom.toFortran(), odb_,
                          ovec.nvars(), ovec.nlocs(), nobsLayer, ovec.toFortran());
}

// -----------------------------------------------------------------------------

// Locations are one-per-(profile,layer), but every layer of a profile
// shares the same lat/lon/time [true for CALIOP], so the model column GetValues
// interpolates for one of them is identical for all the others.
// Requesting only one interpolation path per record (instead of
// one path per location) reduces recomputing that same column once per layer per profile.
ObsExtCoeffProfCRTM::Locations_ ObsExtCoeffProfCRTM::locations() const {
  typedef oops::SampledLocations<ObsTraits> SampledLocations_;

  if (odb_.obs_group_vars().empty()) {
    throw eckit::UserError("ObsExtCoeffProfCRTM::locations(): obs space is not grouped into "
                           "records -- add 'obsgrouping: group variables: [sequenceNumber]' "
                           "to the obs space YAML so flattened (profile,layer) locations are "
                           "grouped into physical profiles.", Here());
  }

  const size_t nlocs = odb_.nlocs();
  const size_t n_Profiles = odb_.nrecs();

  std::vector<float> lons(nlocs), lats(nlocs);
  std::vector<util::DateTime> times(nlocs);
  odb_.get_db("MetaData", "longitude", lons);
  odb_.get_db("MetaData", "latitude", lats);
  odb_.get_db("MetaData", "dateTime", times);

  std::vector<float> pathLons(n_Profiles), pathLats(n_Profiles);
  std::vector<util::DateTime> pathTimes(n_Profiles);
  std::vector<util::Range<size_t>> pathsGroupedByLocation(nlocs);

  size_t profileIdx = 0;
  for (auto it = odb_.recidx_begin(); it != odb_.recidx_end(); ++it) {
    const std::vector<size_t> & locsInRecord = odb_.recidx_vector(it);
    ASSERT_MSG(!locsInRecord.empty(),
              "ObsExtCoeffProfCRTM::locations(): ioda record has no locations");
    const size_t first_loc = locsInRecord[0];
    pathLons[profileIdx] = lons[first_loc];
    pathLats[profileIdx] = lats[first_loc];
    pathTimes[profileIdx] = times[first_loc];
    for (const size_t loc : locsInRecord) {
      pathsGroupedByLocation[loc] = util::Range<size_t>{profileIdx, profileIdx + 1};
    }
    ++profileIdx;
  }
  ASSERT(profileIdx == n_Profiles);

  return SampledLocations_(
      std::make_unique<SampledLocations>(pathLons, pathLats, pathTimes, odb_.distribution(),
                                         std::move(pathsGroupedByLocation)));
}

// -----------------------------------------------------------------------------

void ObsExtCoeffProfCRTM::computeReducedVars(const oops::Variables &reducedVars,
                                             GeoVaLs & geovals) const {
  // Not providing a method for reducing the variable to one profile per location.
  // so when this obs operator is in use, neither it nor any obs
  // filters or bias predictors can request variables in the reduced format.
  if (reducedVars.size() != 0) {
    throw eckit::NotImplemented("ObsExtCoeffProfCRTM is unable to compute the reduced "
                                "representation of GeoVaLs", Here());
  }
}

// -----------------------------------------------------------------------------

void ObsExtCoeffProfCRTM::print(std::ostream & os) const {
  os << "ObsExtCoeffProfCRTM::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
