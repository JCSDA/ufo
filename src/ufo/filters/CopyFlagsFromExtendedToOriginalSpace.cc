/*
 * (C) Copyright 2022 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/CopyFlagsFromExtendedToOriginalSpace.h"

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/DiagnosticFlag.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/profile/ProfileAverageUtils.h"

namespace ufo {

// -----------------------------------------------------------------------------

CopyFlagsFromExtendedToOriginalSpace::CopyFlagsFromExtendedToOriginalSpace
                                    (ioda::ObsSpace &obsdb, const Parameters_ &parameters,
                                     std::shared_ptr<ioda::ObsDataVector<int>> flags,
                                     std::shared_ptr<ioda::ObsDataVector<float>> obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "CopyFlagsFromExtendedToOriginalSpace constructor" << std::endl;

  // Ensure observations have been grouped into profiles.
  if (obsdb_.obs_group_vars().empty()) {
    throw eckit::UserError("Group variables configuration is empty", Here());
  }

  // Check the ObsSpace has been extended. If this is not the case
  // then it will not be possible to access profiles in the original and
  // extended sections of the ObsSpace.
  if (!obsdb_.has("MetaData", "extendedObsSpace")) {
    throw eckit::UserError("The extended obs space has not been produced", Here());
  }

  // Get parameters from options
  allvars_ += parameters_.obsVertCoord;
  allvars_ += parameters_.modelVertCoord;
  oops::Log::trace() << "CopyFlagsFromExtendedToOriginalSpace constructed" << std::endl;
}

// -----------------------------------------------------------------------------

CopyFlagsFromExtendedToOriginalSpace::~CopyFlagsFromExtendedToOriginalSpace() {
  oops::Log::trace() << "CopyFlagsFromExtendedToOriginalSpace destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void CopyFlagsFromExtendedToOriginalSpace::applyFilter(const std::vector<bool> & apply,
                                          const Variables & filtervars,
                                          std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "CopyFlagsFromExtendedToOriginalSpace Filter" << std::endl;
  oops::Log::debug() << "CopyFlagsFromExtendedToOriginalSpace obserr: " << *obserr_ << std::endl;

  // Number of locations.
  const size_t nlocs = obsdb_.nlocs();

  // Get model vertical coordinate (in the extended space)
  std::vector<float> modelVertCoord(nlocs);
  data_.get(parameters_.modelVertCoord, modelVertCoord);

  // Get obs vertical coordinate
  std::vector<float> obsVertCoord(nlocs);
  data_.get(parameters_.obsVertCoord, obsVertCoord);

  // Correspondence between record numbers and indices in the data sample.
  const std::vector<size_t> &unique_recnums = obsdb_.recidx_all_recnums();
  oops::Log::debug() << "Unique record numbers" << std::endl;

  // Number of profiles in the original ObsSpace.
  const size_t nprofs = unique_recnums.size() / 2;

  // number of filter variables
  const size_t nv = filtervars.size();
  // loop over filter variables
  for (size_t jv = 0; jv < nv; ++jv) {
    const std::string varName = filtervars.variable(jv).variable();
    const std::string flagName = filtervars.variable(jv).group();
    oops::Log::debug() << flagName << "/" << varName << std::endl;
    if (flagName.find("DiagnosticFlags") == std::string::npos) {
        throw eckit::UserError("Filter variable must be a Diagnostic Flag.", Here());
    }
    std::vector<DiagnosticFlag> varToCopy;
    obsdb_.get_db(flagName, varName, varToCopy);

    for (std::size_t jprof = 0; jprof < nprofs; ++jprof) {
        oops::Log::debug() << "Profile " << (jprof + 1) << " / " << nprofs << std::endl;
        // Get locations of profile in the original ObsSpace and
        // the corresponding profile in the extended ObsSpace.
        // Assuming the extended ObsSpace has been configured correctly, which is
        // checked above, the profile in the extended ObsSpace is always located
        // nprofs positions further on than the profile in the original ObsSpace.
        const std::vector<size_t> &locsOriginal = obsdb_.recidx_vector(unique_recnums[jprof]);
        const std::vector<size_t> &locsExt = obsdb_.recidx_vector(unique_recnums[jprof + nprofs]);
        oops::Log::debug() << "locsOriginal: " << locsOriginal << std::endl;
        oops::Log::debug() << "locsExt: " << locsExt << std::endl;

        // Create a vector mapping each obs index to a model index (many-to-1):
        const std::vector<size_t> linkObsAndModInds = ufo::linkObsAndModLevIndices
                                                           (locsOriginal,
                                                            locsExt,
                                                            obsVertCoord,
                                                            modelVertCoord,
                                                            nlocs,
                                                            true);  // link obs->mod only
        // compute
        for (size_t jlev : locsOriginal) {
          if (apply[jlev] && !varToCopy[jlev]) {  // don't unset flags that are already set
            varToCopy[jlev] = varToCopy[linkObsAndModInds[jlev]];
          }
        }  // obs location jlev
    }  // profile jprof
    obsdb_.put_db(flagName, varName, varToCopy);
  }
  oops::Log::trace() << "CopyFlagsFromExtendedToOriginalSpace::applyFilter finished" << std::endl;
}

// -----------------------------------------------------------------------------

void CopyFlagsFromExtendedToOriginalSpace::print(std::ostream & os) const {
  os << "CopyFlagsFromExtendedToOriginalSpace: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
