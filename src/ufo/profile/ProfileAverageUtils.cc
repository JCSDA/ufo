/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <set>
#include <utility>

#include "ufo/filters/QCflags.h"

#include "ufo/profile/ProfileAverageUtils.h"
#include "ufo/profile/ProfileDataHolder.h"

namespace ufo {

  void ProfileAverageUtils::passNonMissingAveragedObservations
  (ProfileDataHandler & profileDataHandler,
   std::vector <ProfileDataHolder> & profiles,
   const std::string & flag_name,
   const std::string & obs_name)
  {
    const float missing = util::missingValue(missing);
    const size_t halfnprofs = profileDataHandler.getObsdb().nrecs() / 2;

    // Retrieve filter QC flags.
    ioda::ObsDataVector<int> & filterFlagsAll = profileDataHandler.getFilterFlags();
    ioda::ObsDataRow<int> & filterFlags = filterFlagsAll[flag_name];

    // The averaging routines fill in previously missing values, so the associated flag is updated.
    const std::vector<size_t> &recnums = profileDataHandler.getObsdb().recidx_all_recnums();
    for (std::size_t jprof = halfnprofs; jprof < 2 * halfnprofs; ++jprof) {
      const std::vector<std::size_t> &locsExtended =
        profileDataHandler.getObsdb().recidx_vector(recnums[jprof]);
      auto& profileExtended = profiles[jprof];
      const std::vector <float> &obs = profileExtended.get<float>(obs_name);
      for (size_t idx = 0; idx < locsExtended.size(); ++idx) {
        if (obs[idx] != missing)
          filterFlags[locsExtended[idx]] = QCflags::pass;
      }
    }
  }

  void ProfileAverageUtils::fillValidationData
  (ProfileDataHolder & profile,
   bool extended_obs_space,
   const std::string & average_name,
   const std::string & qcflags_name,
   const std::string & geovals_testreference_name,
   const std::string & geovals_qcflags_name)
  {
    const float missing = util::missingValue(missing);

    // Retrieve, then save, the OPS versions of the variable averaged onto model levels,
    // and the QC flags associated with the averaging process.
    // The quantities retrieved depend on whether the profile lies in the original
    // or averaged sections of the ObsSpace.
    if (extended_obs_space) {
      std::vector<float>& averaged_values =
        profile.getGeoVaLVector(geovals_testreference_name);
      // Ensure all vectors are the correct size to be saved to the ObsSpace.
      const size_t numModelLevels = profile.getNumProfileLevels();
      averaged_values.resize(numModelLevels, missing);
      profile.set<float>("OPS_" + std::string(average_name), std::move(averaged_values));
      // The QC flags are stored as floats but are converted to integers here.
      // Due to the loss of precision, 5 must be added to the missing value.
      const std::vector <float>& average_qcflags_float =
        profile.getGeoVaLVector(geovals_qcflags_name);
      std::vector <int> average_qcflags_int
        (average_qcflags_float.begin(),
         average_qcflags_float.end());
      std::replace(average_qcflags_int.begin(),
                   average_qcflags_int.end(),
                   -2147483648L,
                   -2147483643L);
      // Ensure all vectors are the correct size to be saved to the ObsSpace.
      average_qcflags_int.resize(numModelLevels, 0);
      profile.set<int>("OPS_" + std::string(qcflags_name), std::move(average_qcflags_int));
    } else {
      // Create a copy here because the vector will be used later in the routine.
      std::vector <float> avg = profile.get<float>(average_name);
      profile.set<float>("OPS_" + std::string(average_name), std::move(avg));
    }
  }
}  // namespace ufo
