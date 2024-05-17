/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ObsProfileAverage.h"

#include <algorithm>
#include <ostream>
#include <utility>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsProfileAverage> obsProfileAverageMaker_("ProfileAverage");
// -----------------------------------------------------------------------------

ObsProfileAverage::ObsProfileAverage(const ioda::ObsSpace & odb,
                                     const Parameters_ & parameters)
  : ObsOperatorBase(odb, VariableNameMap(parameters.AliasFile.value())),
    odb_(odb), data_(odb, parameters, nameMap_)
{
  requiredVars_ += data_.requiredVars();
  oops::Log::trace() << "ObsProfileAverage constructed" << std::endl;
}

// -----------------------------------------------------------------------------

ObsProfileAverage::~ObsProfileAverage() {
  oops::Log::trace() << "ObsProfileAverage destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProfileAverage::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                    ObsDiagnostics & diags,
                                    const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsProfileAverage: simulateObs started" << std::endl;

  // Cache the input GeoVaLs for use in the slant path location algorithm.
  data_.cacheGeoVaLs(gv);

  // Get correspondence between record numbers and indices in the total sample.
  const std::vector<std::size_t> &recnums = odb_.recidx_all_recnums();

  // Number of profiles in the original ObsSpace.
  const std::size_t nprofs = recnums.size() / 2;

  // Loop over profiles.
  for (std::size_t jprof = 0; jprof < nprofs; ++jprof) {
    oops::Log::debug() << "Profile " << (jprof + 1) << " / " << nprofs << std::endl;

    // Get locations of profile in the original ObsSpace and
    // the corresponding profile in the extended ObsSpace.
    // Assuming the extended ObsSpace has been configured correctly, which is
    // checked in the constructor of the data_ member variable,
    // the profile in the extended ObsSpace is always located
    // nprofs positions further on than the profile in the original ObsSpace.
    const std::vector<std::size_t> &locsOriginal = odb_.recidx_vector(recnums[jprof]);
    const std::vector<std::size_t> &locsExtended = odb_.recidx_vector(recnums[jprof + nprofs]);

    // Retrieve slant path locations.
    const std::vector<std::size_t>& slant_path_location =
      data_.getSlantPathLocations(locsOriginal, locsExtended);

    // Fill H(x) vector for each variable.
    for (int jvar : data_.operatorVarIndices()) {
      const auto& variable = nameMap_.convertName(ovec.varnames().variables()[jvar]);
      // Number of levels for this variable.
      const std::size_t nlevs_var = gv.nlevs(variable);
      // GeoVaL vector for this variable.
      std::vector<double> var_gv(nlevs_var);
      // For each level:
      // - get the relevant slant path location,
      // - retrieve the GeoVaL at that location,
      // - fill H(x) with the relevant level in the GeoVaL.
      for (std::size_t mlev = 0; mlev < nlevs_var; ++mlev) {
        const std::size_t jloc = slant_path_location[mlev];
        gv.getAtLocation(var_gv, variable, jloc);
        if (data_.geovalsObsSameDir()) {  // geovals and observations are the same way round:
          ovec[locsExtended[mlev] * ovec.nvars() + jvar] = var_gv[mlev];
        } else {  // reverse geovals so they're the same way round in extended space as
          // observations / H(x) in original space:
          ovec[locsExtended[mlev] * ovec.nvars() + jvar] = var_gv[nlevs_var - 1 - mlev];
        }
      }
    }
  }

  oops::Log::trace() << "ObsProfileAverage: simulateObs finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsProfileAverage::print(std::ostream & os) const {
  data_.print(os);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
