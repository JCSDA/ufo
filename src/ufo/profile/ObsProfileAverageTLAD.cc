/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ObsProfileAverageTLAD.h"

#include <iomanip>
#include <ostream>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsProfileAverageTLAD> obsProfileAverageMaker_("ProfileAverage");
// -----------------------------------------------------------------------------

ObsProfileAverageTLAD::ObsProfileAverageTLAD(const ioda::ObsSpace & odb,
                                             const Parameters_ & parameters)
  : LinearObsOperatorBase(odb, VariableNameMap(parameters.AliasFile.value())),
    odb_(odb), data_(odb, parameters, nameMap_)
{
  requiredVars_ += data_.requiredVars();
  oops::Log::trace() << "ObsProfileAverageTLAD constructed" << std::endl;
}

// -----------------------------------------------------------------------------

ObsProfileAverageTLAD::~ObsProfileAverageTLAD() {
  oops::Log::trace() << "ObsProfileAverageTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProfileAverageTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                          const QCFlags_t & qc_flags) {
  // Cache the model trajectory for use in the slant path location algorithm.
  data_.cacheGeoVaLs(geovals);
  oops::Log::trace() << "ObsProfileAverageTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProfileAverageTLAD::simulateObsTL(const GeoVaLs & dx, ioda::ObsVector & dy,
                                          const QCFlags_t &  qc_flags) const {
  oops::Log::trace() << "ObsProfileAverageTLAD: simulateObsTL started" << std::endl;

  // Get correspondence between record numbers and indices in the total sample.
  const std::vector<std::size_t> &recnums = odb_.recidx_all_recnums();

  // Number of profiles in the original ObsSpace.
  const std::size_t nprofs = recnums.size() / 2;

  // Loop over profiles.
  for (std::size_t jprof = 0; jprof < nprofs; ++jprof) {
    const std::vector<std::size_t> &locsOriginal = odb_.recidx_vector(recnums[jprof]);
    const std::vector<std::size_t> &locsExtended = odb_.recidx_vector(recnums[jprof + nprofs]);

    // Retrieve slant path locations.
    const std::vector<std::size_t>& slant_path_location =
      data_.getSlantPathLocations(locsOriginal, locsExtended);

    for (int jvar : data_.operatorVarIndices()) {
      const auto& variable = nameMap_.convertName(dy.varnames().variables()[jvar]);
      const std::size_t nlevs_var = dx.nlevs(variable);
      std::vector<double> var_gv(nlevs_var);
      for (std::size_t mlev = 0; mlev < nlevs_var; ++mlev) {
        const std::size_t jloc = slant_path_location[mlev];
        dx.getAtLocation(var_gv, variable, jloc);
        if (data_.geovalsObsSameDir()) {  // geovals and observations are the same way round:
          dy[locsExtended[mlev] * dy.nvars() + jvar] = var_gv[mlev];
        } else {  // reverse geovals so they're the same way round in extended space as
          // observations / H(x) in original space:
          dy[locsExtended[mlev] * dy.nvars() + jvar] = var_gv[nlevs_var - 1 - mlev];
        }
      }
    }
  }

  oops::Log::trace() << "ObsProfileAverageTLAD: simulateObsTL finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsProfileAverageTLAD::simulateObsAD(GeoVaLs & dx, const ioda::ObsVector & dy,
                                          const QCFlags_t &  qc_flags) const {
  oops::Log::trace() << "ObsProfileAverageTLAD: simulateObsAD started" << std::endl;

  const double missing = util::missingValue<double>();

  // Get correspondence between record numbers and indices in the total sample.
  const std::vector<std::size_t> &recnums = odb_.recidx_all_recnums();

  // Number of profiles in the original ObsSpace.
  const std::size_t nprofs = recnums.size() / 2;

  // Loop over profiles.
  for (std::size_t jprof = 0; jprof < nprofs; ++jprof) {
    const std::vector<std::size_t> &locsOriginal = odb_.recidx_vector(recnums[jprof]);
    const std::vector<std::size_t> &locsExtended = odb_.recidx_vector(recnums[jprof + nprofs]);

    // Retrieve slant path locations.
    const std::vector<std::size_t>& slant_path_location =
      data_.getSlantPathLocations(locsOriginal, locsExtended);

    for (int jvar : data_.operatorVarIndices()) {
      const auto& variable = nameMap_.convertName(dy.varnames().variables()[jvar]);
      const std::size_t nlevs_var = dx.nlevs(variable);
      std::vector<double> var_gv(nlevs_var);
      for (std::size_t mlev = 0; mlev < nlevs_var; ++mlev) {
        const std::size_t jloc = slant_path_location[mlev];
        // Get the current value of dx.
        dx.getAtLocation(var_gv, variable, jloc);
        const std::size_t idx = locsExtended[mlev] * dy.nvars() + jvar;
        if (dy[idx] != missing) {
          if (data_.geovalsObsSameDir()) {
            // geovals and observations are the same way round:
            var_gv[mlev] += dy[idx];
          } else {
            // geovals reversed relative to obs
            // write dy into correct index of geoval
            var_gv[nlevs_var - 1 - mlev] += dy[idx];
          }
        }
        // Store the new value of dx.
        dx.putAtLocation(var_gv, variable, jloc);
      }
    }
  }

  oops::Log::trace() << "ObsProfileAverageTLAD: simulateObsAD finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsProfileAverageTLAD::print(std::ostream & os) const {
  os << "ObsProfileAverageTLAD operator" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
