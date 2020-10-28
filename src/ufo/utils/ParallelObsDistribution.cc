/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/ParallelObsDistribution.h"
#include "oops/mpi/mpi.h"

namespace ufo {

ParallelObsDistribution::ParallelObsDistribution(const ioda::ObsSpace &obsspace)
  : globalObsCount_(obsspace.gnlocs())
{
  const size_t numProcs = obsspace.comm().size();

  localObsCounts_.resize(numProcs);
  obsspace.comm().allGather(static_cast<int>(obsspace.nlocs()),
                            localObsCounts_.begin(), localObsCounts_.end());

  localObsIdDisplacements_.resize(numProcs);
  localObsIdDisplacements_[0] = 0;
  for (size_t i = 1; i < numProcs; ++i)
    localObsIdDisplacements_[i] = localObsIdDisplacements_[i - 1] + localObsCounts_[i - 1];
}

// Generic implementation
template <typename T>
std::vector<T> getGlobalVariableValues(const ioda::ObsSpace &obsspace,
                                       const ParallelObsDistribution &obsDistribution,
                                       const std::string &variable,
                                       const std::string &group) {
  std::vector<T> localValues(obsspace.nlocs());
  obsspace.get_db(group, variable, localValues);

  std::vector<T> globalValues(obsDistribution.globalObsCount());
  obsspace.comm().allGatherv(localValues.begin(), localValues.end(),
                             globalValues.begin(),
                             obsDistribution.localObsCounts().data(),
                             obsDistribution.localObsIdDisplacements().data());
  return globalValues;
}

// Specialisation for date/time variables.
template <>
std::vector<util::DateTime> getGlobalVariableValues(const ioda::ObsSpace &obsspace,
                                                    const ParallelObsDistribution &obsDistribution,
                                                    const std::string &variable,
                                                    const std::string &group) {
  std::vector<util::DateTime> localValues(obsspace.nlocs());
  obsspace.get_db(group, variable, localValues);

  std::vector<util::DateTime> globalValues(obsDistribution.globalObsCount());
  oops::mpi::allGathervUsingSerialize(obsspace.comm(),
                                      localValues.begin(), localValues.end(), globalValues.begin());
  return globalValues;
}

// Explicit instantiations for the variable types supported by ioda::ObsSpace.
#define INSTANTIATE_GET_GLOBAL_VARIABLE_VALUES(TYPE) \
  template std::vector<TYPE> getGlobalVariableValues(const ioda::ObsSpace &, \
                                                     const ParallelObsDistribution &, \
                                                     const std::string &, \
                                                     const std::string &)
INSTANTIATE_GET_GLOBAL_VARIABLE_VALUES(int);
INSTANTIATE_GET_GLOBAL_VARIABLE_VALUES(float);
INSTANTIATE_GET_GLOBAL_VARIABLE_VALUES(double);
// It's unnecessary to instantiate the template function for T = util::DateTime, since a
// specialization exists for that value of T.

}  // namespace ufo
