/*
 * (C) Crown copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ioda/ObsDataVector.h"

#include "oops/mpi/mpi.h"

#include "ufo/filters/obsfunctions/MPIRank.h"

namespace ufo {

static ObsFunctionMaker<MPIRank> makerMPIRank_("MPIRank");

// -----------------------------------------------------------------------------

MPIRank::MPIRank(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::trace() << "MPIRank constructor" << std::endl;
  // Validate and deserialize options
  options_.validateAndDeserialize(conf);
}

// -----------------------------------------------------------------------------

MPIRank::~MPIRank() {
  oops::Log::trace() << "MPIRank destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void MPIRank::compute(const ObsFilterData & in,
                      ioda::ObsDataVector<int> & out) const {
  oops::Log::trace() << "MPIRank compute start" << std::endl;

  const ioda::ObsSpace & obsdb = in.obsspace();
  const std::size_t nlocs = obsdb.nlocs();
  const std::size_t gnlocs = obsdb.globalNumLocs();
  const int rank = static_cast<int>(obsdb.comm().rank());

  // Vector signifiying whether each location on this rank is a 'patch' observation.
  std::vector<bool> patchObsVec(nlocs);
  obsdb.distribution()->patchObs(patchObsVec);

  // Fill vector of ranks of patch observations.
  std::vector<int> patchRanks;
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    if (patchObsVec[jloc]) {
      patchRanks.push_back(rank);
    }
  }

  // Gather ranks of patch observations.
  // The use of oops::mpi::allGatherv() ensures the gather procedure
  // works for subclasses of ioda::Distribution whose implementation of allGatherv() may
  // not produce the desired result.
  std::vector<int> globalPatchRanks(patchRanks);
  oops::mpi::allGatherv(obsdb.comm(), globalPatchRanks);
  ASSERT(globalPatchRanks.size() == gnlocs);

  // Indices of locations on this rank in the global location vector.
  std::vector<std::size_t> index(nlocs);
  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    index[jloc] =
      obsdb.distribution()->globalUniqueConsecutiveLocationIndex(jloc);
  }

  // Write out the ranks of the patch observations that correspond to locations on this rank.
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    out[0][jloc] = globalPatchRanks[index[jloc]];
  }

  oops::Log::trace() << "MPIRank compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & MPIRank::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
