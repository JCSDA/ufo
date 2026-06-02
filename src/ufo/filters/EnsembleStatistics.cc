/* * (C) Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/EnsembleStatistics.h"

#include <algorithm>
#include <set>
#include <utility>

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Buffer.h"
#include "eckit/mpi/Comm.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/mpi/mpi.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Logger.h"

namespace ufo {

namespace {

eckit::mpi::Comm & getOrCreateEnsembleComm(size_t geoCommSize, size_t timeCommSize,
                                           const std::string &name) {
  const size_t ranks = oops::mpi::world().size();
  const size_t ranksPerMember = geoCommSize * timeCommSize;
  const size_t ensembleSize = ranks / ranksPerMember;
  if (ranksPerMember * ensembleSize != ranks)
    throw eckit::UserError("MPI processes must be split evenly across ensemble members");
  const size_t globalRank = oops::mpi::world().rank();
  const size_t localRank = globalRank % ranksPerMember;
  const size_t myMember = globalRank / ranksPerMember;

  const std::string ensembleCommName = name + '_' + std::to_string(localRank);
  // Eckit does not allow creating multiple communicators with the same name (even if they were
  // identical); a SeriousBug exception is thrown on an attempt to do so.
  // So if a communicator with the specified name already exists (perhaps because it's already been
  // created by another instance of the EnsembleStatistics filter, another ObsSpace, or another
  // iteration of the assimilation loop), reuse it instead of creating a new one.
  eckit::mpi::Comm &comm = eckit::mpi::hasComm(ensembleCommName) ?
        eckit::mpi::comm(ensembleCommName) :
        oops::mpi::world().split(localRank, ensembleCommName);
  ASSERT(comm.size() == ensembleSize);
  ASSERT(comm.rank() == myMember);
  return comm;
}

}  // namespace

// -----------------------------------------------------------------------------

EnsembleStatistics::EnsembleStatistics(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                                       std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                       std::shared_ptr<ioda::ObsDataVector<float> > obserr)
: ObsProcessorBase(obsdb, true /*deferToPost?*/, std::move(flags), std::move(obserr)),
  parameters_(parameters)
{
  oops::Log::trace() << "EnsembleStatistics constructor start" << std::endl;
  if (parameters.filterVariables.value().has_value()) {
    const oops::ObsVariables &observedVariables = obsdb.obsvariables();
    const oops::ObsVariables filterVariables = getFilterVariables();
    for (size_t i = 0; i < filterVariables.size(); ++i)
      if (!observedVariables.has(filterVariables[i])) {
        throw eckit::UserError("Filter variable '" + filterVariables[i] +
                               "' is not an observed variable", Here());
      }
  }
  oops::Log::trace() << "EnsembleStatistics constructor complete" << std::endl;
}

// -----------------------------------------------------------------------------

EnsembleStatistics::~EnsembleStatistics() {
  oops::Log::trace() << "EnsembleStatistics destructor" << std::endl;
}

// -----------------------------------------------------------------------------

oops::ObsVariables EnsembleStatistics::getFilterVariables() const {
  if (parameters_.filterVariables.value().has_value()) {
    ufo::Variables vars;
    for (const Variable &var : *parameters_.filterVariables.value())
      vars += var;
    return vars.toOopsObsVariables();
  } else {
    return obsdb_.obsvariables();
  }
}

// -----------------------------------------------------------------------------

void EnsembleStatistics::doFilter() {
  oops::Log::trace() << "EnsembleStatistics doFilter start" << std::endl;

  const std::set<EnsembleStatistic> statistics(parameters_.statistics.value().begin(),
                                               parameters_.statistics.value().end());
  const oops::ObsVariables filterVariables = getFilterVariables();

  if (!statistics.empty() && filterVariables.size() > 0) {
    const eckit::mpi::Comm &ensembleComm = getOrCreateEnsembleComm(
          obsdb_.comm().size(),  obsdb_.commTime().size(), "comm_geo_time");
    const size_t ensembleSize = ensembleComm.size();

    ioda::ObsDataVector<float> hofx(obsdb_, filterVariables, "");
    // By a quirk of the implementation, this will retrieve the values of all variables in `hofx`.
    data_.get(Variable("HofX/" + filterVariables[0]), hofx);

    const size_t nlocs = hofx.nlocs();
    oops::Log::info() << "nlocs = " << nlocs << std::endl;
    const size_t nvars = hofx.nvars();

    // Concatenate the local H(x) vectors of all filter variables so that they can be transferred to
    // other MPI processes in one go.
    std::vector<float> hofxValues;
    hofxValues.reserve(nlocs * nvars);
    for (size_t i = 0; i < nvars; ++i)
      hofxValues.insert(hofxValues.end(), hofx[i].begin(), hofx[i].end());

    std::vector<float> meanHofxValues;
    if (oops::contains(statistics, EnsembleStatistic::HOFX_MEAN) ||
        oops::contains(statistics, EnsembleStatistic::HOFX_STDDEV)) {
      // Calculate the ensemble mean of H(x) vectors. (It is needed also for the calculation of the
      // ensemble spread.)
      meanHofxValues.resize(hofxValues.size());

      // Step 1: Sum over ensemble members.
      ensembleComm.allReduce(hofxValues, meanHofxValues, eckit::mpi::Operation::SUM);

      // Step 2: Divide by the number of members to obtain the mean.
      const double scale = 1.0f / ensembleSize;
      for (float & x : meanHofxValues)
        x *= scale;

      if (oops::contains(statistics, EnsembleStatistic::HOFX_MEAN)) {
        // Save the calculated ensemble means to the ObsSpace.
        std::vector<float> statValues(hofx.nlocs());
        for (size_t i = 0; i < hofx.nvars(); ++i) {
          std::copy(meanHofxValues.begin() + i * nlocs,
                    meanHofxValues.begin() + (i + 1) * nlocs,
                    statValues.begin());
          obsdb_.put_db("HofXEnsembleMean", hofx.varnames()[i], statValues);
        }
      }
    }

    if (oops::contains(statistics, EnsembleStatistic::HOFX_STDDEV)) {
      // Calculate the ensemble spread of H(x) vectors.

      std::vector<float> hofxStdDevValues(hofxValues.size());

      // Step 1: Calculate local squared deviations from the mean.
      for (size_t i = 0; i < hofxValues.size(); ++i) {
        float deviation = hofxValues[i] - meanHofxValues[i];
        hofxStdDevValues[i] = deviation * deviation;
      }

      // Step 2: Sum over ensemble members.
      ensembleComm.allReduceInPlace(hofxStdDevValues.begin(),
                                    hofxStdDevValues.end(), eckit::mpi::Operation::SUM);

      // Step 3: Divide by the number of members and take the square root.
      const float scale = 1.0f / ensembleSize;
      for (float & x : hofxStdDevValues)
        x = std::sqrt(x * scale);

      // Save the calculated ensemble spreads to the ObsSpace.
      std::vector<float> statValues(hofx.nlocs());
      for (size_t i = 0; i < hofx.nvars(); ++i) {
        std::copy(hofxStdDevValues.begin() + i * nlocs,
                  hofxStdDevValues.begin() + (i + 1) * nlocs,
                  statValues.begin());
        obsdb_.put_db("HofXEnsembleStdDev", hofx.varnames()[i], statValues);
      }
    }
  }

  oops::Log::trace() << "EnsembleStatistics doFilter complete" << std::endl;
}

// -----------------------------------------------------------------------------

void EnsembleStatistics::print(std::ostream & os) const {
  os << "Ensemble Statistics: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
