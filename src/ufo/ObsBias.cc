/*
 * (C) Copyright 2017-2024 UCAR
 * (C) Crown Copyright 2024, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsBias.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/distribution/Accumulator.h"
#include "ioda/distribution/Distribution.h"
#include "ioda/Engines/EngineUtils.h"
#include "ioda/Engines/HH.h"
#include "ioda/Layout.h"
#include "ioda/ObsGroup.h"
#include "ioda/ObsSpace.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ioda/ObsVector.h"

#include "ufo/ObsBiasIncrement.h"
#include "ufo/utils/IodaGroupIndices.h"
#include "ufo/utils/SaveBiasCoeffs.h"

namespace ufo {

namespace {

/// \brief Mode of the histogram of uncorrected departures `y - H(x)`, per (record, variable).
///
/// The numerical core of the VarBC cold start, following step [1] of WRFDA's
/// `da_varbc_coldstart`. Only entries flagged in \p wanted are accumulated.
std::vector<double> departureModes(const ioda::ObsSpace & odb, const ioda::ObsVector & hofx,
                                   const ioda::ObsDataVector<int> & qcFlags,
                                   const std::vector<bool> & wanted, const size_t nrecs,
                                   const size_t nbins, const double halfWidth,
                                   const size_t minObs, const eckit::mpi::Comm & commTime) {
  const std::vector<std::string> varnames = odb.assimvariables().variables();
  const size_t nvars = varnames.size();
  const size_t nlocs = odb.nlocs();
  const double missing = util::missingValue<double>();
  ASSERT(nbins > 0 && halfWidth > 0.0);
  const double binWidth = 2.0 * halfWidth / static_cast<double>(nbins);

  // Histogram only the wanted entries, so a hyperspectral sounder with a few new channels does
  // not allocate a histogram per channel.
  std::vector<int> slot(wanted.size(), -1);
  size_t ncold = 0;
  for (size_t je = 0; je < wanted.size(); ++je)
    if (wanted[je]) slot[je] = static_cast<int>(ncold++);

  std::vector<double> modes(wanted.size(), missing);
  if (ncold == 0) return modes;

  // The Accumulator avoids double-counting observations duplicated across MPI tasks.
  auto hist = odb.distribution()->createAccumulator<std::size_t>(ncold * nbins);
  auto count = odb.distribution()->createAccumulator<std::size_t>(ncold);
  const std::vector<std::size_t> recnums = odb.recidx_all_recnums();

  std::vector<double> obsValues(nlocs);
  for (size_t jvar = 0; jvar < nvars; ++jvar) {
    odb.get_db("ObsValue", varnames[jvar], obsValues);
    for (size_t jloc = 0; jloc < nlocs; ++jloc) {
      if (qcFlags[jvar][jloc] != 0) continue;
      const double obs = obsValues[jloc];
      const double sim = hofx[jloc * nvars + jvar];
      if (obs == missing || sim == missing) continue;

      size_t jrec = 0;
      if (nrecs > 1)
        jrec = std::find(recnums.begin(), recnums.end(), odb.recnum()[jloc]) - recnums.begin();
      const int icold = slot[jrec * nvars + jvar];
      if (icold < 0) continue;

      const std::int64_t ibin = std::llround((obs - sim + halfWidth) / binWidth);
      if (ibin <= 0 || ibin > static_cast<std::int64_t>(nbins)) continue;
      hist->addTerm(jloc, icold * nbins + (ibin - 1), 1);
      count->addTerm(jloc, icold, 1);
    }
  }

  std::vector<std::size_t> h = hist->computeResult();
  std::vector<std::size_t> n = count->computeResult();
  commTime.allReduceInPlace(h.begin(), h.end(), eckit::mpi::sum());   // 4D sub-windows
  commTime.allReduceInPlace(n.begin(), n.end(), eckit::mpi::sum());

  for (size_t je = 0; je < wanted.size(); ++je) {
    const int icold = slot[je];
    if (icold < 0) continue;
    if (n[icold] < minObs) {
      // Leave the coefficients alone rather than wiping them: under `force` that keeps the
      // previous cycle's values, which beats falling back to zero.
      oops::Log::info() << "ObsBias: cold start skipped for " << varnames[je % nvars] << " ("
                        << n[icold] << " valid departures, minimum is " << minObs
                        << "); coefficients left unchanged" << std::endl;
      continue;
    }
    const auto first = h.begin() + icold * nbins;
    const auto peak = std::max_element(first, first + nbins);
    if (*peak == 0) continue;
    modes[je] = static_cast<double>(peak - first + 1) * binWidth - halfWidth;
  }
  return modes;
}

}  // namespace

// -----------------------------------------------------------------------------

ObsBias::ObsBias(ioda::ObsSpace & odb, const eckit::Configuration & config)
  : numStaticPredictors_(0), numVariablePredictors_(0), byRecord_(),
    vars_(odb.assimvariables()), rank_(odb.distribution()->rank()),
    comm_(odb.comm()), commTime_(odb.commTime()),
    coldStartDone_(false) {
  oops::Log::trace() << "ObsBias::create starting." << std::endl;

  ObsBiasParameters params;
  params.validateAndDeserialize(config);
  byRecord_ = params.BiasCorrectionByRecord;
  const ObsBiasColdStartParameters & csParams = params.coldStart.value();
  coldStartEnabled_ = csParams.enable.value();
  coldStartForce_ = csParams.force.value();
  coldStartBins_ = csParams.bins.value();
  coldStartHalfWidth_ = csParams.halfWidth.value();
  coldStartMinObs_ = csParams.minimumObsNumber.value();

  // Predictor factory
  for (const PredictorParametersWrapper &wrapper :
       params.staticBC.value().predictors.value()) {
    initPredictor(wrapper);
    ++numStaticPredictors_;
  }
  for (const PredictorParametersWrapper &wrapper :
       params.variationalBC.value().predictors.value()) {
    initPredictor(wrapper);
    ++numVariablePredictors_;
  }

  nrecs_ = (byRecord_ && odb.obs_group_vars().size() > 0) ? odb.nrecs() : 1;
  if (byRecord_ && odb.obs_group_vars().size() == 0) {
    throw eckit::BadParameter("ObsBiasParameters: BiasCorrectionByRecord is turned on, "
                              "but the observations are not grouped into records.");
  }
  if (byRecord_ && coldStartEnabled_) {
    // The cold start reduces its histograms over all MPI tasks, but with by-record bias
    // correction each task owns a different, non-overlapping set of records (see the comment in
    // ObsBiasCovariance::linearize). The per-task record count would then size the reduction
    // buffers differently on each rank, so the combination is rejected rather than left to
    // deadlock or return nonsense.
    throw eckit::BadParameter("ObsBiasParameters: 'cold start' is not supported together with "
                              "'bc by record'.");
  }
  ASSERT(nrecs_ > 0);

  oops::ObsVariables varsNoBC = params.variablesNoBC;
  varsNoBC.intersection(vars_);  // Safeguard to make sure that varsNoBC is a subset of vars_
  for (size_t ii = 0; ii < varsNoBC.size(); ++ii) {
    size_t index = vars_.find(varsNoBC[ii]);
    varIndexNoBC_.push_back(index);
  }

  // save record IDs for matching
  if (byRecord_) {
    recIds_.resize(nrecs_);
    // get all ids and obs types (for the hack to be removed)
    std::vector<std::string> allids;
    odb.get_db("MetaData", "stationIdentification", allids);
    // save station ids for all records
    size_t jrec = 0;
    for (auto irec = odb.recidx_begin(); irec != odb.recidx_end(); ++irec, ++jrec) {
      // all the identifiers will be the same for the same record, use the first one
      const size_t iloc = odb.recidx_vector(irec)[0];
      // remove trailing whitespaces (should really be done in files)
      const size_t strEnd = allids[iloc].find_last_not_of(" \t");
      recIds_[jrec] = allids[iloc].substr(0, strEnd+1);
    }
  }

  if (prednames_.size() == 0) {
    oops::Log::info() << "No bias-correction is performed for this ObsSpace." << std::endl;
  } else if (varIndexNoBC_.empty()) {
    oops::Log::info() << "All variables / channels for this ObsSpace are bias-corrected."
                      << std::endl;
  } else if (varsNoBC == vars_) {
    oops::Log::warning() << "None of the variables / channels for this ObsSpace is bias-corrected."
                         << std::endl;
  } else {
    oops::Log::info()
            << "The following variables / channels for this ObsSpace are not bias-corrected: "
            << varsNoBC << std::endl;
  }

  if (vars_.size() * prednames_.size() > 0) {
    // Initialize the coefficients of variable predictors to 0. (Coefficients of static predictors
    // are not stored; they are always equal to 1.)
    biascoeffs_ = Eigen::VectorXd::Zero(nrecs_ * vars_.size() * numVariablePredictors_);
    // Read or initialize bias coefficients
    this->read(config);
  }

  oops::Log::trace() << "ObsBias::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBias::ObsBias(const ObsBias & other, const bool copy)
  : predictors_(other.predictors_),
    prednames_(other.prednames_),
    numStaticPredictors_(other.numStaticPredictors_),
    numVariablePredictors_(other.numVariablePredictors_),
    byRecord_(other.byRecord_),
    nrecs_(other.nrecs_),
    vars_(other.vars_), varIndexNoBC_(other.varIndexNoBC_),
    geovars_(other.geovars_), hdiags_(other.hdiags_), rank_(other.rank_),
    comm_(other.comm_), commTime_(other.commTime_),
    coldStartEnabled_(other.coldStartEnabled_),
    coldStartForce_(other.coldStartForce_),
    coldStartDone_(other.coldStartDone_),
    coldStartBins_(other.coldStartBins_),
    coldStartHalfWidth_(other.coldStartHalfWidth_),
    coldStartMinObs_(other.coldStartMinObs_) {
  oops::Log::trace() << "ObsBias::copy ctor starting." << std::endl;

  // Initialize the biascoeffs
  biascoeffs_ = Eigen::VectorXd::Zero(nrecs_ * vars_.size() * numVariablePredictors_);

  // Copy the bias coeff data
  if (copy && biascoeffs_.size() > 0) *this = other;

  oops::Log::trace() << "ObsBias::copy ctor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBias & ObsBias::operator+=(const ObsBiasIncrement & dx) {
  biascoeffs_ += dx.data();
  return *this;
}

// -----------------------------------------------------------------------------

ObsBias & ObsBias::operator=(const ObsBias & rhs) {
  if (rhs.size() > 0 && this->size() == rhs.size()) {
    biascoeffs_ = rhs.biascoeffs_;
    predictors_  = rhs.predictors_;
    prednames_  = rhs.prednames_;
    numStaticPredictors_ = rhs.numStaticPredictors_;
    numVariablePredictors_ = rhs.numVariablePredictors_;
    inputBiasCoeffs_ = rhs.inputBiasCoeffs_;
    inputPredictors_ = rhs.inputPredictors_;
    inputRecords_ = rhs.inputRecords_;
    byRecord_   = rhs.byRecord_;
    nrecs_      = rhs.nrecs_;
    vars_       = rhs.vars_;
    geovars_    = rhs.geovars_;
    hdiags_     = rhs.hdiags_;
    rank_       = rhs.rank_;
    coldStartEnabled_    = rhs.coldStartEnabled_;
    coldStartForce_      = rhs.coldStartForce_;
    coldStartDone_       = rhs.coldStartDone_;
    coldStartBins_       = rhs.coldStartBins_;
    coldStartHalfWidth_  = rhs.coldStartHalfWidth_;
    coldStartMinObs_     = rhs.coldStartMinObs_;
  }
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBias::read(const eckit::Configuration & config) {
  oops::Log::trace() << "ObsBias::read and initialize from file, starting "<< std::endl;

  Parameters_ params;
  params.validateAndDeserialize(config);

  if (params.inputFile.value() != boost::none) {
    // Open an hdf5 file with bias coefficients, read only
    ioda::Engines::BackendNames  backendName = ioda::Engines::BackendNames::Hdf5File;
    ioda::Engines::BackendCreationParameters backendParams;
    backendParams.fileName = *params.inputFile.value();
    backendParams.action   = ioda::Engines::BackendFileActions::Open;
    backendParams.openMode = ioda::Engines::BackendOpenModes::Read_Only;

    // Create the backend and attach it to an ObsGroup
    // Use the None DataLyoutPolicy for now to accommodate the current file format
    ioda::Group backend = constructBackend(backendName, backendParams);
    ioda::ObsGroup obsgroup = ioda::ObsGroup(backend,
                   ioda::detail::DataLayoutPolicy::generate(
                         ioda::detail::DataLayoutPolicy::Policies::None));

    // setup variables
    std::vector<Eigen::ArrayXXf> allbiascoeffs;
    std::vector<std::string> predictors;

    // loop through list of coefficients, read them, and store in vector
    for (size_t jpred = numStaticPredictors_; jpred < prednames_.size(); ++jpred) {
      ioda::Variable coeffvar = obsgroup.vars["BiasCoefficients/"+prednames_[jpred]];
      Eigen::ArrayXXf biascoeffs;
      coeffvar.readWithEigenRegular(biascoeffs);
      allbiascoeffs.push_back(biascoeffs);
      predictors.push_back(prednames_[jpred]);
    }

    // Read all record names into the Eigen array
    std::vector<std::string> allrecords;
    if (obsgroup.vars.exists("stationIdentification")) {
      ioda::Variable recvar = obsgroup.vars.open("stationIdentification");
      recvar.read<std::string>(allrecords);
    }

    // If by record then store the read in data
    if (byRecord_) {
      inputBiasCoeffs_ = allbiascoeffs;
      inputPredictors_ = predictors;
      inputRecords_ = allrecords;
    }

    // TODO(corymartin-noaa) read in timestamp of last update

    // Find indices of variables/channels that we need in the data read from the file
    // Don't throw an exception if the variable is not in the file if it does not need to be
    // bias-corrected.
    // A variable or channel that is absent from the input file is normally a fatal error. Two
    // cases are tolerated. It may be explicitly excluded from bias correction via
    // `variables without bc`; or a cold start may be enabled, in which case the variable simply
    // has no prior value yet -- the usual situation when a channel is added to an experiment that
    // is already cycling. Its coefficients are left at zero below, and
    // flagEntriesNeedingColdStart() then seeds them from the first-guess departures. Channels
    // that *are* present in the file still load normally, so a genuinely broken file is not
    // masked: only the missing entries are cold-started, and each one is logged.
    const bool tolerateMissingVars = (varIndexNoBC_.size() > 0) || coldStartEnabled_;
    const std::vector<int> var_idx = getRequiredVarOrChannelIndices(obsgroup, vars_,
                                                                    !tolerateMissingVars);
    // sanity check
    const std::vector<std::string> varnames = vars_.variables();
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      if (var_idx[jvar] == -1) {
        const bool excludedFromBC =
            std::find(varIndexNoBC_.begin(), varIndexNoBC_.end(), jvar) != varIndexNoBC_.end();
        ASSERT(excludedFromBC || coldStartEnabled_);
        if (!excludedFromBC) {
          oops::Log::info() << "ObsBias: " << varnames[jvar] << " has no entry in the input file; "
                            << "it will be cold-started from the first-guess departures."
                            << std::endl;
        }
      }
    }
    // Find indices of predictors that we need in the data read from the file
    const std::vector<int> pred_idx = getAllStrIndices(predictors,
                                      prednames_.begin() + numStaticPredictors_, prednames_.end());
    // Determine if the records are in the input file, if not, add it to the list
    std::vector<int> rec_idx;
    if (byRecord_) {
      bool throwexception = false;
      rec_idx = getAllStrIndices(allrecords,
                recIds_.begin(), recIds_.end(), throwexception);
    } else {
      rec_idx.push_back(0);
    }

    // Filter predictors and channels that we need
    for (size_t jpred = 0; jpred < pred_idx.size(); ++jpred) {
      for (size_t jvar = 0; jvar < var_idx.size(); ++jvar) {
        for (size_t jrec = 0; jrec < nrecs_; ++jrec) {
          if (rec_idx[jrec] == -1) {
            // coeffs are set to 0 if record not in input file
            biascoeffs_[index(jrec, jvar, jpred)] = 0.0;
          } else if (var_idx[jvar] == -1) {
            // coeffs are set to 0 if variable not in input file and
            // does not need to be bias corrected
            biascoeffs_[index(jrec, jvar, jpred)] = 0.0;
          } else {
            // use value from input file
            biascoeffs_[index(jrec, jvar, jpred)] =
                        allbiascoeffs[pred_idx[jpred]](rec_idx[jrec], var_idx[jvar]);
          }
        }
      }
    }
  } else {
    if (numVariablePredictors_ > 0)
      oops::Log::warning() << "ObsBias::prior file is NOT available, starting from ZERO"
                           << std::endl;
  }

  oops::Log::trace() << "ObsBias::read and initilization done " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBias::coldStart(const ioda::ObsSpace & odb, const ioda::ObsVector & hofx,
                        const ioda::ObsDataVector<int> & qcFlags) {
  if (!coldStartEnabled_ || coldStartDone_) return;

  // Only the constant term is ever cold-started.
  const auto constant = std::find(prednames_.begin() + numStaticPredictors_, prednames_.end(),
                                  "constant");
  if (constant == prednames_.end()) {
    oops::Log::warning() << "ObsBias: cold start requested, but 'constant' is not among the "
                         << "variational bias predictors; ignoring." << std::endl;
    return;
  }
  const size_t jconst = constant - (prednames_.begin() + numStaticPredictors_);

  // A (record, variable) is cold-started when read() left all of its coefficients at zero: no
  // input file, or no entry for it in the file.
  // `force` bypasses the test and re-derives everything.
  const size_t nvars = vars_.size();
  std::vector<bool> wanted(nrecs_ * nvars, false);
  size_t ncold = 0;
  for (size_t je = 0; je < wanted.size(); ++je) {
    const size_t jrec = je / nvars, jvar = je % nvars;
    if (std::find(varIndexNoBC_.begin(), varIndexNoBC_.end(), jvar) != varIndexNoBC_.end())
      continue;   // excluded from bias correction altogether
    bool allZero = true;
    for (size_t jp = 0; jp < numVariablePredictors_ && allZero; ++jp)
      allZero = (biascoeffs_[index(jrec, jvar, jp)] == 0.0);
    if (coldStartForce_ || allZero) {
      wanted[je] = true;
      ++ncold;
    }
  }
  if (ncold == 0) return;
  oops::Log::info() << "ObsBias: VarBC cold start for " << ncold << " of " << nrecs_ * nvars
                    << " (record, variable) pairs." << std::endl;

  const std::vector<double> modes = departureModes(odb, hofx, qcFlags, wanted, nrecs_,
                                                   coldStartBins_, coldStartHalfWidth_,
                                                   coldStartMinObs_, commTime_);

  const double missing = util::missingValue<double>();
  const std::vector<std::string> varnames = vars_.variables();
  for (size_t jrec = 0; jrec < nrecs_; ++jrec) {
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      const double mode = modes[jrec * nvars + jvar];
      if (mode == missing) continue;   // not wanted, or too few departures (logged already)
      for (size_t jp = 0; jp < numVariablePredictors_; ++jp)
        biascoeffs_[index(jrec, jvar, jp)] = 0.0;
      biascoeffs_[index(jrec, jvar, jconst)] = mode;
      oops::Log::info() << "ObsBias: cold-starting " << varnames[jvar] << " --> " << mode
                        << std::endl;
    }
  }
  coldStartDone_ = true;
}

// -----------------------------------------------------------------------------


void ObsBias::write(const eckit::Configuration & config) const {
  oops::Log::trace() << "ObsBias::write start" << std::endl;

  std::vector<std::string> globalRecordIds;
  std::vector<double> globalBiasCoeffs;
  if (byRecord_) {
    // gather the records from all MPI threads
    globalRecordIds = recIds_;
    oops::mpi::allGatherv(comm_, globalRecordIds);

    // gather the bias coefficients from all MPI threads to the zeroth thread
    const std::vector<double> localcoeffs(
                biascoeffs_.data(), biascoeffs_.data() + biascoeffs_.size());
    oops::mpi::gather(comm_, localcoeffs, globalBiasCoeffs, 0);
  }

  // only write files out on the task with MPI rank 0
  if (rank_ != 0 || commTime_.rank() != 0) return;

  Parameters_ params;
  params.validateAndDeserialize(config);

  if (params.outputFile.value() != boost::none) {
    // Create a file, overwrite if exists
    const std::string output_filename = *params.outputFile.value();
    ioda::Group group = ioda::Engines::HH::createFile(output_filename,
                        ioda::Engines::BackendCreateModes::Truncate_If_Exists);

    // put only variable bias predictors into the predictors vector
    const std::vector<std::string> predictors(prednames_.begin() + numStaticPredictors_,
                                              prednames_.end());

    // map coefficients to 2D for saving
    if (byRecord_) {
      // Get global record indices and work out if there are new records
      bool throwexception = false;
      const std::vector<int> rec_idx = getAllStrIndices(
                  inputRecords_, globalRecordIds.begin(), globalRecordIds.end(), throwexception);
      const int nnewrecs = std::count_if(rec_idx.begin(), rec_idx.end(), [](int x) {
          return x < 0; });

      // Get used predictor indices
      const std::vector<int> pred_idx = getAllStrIndices(predictors,
                  prednames_.begin() + numStaticPredictors_, prednames_.end());

      // Setup matrix for output
      const size_t nrecs = inputRecords_.size() + nnewrecs;
      const size_t npreds = inputPredictors_.size();
      const size_t nvars = vars_.size();
      Eigen::VectorXd finalcoeffs = Eigen::VectorXd::Zero(nrecs * nvars * npreds);
      std::vector<std::string> finalrecords(nrecs);

      // Loop over the matrix and populate with data from the input file
      for (size_t jpred = 0; jpred < npreds; ++jpred) {
        for (size_t jvar = 0; jvar < nvars; ++jvar) {
          for (size_t jrec = 0; jrec < nrecs - nnewrecs; ++jrec) {
            finalcoeffs[index(jrec, jvar, jpred)] = inputBiasCoeffs_[jpred](jrec, jvar);
            if (jpred == 0 && jvar == 0) finalrecords[jrec] = inputRecords_[jrec];
          }
        }
      }

      // Add new records to output list and assign indices
      std::vector<int> outrec_idx = rec_idx;
      int newRecords = 0;
      for (size_t jrec = 0; jrec < outrec_idx.size(); ++jrec) {
        if (outrec_idx[jrec] == -1) {
          // Add to the end of output
          const size_t recindex = nrecs - nnewrecs + newRecords;
          finalrecords[recindex] = globalRecordIds[jrec];
          outrec_idx[jrec] = recindex;
          newRecords += 1;
        }
      }

      // Update coefficients with active records
      for (size_t jpred = 0; jpred < pred_idx.size(); ++jpred) {
        for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
          for (size_t jrec = 0; jrec < outrec_idx.size(); ++jrec) {
              finalcoeffs[index(outrec_idx[jrec], jvar, pred_idx[jpred])] =
                      globalBiasCoeffs[index(jrec, jvar, jpred)];
          }
        }
      }

      // Convert coefficients and send off for writing
      const Eigen::Map<const Eigen::MatrixXd> coeffs(finalcoeffs.data(), npreds, nrecs * nvars);
      saveBiasCoeffsWithRecords(group, inputPredictors_, finalrecords, vars_.variables(), coeffs);
    } else {
      Eigen::Map<const Eigen::MatrixXd>
          coeffs(biascoeffs_.data(), numVariablePredictors_, nrecs_ * vars_.size());
      saveBiasCoeffsWithChannels(group, predictors, vars_.channels(), coeffs);
    }
    // map coefficients to 2D for saving
  } else {
    if (numVariablePredictors_ > 0) {
      oops::Log::warning() << "obs bias.output file is NOT available, bias coefficients "
                           << "will not be saved." << std::endl;
    }
  }
  oops::Log::trace() << "ObsBias::write end" << std::endl;
}

// -----------------------------------------------------------------------------

double ObsBias::norm() const {
  oops::Log::trace() << "ObsBias::norm starting." << std::endl;
  double zz = 0.0;

  // Static predictors
  const int numUnitCoeffs = nrecs_ * vars_.size() * numStaticPredictors_;
  zz += numUnitCoeffs;

  // Variable predictors
  zz += biascoeffs_.squaredNorm();

  // Compute average and take square root
  const int numCoeffs = numUnitCoeffs + biascoeffs_.size();
  if (numCoeffs > 0) zz = std::sqrt(zz/numCoeffs);

  oops::Log::trace() << "ObsBias::norm done." << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

void ObsBias::zero() {
  biascoeffs_ = Eigen::VectorXd::Zero(nrecs_ * vars_.size() * numVariablePredictors_);
}

// -----------------------------------------------------------------------------

std::vector<std::shared_ptr<const PredictorBase>> ObsBias::variablePredictors() const {
  return std::vector<std::shared_ptr<const PredictorBase>>(
    predictors_.begin() + numStaticPredictors_, predictors_.end());
}

// -----------------------------------------------------------------------------

void ObsBias::print(std::ostream & os) const {
  if (this->size() > 0) {
    // map bias coeffs to eigen matrix
    Eigen::Map<const Eigen::MatrixXd>
      coeffs(biascoeffs_.data(), numVariablePredictors_, nrecs_ * vars_.size());
    os << std::endl << "Obs bias coefficients: " << std::endl;
    os << "---------------------------------------------------------------" << std::endl;
    for (std::size_t p = 0; p < numStaticPredictors_; ++p) {
      os << std::setw(20) << prednames_[p]
         << ":  Min= " << std::setw(15) << 1.0f
         << ",  Max= " << std::setw(15) << 1.0f
         << ",  Norm= " << std::setw(15) << std::sqrt(static_cast<double>(nrecs_ * vars_.size()))
         << std::endl;
    }
    for (std::size_t p = 0; p < numVariablePredictors_; ++p) {
      os << std::setw(20) << prednames_[numStaticPredictors_ + p]
         << ":  Min= " << std::setw(15) << coeffs.row(p).minCoeff()
         << ",  Max= " << std::setw(15) << coeffs.row(p).maxCoeff()
         << ",  Norm= " << std::setw(15) << coeffs.row(p).norm()
         << std::endl;
    }
    os << "---------------------------------------------------------------";
  }
}

// -----------------------------------------------------------------------------

void ObsBias::initPredictor(const PredictorParametersWrapper &params) {
  std::shared_ptr<PredictorBase> pred(PredictorFactory::create(params.predictorParameters, vars_));
  predictors_.push_back(pred);
  prednames_.push_back(pred->name());
  geovars_ += pred->requiredGeovars();
  hdiags_ += pred->requiredHdiagnostics();

  // Reserve the space for ObsBiasTerm for predictor
  if (vars_.channels().size() > 0) {
    // At present we can label predictors with either the channel number or the variable
    // name, but not both. So make sure there's only one multi-channel variable.
    ASSERT(vars_.size() == vars_.channels().size());
    for (auto & job : vars_.channels()) {
      hdiags_ += oops::ObsVariables({prednames_.back() + "_" + std::to_string(job)});
    }
  } else {
    for (const std::string & variable : vars_.variables())
      hdiags_ += oops::ObsVariables({prednames_.back() + "_" + variable});
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
