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
#include <fstream>
#include <iomanip>
#include <set>

#include "eckit/config/Configuration.h"

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

#include "ufo/ObsBiasIncrement.h"
#include "ufo/utils/IodaGroupIndices.h"
#include "ufo/utils/SaveBiasCoeffs.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBias::ObsBias(ioda::ObsSpace & odb, const eckit::Configuration & config)
  : numStaticPredictors_(0), numVariablePredictors_(0), byRecord_(),
    vars_(odb.assimvariables()), rank_(odb.distribution()->rank()),
    comm_(odb.comm()), commTime_(odb.commTime()) {
  oops::Log::trace() << "ObsBias::create starting." << std::endl;

  ObsBiasParameters params;
  params.validateAndDeserialize(config);
  byRecord_ = params.BiasCorrectionByRecord;

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
    comm_(other.comm_), commTime_(other.commTime_) {
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
    bool throwexception = (varIndexNoBC_.size() == 0) ? true : false;
    const std::vector<int> var_idx = getRequiredVarOrChannelIndices(obsgroup, vars_,
                                                                    throwexception);
    // sanity check
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      if (var_idx[jvar] == -1) {
        ASSERT(std::find(varIndexNoBC_.begin(), varIndexNoBC_.end(), jvar) != varIndexNoBC_.end());
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
      os << std::fixed << std::setw(20) << prednames_[p]
         << ":  Min= " << std::setw(15) << 1.0f
         << ",  Max= " << std::setw(15) << 1.0f
         << ",  Norm= " << std::setw(15) << std::sqrt(static_cast<double>(nrecs_ * vars_.size()))
         << std::endl;
    }
    for (std::size_t p = 0; p < numVariablePredictors_; ++p) {
      os << std::fixed << std::setw(20) << prednames_[numStaticPredictors_ + p]
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
