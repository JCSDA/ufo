/*
 * (C) Copyright 2017-2021 UCAR
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

#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

#include "ufo/ObsBiasIncrement.h"
#include "ufo/utils/IodaGroupIndices.h"
#include "ufo/utils/SaveBiasCoeffs.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBias::ObsBias(ioda::ObsSpace & odb, const eckit::Configuration & config)
  : numStaticPredictors_(0), numVariablePredictors_(0), vars_(odb.assimvariables()),
    rank_(odb.distribution()->rank()), commTime_(odb.commTime())
{
  oops::Log::trace() << "ObsBias::create starting." << std::endl;

  ObsBiasParameters params;
  params.validateAndDeserialize(config);

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

  nrecs_ = params.BiasCorrectionByRecord ? odb.nrecs() : 1;
  ASSERT(nrecs_ > 0);

  if (vars_.size() * prednames_.size() > 0) {
    // Initialize the coefficients of variable predictors to 0. (Coefficients of static predictors
    // are not stored; they are always equal to 1.)
    biascoeffs_ = Eigen::VectorXd::Zero(nrecs_ * vars_.size() * numVariablePredictors_);
    // Read or initialize bias coefficients
    this->read(config);
  }

  oops::Variables varsNoBC = params.variablesNoBC;
  varsNoBC.intersection(vars_);  // Safeguard to make sure that varsNoBC is a subset of vars_
  for (size_t ii = 0; ii < varsNoBC.size(); ++ii) {
    size_t index = vars_.find(varsNoBC[ii]);
    varIndexNoBC_.push_back(index);
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

  oops::Log::trace() << "ObsBias::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBias::ObsBias(const ObsBias & other, const bool copy)
  : predictors_(other.predictors_),
    prednames_(other.prednames_),
    numStaticPredictors_(other.numStaticPredictors_),
    numVariablePredictors_(other.numVariablePredictors_),
    nrecs_(other.nrecs_),
    vars_(other.vars_), varIndexNoBC_(other.varIndexNoBC_),
    geovars_(other.geovars_), hdiags_(other.hdiags_), rank_(other.rank_),
    commTime_(other.commTime_) {
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

    // Read all coefficients into the Eigen array
    ioda::Variable coeffvar = obsgroup.vars["bias_coefficients"];
    Eigen::ArrayXXf allbiascoeffs;
    coeffvar.readWithEigenRegular(allbiascoeffs);

    // Find indices of predictors and variables/channels that we need in the data read from the file
    const std::vector<int> pred_idx = getRequiredVariableIndices(obsgroup, "predictors",
                                      prednames_.begin() + numStaticPredictors_, prednames_.end());
    const std::vector<int> var_idx = getRequiredVarOrChannelIndices(obsgroup, vars_);

    // Filter predictors and channels that we need
    // FIXME: may be possible by indexing allbiascoeffs(pred_idx, chan_idx) when Eigen 3.4
    // is available

    for (size_t jpred = 0; jpred < pred_idx.size(); ++jpred) {
      for (size_t jvar = 0; jvar < var_idx.size(); ++jvar) {
        for (size_t jrec = 0; jrec < nrecs_; ++jrec) {
          // TODO(Someone): fix RHS to handle records
          biascoeffs_[index(jrec, jvar, jpred)] = allbiascoeffs(pred_idx[jpred], var_idx[jvar]);
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
  Parameters_ params;
  params.validateAndDeserialize(config);

  // only write files out on the task with MPI rank 0
  if (rank_ != 0 || commTime_.rank() != 0) return;

  if (params.outputFile.value() != boost::none) {
    // FIXME: only implemented for channels currently
    if (vars_.channels().size() == 0) {
      throw eckit::NotImplemented("ObsBias::write not implemented for variables without channels",
                                  Here());
    }
    // Create a file, overwrite if exists
    const std::string output_filename = *params.outputFile.value();
    ioda::Group group = ioda::Engines::HH::createFile(output_filename,
                        ioda::Engines::BackendCreateModes::Truncate_If_Exists);

    // put only variable bias predictors into the predictors vector
    std::vector<std::string> predictors(prednames_.begin() + numStaticPredictors_,
                                        prednames_.end());
    // map coefficients to 2D for saving
    Eigen::Map<const Eigen::MatrixXd>
        coeffs(biascoeffs_.data(), numVariablePredictors_, nrecs_ * vars_.size());
    // TODO(Someone): change to relax the channels assumption and handle records when needed.
    saveBiasCoeffsWithChannels(group, predictors, vars_.channels(), coeffs);
  } else {
    if (numVariablePredictors_ > 0) {
      oops::Log::warning() << "obs bias.output file is NOT available, bias coefficients "
                           << "will not be saved." << std::endl;
    }
  }
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
    os << "Obs bias coefficients: " << std::endl;
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
      hdiags_ += oops::Variables({prednames_.back() + "_" + std::to_string(job)});
    }
  } else {
    for (const std::string & variable : vars_.variables())
      hdiags_ += oops::Variables({prednames_.back() + "_" + variable});
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
