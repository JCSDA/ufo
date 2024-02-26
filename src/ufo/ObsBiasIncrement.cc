/*
 * (C) Copyright 2018-2021 UCAR
 * (C) Crown Copyright 2024, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsBiasIncrement.h"

#include <iomanip>
#include <memory>

#include "eckit/config/Configuration.h"

#include "ioda/Engines/HH.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include "ufo/ObsBias.h"
#include "ufo/utils/IodaGroupIndices.h"
#include "ufo/utils/SaveBiasCoeffs.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : vars_(odb.assimvariables()), outputFile_(),
    rank_(odb.distribution()->rank()), commTime_(odb.commTime())
{
  oops::Log::trace() << "ufo::ObsBiasIncrement::create starting." << std::endl;

  Parameters_ params;
  params.validateAndDeserialize(config);

  // Predictor factory
  for (const PredictorParametersWrapper &wrapper :
       params.variationalBC.value().predictors.value()) {
    std::unique_ptr<PredictorBase> predictor(PredictorFactory::create(wrapper.predictorParameters,
                                                                      vars_));
    prednames_.push_back(predictor->name());
  }

  nrecs_ = params.BiasCorrectionByRecord ? odb.nrecs() : 1;
  ASSERT(nrecs_ > 0);

  // initialize bias coefficient perturbations
  biascoeffsinc_ = Eigen::VectorXd::Zero(nrecs_ * vars_.size() * prednames_.size());

  if (params.outputFileInc.value() != boost::none) outputFile_ = *params.outputFileInc.value();

  oops::Log::trace() << "ufo::ObsBiasIncrement::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other, const bool copy)
  : prednames_(other.prednames_), nrecs_(other.nrecs_), vars_(other.vars_),
    outputFile_(other.outputFile_), rank_(other.rank_), commTime_(other.commTime_) {
  oops::Log::trace() << "ufo::ObsBiasIncrement::copy ctor starting" << std::endl;

  // Copy the bias coefficients data, or fill in with zeros
  if (copy) {
    biascoeffsinc_ = other.biascoeffsinc_;
  } else {
    biascoeffsinc_ = Eigen::VectorXd::Zero(nrecs_ * vars_.size() * prednames_.size());
  }

  oops::Log::trace() << "ufo::ObsBiasIncrement::copy ctor done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::diff(const ObsBias & b1, const ObsBias & b2) {
  biascoeffsinc_ = b1.data() - b2.data();
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::zero() {
  biascoeffsinc_ = Eigen::VectorXd::Zero(nrecs_ * vars_.size() * prednames_.size());
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator=(const ObsBiasIncrement & rhs) {
  if (rhs) {
    prednames_     = rhs.prednames_;
    nrecs_         = rhs.nrecs_;
    vars_          = rhs.vars_;
    outputFile_    = rhs.outputFile_;
    rank_          = rhs.rank_;
    biascoeffsinc_ = rhs.biascoeffsinc_;
  }
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator+=(const ObsBiasIncrement & rhs) {
  biascoeffsinc_ += rhs.biascoeffsinc_;
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator-=(const ObsBiasIncrement & rhs) {
  biascoeffsinc_ -= rhs.biascoeffsinc_;
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator*=(const double fact) {
  biascoeffsinc_ *= fact;
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::axpy(const double fact, const ObsBiasIncrement & rhs) {
  biascoeffsinc_ += fact * rhs.biascoeffsinc_;
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::dot_product_with(const ObsBiasIncrement & rhs) const {
  return biascoeffsinc_.dot(rhs.biascoeffsinc_);
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::read(const eckit::Configuration & conf) {
  oops::Log::trace() << "ObsBiasIncrement::read and initialize from file, starting "<< std::endl;

  if (conf.has("test input file")) {
    // Open an hdf5 file with bias coefficients, read only
    ioda::Engines::BackendNames  backendName = ioda::Engines::BackendNames::Hdf5File;
    ioda::Engines::BackendCreationParameters backendParams;
    backendParams.fileName = conf.getString("test input file");
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
                                      prednames_.begin(), prednames_.end());
    const std::vector<int> var_idx = getRequiredVarOrChannelIndices(obsgroup, vars_);

    // Filter predictors and channels that we need
    // FIXME: may be possible by indexing allbiascoeffs(pred_idx, chan_idx) when Eigen 3.4
    // is available

    for (size_t jpred = 0; jpred < pred_idx.size(); ++jpred) {
      for (size_t jvar = 0; jvar < var_idx.size(); ++jvar) {
        for (size_t jrec = 0; jrec < nrecs_; ++jrec) {
          // TODO(Someone): fix RHS to handle records
          biascoeffsinc_[index(jrec, jvar, jpred)] = allbiascoeffs(pred_idx[jpred], var_idx[jvar]);
        }
      }
    }
  } else {
    oops::Log::warning() << "ObsBiasIncrement::read is not designed to be called outside tests"
                         << std::endl;
  }

  oops::Log::trace() << "ObsBiasIncrement::read and initilization done " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::write(const eckit::Configuration &) const {
  // only write files out on the task with MPI rank 0
  if (rank_ != 0 || commTime_.rank() != 0) return;

  if (!outputFile_.empty()) {
    // FIXME: only implemented for channels currently
    if (vars_.channels().size() == 0) {
      throw eckit::NotImplemented("ObsBias::write not implemented for variables without channels",
                                  Here());
    }
    // Create a file, overwrite if exists
    ioda::Group group = ioda::Engines::HH::createFile(outputFile_,
                        ioda::Engines::BackendCreateModes::Truncate_If_Exists);

    // put only variable bias predictors into the predictors vector
    std::vector<std::string> predictors(prednames_.begin(), prednames_.end());
    // map coefficients to 2D for saving
    Eigen::Map<const Eigen::MatrixXd>
        coeffs(biascoeffsinc_.data(), prednames_.size(), nrecs_ * vars_.size());
    // TODO(Someone): change to relax the channels assumption and handle records when needed.
    saveBiasCoeffsWithChannels(group, predictors, vars_.channels(), coeffs);
  } else {
    oops::Log::warning() << "obs bias.increment output file is NOT available, bias coefficients "
                         << "will not be saved." << std::endl;
  }
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::norm() const {
  double zz = 0.0;
  if (biascoeffsinc_.size() > 0) {
    zz = biascoeffsinc_.norm()/std::sqrt(biascoeffsinc_.size());
  }
  return zz;
}

// -----------------------------------------------------------------------------

std::vector<double> ObsBiasIncrement::coefficients(size_t jpred) const {
  std::vector<double> coeffs(nrecs_ * vars_.size());
  for (size_t jrec = 0; jrec < nrecs_; ++jrec) {
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      const size_t jrecvar = jrec * vars_.size() + jvar;
      coeffs[jrecvar] = biascoeffsinc_(jrecvar * prednames_.size() + jpred);
    }
  }
  return coeffs;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::updateCoeff(size_t jpred, const std::vector<double> & coeffs) {
  for (size_t jrec = 0; jrec < nrecs_; ++jrec) {
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      const size_t jrecvar = jrec * vars_.size() + jvar;
      biascoeffsinc_[jrecvar * prednames_.size() + jpred] += coeffs[jrecvar];
    }
  }
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::serialize(std::vector<double> & vect) const {
  std::vector<double> vec_obs(biascoeffsinc_.data(),
                              biascoeffsinc_.data() + biascoeffsinc_.size());
  vect.insert(vect.end(), vec_obs.begin(), vec_obs.end());
  oops::Log::trace() << "ObsBiasIncrement::serialize done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::deserialize(const std::vector<double> & vect, std::size_t & index) {
  for (unsigned int jj = 0; jj < biascoeffsinc_.size(); ++jj) {
    biascoeffsinc_[jj] = vect[index];
    ++index;
  }
  oops::Log::trace() << "ObsBiasIncrement::deserialize done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::print(std::ostream & os) const {
  if (this->serialSize() > 0) {
    // map bias coeffs to eigen matrix
    Eigen::Map<const Eigen::MatrixXd>
      coeffs(biascoeffsinc_.data(), prednames_.size(), nrecs_ * vars_.size());
    os << "ufo::ObsBiasIncrement::print " << std::endl;
    os << "---------------------------------------------------------------" << std::endl;
    for (std::size_t p = 0; p < prednames_.size(); ++p) {
      os << std::fixed << std::setw(20) << prednames_[p]
         << ":  Min= " << std::setw(15) << coeffs.row(p).minCoeff()
         << ",  Max= " << std::setw(15) << coeffs.row(p).maxCoeff()
         << ",  Norm= " << std::setw(15) << coeffs.row(p).norm()
         << std::endl;
    }
    os << "---------------------------------------------------------------" << std::endl;
    os << "ufo::ObsBiasIncrement::print done" << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
