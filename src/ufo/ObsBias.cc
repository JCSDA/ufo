/*
 * (C) Copyright 2017-2021 UCAR
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

#include "ioda/Engines/Factory.h"
#include "ioda/Engines/HH.h"
#include "ioda/Layout.h"
#include "ioda/ObsGroup.h"
#include "ioda/ObsSpace.h"

#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/ObsBiasIncrement.h"
#include "ufo/utils/IodaGroupIndices.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBias::ObsBias(ioda::ObsSpace & odb, const ObsBiasParameters & params)
  : numStaticPredictors_(0), numVariablePredictors_(0), vars_(odb.assimvariables()),
    rank_(odb.distribution()->rank()) {
  oops::Log::trace() << "ObsBias::create starting." << std::endl;

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

  if (prednames_.size() * vars_.size() > 0) {
    // Initialize the coefficients of variable predictors to 0. (Coefficients of static predictors
    // are not stored; they are always equal to 1.)
    biascoeffs_ = Eigen::VectorXd::Zero(numVariablePredictors_ * vars_.size());
    // Read or initialize bias coefficients
    this->read(params);
  }

  if (params.channelsNoBC.value() != boost::none) {
    std::set<int> chNoBC = oops::parseIntSet(*params.channelsNoBC.value());
    std::copy(chNoBC.begin(), chNoBC.end(), std::back_inserter(chlistNoBC_));
  } else {
    oops::Log::trace() << "ObsBias::all channels are bias corrected" << std::endl;
  }

  oops::Log::trace() << "ObsBias::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBias::ObsBias(const ObsBias & other, const bool copy)
  : predictors_(other.predictors_),
    prednames_(other.prednames_),
    numStaticPredictors_(other.numStaticPredictors_),
    numVariablePredictors_(other.numVariablePredictors_),
    vars_(other.vars_),
    geovars_(other.geovars_), hdiags_(other.hdiags_), rank_(other.rank_) {
  oops::Log::trace() << "ObsBias::copy ctor starting." << std::endl;

  // Initialize the biascoeffs
  biascoeffs_ = Eigen::VectorXd::Zero(numVariablePredictors_ * vars_.size());

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
    vars_       = rhs.vars_;
    geovars_    = rhs.geovars_;
    hdiags_     = rhs.hdiags_;
    rank_       = rhs.rank_;
  }
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBias::read(const Parameters_ & params) {
  oops::Log::trace() << "ObsBias::read and initialize from file, starting "<< std::endl;

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
         biascoeffs_[index(jpred, jvar)] = allbiascoeffs(pred_idx[jpred], var_idx[jvar]);
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
/// Create ObsGroup with dimensions npredictors = size(predictors) and
/// nchannels = size(channels), variables predictors, channels and
/// bias_cooefficients (npredictors x nchannels)
ioda::ObsGroup saveBiasCoeffsWithChannels(ioda::Group & parent,
                                          const std::vector<std::string> & predictors,
                                          const std::vector<int> & channels,
                                          const Eigen::MatrixXd & coeffs) {
    // dimensions
    ioda::NewDimensionScales_t dims {
        ioda::NewDimensionScale<int>("npredictors", predictors.size()),
        ioda::NewDimensionScale<int>("nchannels", channels.size())
    };
    // new ObsGroup
    ioda::ObsGroup ogrp = ioda::ObsGroup::generate(parent, dims);

    // save the predictors
    ioda::Variable predsVar = ogrp.vars.createWithScales<std::string>(
                              "predictors", {ogrp.vars["npredictors"]});
    predsVar.write(predictors);
    // and the variables
    ioda::Variable chansVar = ogrp.vars.createWithScales<int>("channels", {ogrp.vars["nchannels"]});
    chansVar.write(channels);

    // Set up the creation parameters for the bias coefficients variable
    ioda::VariableCreationParameters float_params;
    float_params.chunk = true;               // allow chunking
    float_params.compressWithGZIP();         // compress using gzip
    float missing_value = util::missingValue(missing_value);
    float_params.setFillValue<float>(missing_value);

    // Create a variable for bias coefficients, save bias coeffs to the variable
    ioda::Variable biasVar = ogrp.vars.createWithScales<float>("bias_coefficients",
                       {ogrp.vars["npredictors"], ogrp.vars["nchannels"]}, float_params);
    biasVar.writeWithEigenRegular(coeffs);
    return ogrp;
}

// -----------------------------------------------------------------------------

void ObsBias::write(const Parameters_ & params) const {
  // only write files out on the task with MPI rank 0
  if (rank_ != 0) return;

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
        coeffs(biascoeffs_.data(), numVariablePredictors_, vars_.size());

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
  const int numUnitCoeffs = numStaticPredictors_ * vars_.size();
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
  biascoeffs_ = Eigen::VectorXd::Zero(numVariablePredictors_ * vars_.size());
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
      coeffs(biascoeffs_.data(), numVariablePredictors_, vars_.size());
    os << "Obs bias coefficients: " << std::endl;
    os << "---------------------------------------------------------------" << std::endl;
    for (std::size_t p = 0; p < numStaticPredictors_; ++p) {
      os << std::fixed << std::setw(20) << prednames_[p]
         << ":  Min= " << std::setw(15) << std::setprecision(8)
         << 1.0f
         << ",  Max= " << std::setw(15) << std::setprecision(8)
         << 1.0f
         << ",  Norm= " << std::setw(15) << std::setprecision(8)
         << std::sqrt(static_cast<double>(vars_.size()))
         << std::endl;
    }
    for (std::size_t p = 0; p < numVariablePredictors_; ++p) {
      os << std::fixed << std::setw(20) << prednames_[numStaticPredictors_ + p]
         << ":  Min= " << std::setw(15) << std::setprecision(8)
         << coeffs.row(p).minCoeff()
         << ",  Max= " << std::setw(15) << std::setprecision(8)
         << coeffs.row(p).maxCoeff()
         << ",  Norm= " << std::setw(15) << std::setprecision(8)
         << coeffs.row(p).norm()
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
