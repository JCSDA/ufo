/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsBias.h"

#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <set>

#include "ioda/Engines/Factory.h"
#include "ioda/Layout.h"
#include "ioda/ObsGroup.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBias::ObsBias(ioda::ObsSpace & odb, const eckit::Configuration & conf)
  : predictors_(0), vars_(odb.obsvariables()) {
  oops::Log::trace() << "ObsBias::create starting." << std::endl;

  // Predictor factory
  if (conf.has("predictors")) {
    std::vector<eckit::LocalConfiguration> confs;
    conf.get("predictors", confs);
    for (std::size_t j = 0; j < confs.size(); ++j) {
      std::shared_ptr<PredictorBase> pred(PredictorFactory::create(confs[j], vars_));
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
  }

  if (prednames_.size() * vars_.size() > 0) {
    // Initialize the biascoeffs to ZERO
    biascoeffs_ = Eigen::MatrixXf::Zero(prednames_.size(), vars_.size());
    // Read or initialize bias coefficients
    this->read(conf);
  }

  oops::Log::trace() << "ObsBias::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBias::ObsBias(const ObsBias & other, const bool copy)
  : predictors_(other.predictors_),
    prednames_(other.prednames_), vars_(other.vars_),
    geovars_(other.geovars_), hdiags_(other.hdiags_) {
  oops::Log::trace() << "ObsBias::copy ctor starting." << std::endl;

  // Initialize the biascoeffs
  biascoeffs_ = Eigen::MatrixXf::Zero(prednames_.size(), vars_.size());

  // Copy the bias coeff data
  if (copy && biascoeffs_.size() > 0) *this = other;

  oops::Log::trace() << "ObsBias::copy ctor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBias & ObsBias::operator+=(const ObsBiasIncrement & dx) {
  for (std::size_t jj = 0; jj < biascoeffs_.size(); ++jj)
    biascoeffs_(jj) += dx[jj];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBias & ObsBias::operator=(const ObsBias & rhs) {
  if (rhs.size() > 0 && this->size() == rhs.size()) {
    biascoeffs_ = rhs.biascoeffs_;
    predictors_  = rhs.predictors_;
    prednames_  = rhs.prednames_;
    vars_       = rhs.vars_;
    geovars_    = rhs.geovars_;
    hdiags_     = rhs.hdiags_;
  }
  return *this;
}

// -----------------------------------------------------------------------------
/// Returns indices of all \p elements_to_look_for in the \p all_elements vector,
/// in the same order that \p elements_to_look_for are in.
/// Throws an exception if at least one of \p elements_to_look_for is missing
/// from \p all_elements.
template<typename T>
std::vector<int> getAllIndices(const std::vector<T> & all_elements,
                                 const std::vector<T> & elements_to_look_for) {
  std::vector<int> result;
  for (const T & element : elements_to_look_for) {
    const auto it = std::find(all_elements.begin(), all_elements.end(), element);
    if (it != all_elements.end()) {
      result.push_back(std::distance(all_elements.begin(), it));
    } else {
      const std::string errormsg = "getAllIndices: Can't find element in the vector";
      throw eckit::BadParameter(errormsg, Here());
    }
  }
  return result;
}

// -----------------------------------------------------------------------------

void ObsBias::read(const eckit::Configuration & conf) {
  oops::Log::trace() << "ObsBias::read and initialize from file, starting "<< std::endl;

  // Bias coefficients input file name
  std::string input_filename;
  if (conf.has("input file")) {
    input_filename = conf.getString("input file");
    // Open an hdf5 file, read only
    ioda::Engines::BackendNames  backendName = ioda::Engines::BackendNames::Hdf5File;
    ioda::Engines::BackendCreationParameters backendParams;
    backendParams.fileName = input_filename;
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
    const std::vector<int> pred_idx = getRequiredPredictorIndices(obsgroup);
    const std::vector<int> var_idx = getRequiredVarOrChannelIndices(obsgroup);

    // Filter predictors and channels that we need
    // FIXME: may be possible by indexing allbiascoeffs(pred_idx, chan_idx) when Eigen 3.4
    // is available
    for (size_t jpred = 0; jpred < pred_idx.size(); ++jpred) {
      for (size_t jvar = 0; jvar < var_idx.size(); ++jvar) {
         biascoeffs_(jpred, jvar) = allbiascoeffs(pred_idx[jpred], var_idx[jvar]);
      }
    }
  } else {
    oops::Log::warning() << "ObsBias::prior file is NOT available, starting from ZERO"
                         << std::endl;
  }

  oops::Log::trace() << "ObsBias::read and initilization done " << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<int> ObsBias::getRequiredPredictorIndices(const ioda::ObsGroup &obsgroup) const {
  ioda::Variable predictorsvar = obsgroup.vars.open("predictors");
  std::vector<std::string> predictors;
  predictorsvar.read<std::string>(predictors);
  return getAllIndices(predictors, prednames_);
}

// -----------------------------------------------------------------------------

std::vector<int> ObsBias::getRequiredVarOrChannelIndices(const ioda::ObsGroup &obsgroup) const {
  if (vars_.channels().empty()) {
    // Read all variables from the file into std vector
    ioda::Variable variablesvar = obsgroup.vars.open("variables");
    std::vector<std::string> variables;
    variablesvar.read<std::string>(variables);

    // Find the indices of the ones we need
    return getAllIndices(variables, vars_.variables());
  } else {
    // At present we can label predictors with either the channel number or the variable
    // name, but not both. So make sure there's only one multi-channel variable.
    ASSERT(vars_.variables().size() == vars_.channels().size());

    // Read all channels from the file into std vector
    ioda::Variable channelsvar = obsgroup.vars.open("channels");
    std::vector<int> channels;
    channelsvar.read<int>(channels);

    // Find the indices of the ones we need
    return getAllIndices(channels, vars_.channels());
  }
}

// -----------------------------------------------------------------------------

void ObsBias::write(const eckit::Configuration & conf) const {
  oops::Log::trace() << "ObsBias::write to file not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

double ObsBias::norm() const {
  oops::Log::trace() << "ObsBias::norm starting." << std::endl;
  double zz = 0.0;
  for (std::size_t jj = 0; jj < biascoeffs_.size(); ++jj) {
    zz += biascoeffs_(jj) * biascoeffs_(jj);
  }
  if (biascoeffs_.size() > 0) zz = std::sqrt(zz/biascoeffs_.size());
  oops::Log::trace() << "ObsBias::norm done." << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

void ObsBias::print(std::ostream & os) const {
  if (this->size() > 0) {
    // map bias coeffs to eigen matrix (writable)
    os << "ObsBias::print " << std::endl;
    os << "---------------------------------------------------------------" << std::endl;
    for (std::size_t p = 0; p < prednames_.size(); ++p) {
      os << std::fixed << std::setw(20) << prednames_[p]
         << ":  Min= " << std::setw(15) << std::setprecision(8)
         << biascoeffs_.row(p).minCoeff()
         << ",  Max= " << std::setw(15) << std::setprecision(8)
         << biascoeffs_.row(p).maxCoeff()
         << ",  Norm= " << std::setw(15) << std::setprecision(8)
         << biascoeffs_.row(p).norm()
         << std::endl;
    }
    os << "---------------------------------------------------------------" << std::endl;
    os << "ObsBias::print done" << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
