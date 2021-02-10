/*
 * (C) Copyright 2017-2019 UCAR
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
  : predbases_(0), jobs_(0), odb_(odb) {
  oops::Log::trace() << "ObsBias::create starting." << std::endl;

  // Get the jobs(channels)
  if (conf.has("jobs")) {
    const std::set<int> jobs = oops::parseIntSet(conf.getString("jobs"));
    jobs_.assign(jobs.begin(), jobs.end());
  }

  // Predictor factory
  if (conf.has("predictors")) {
    std::vector<eckit::LocalConfiguration> confs;
    conf.get("predictors", confs);
    for (std::size_t j = 0; j < confs.size(); ++j) {
      std::shared_ptr<PredictorBase> pred(PredictorFactory::create(confs[j], jobs_));
      predbases_.push_back(pred);
      prednames_.push_back(pred->name());
      geovars_ += pred->requiredGeovars();
      hdiags_ += pred->requiredHdiagnostics();

      // Reserve the space for ObsBiasTerm for predictor
      if (jobs_.size() > 0) {
        for (auto & job : jobs_) {
          hdiags_ += oops::Variables({prednames_.back() + "_" + std::to_string(job)});
        }
      } else {
        hdiags_ += oops::Variables({prednames_.back()});
      }
    }
  }

  if (prednames_.size() * jobs_.size() > 0) {
    // Initialize the biascoeffs to ZERO
    biascoeffs_ = Eigen::MatrixXf::Zero(prednames_.size(), jobs_.size());
    // Read or initialize bias coefficients
    this->read(conf);
  }

  oops::Log::trace() << "ObsBias::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBias::ObsBias(const ObsBias & other, const bool copy)
  : odb_(other.odb_), predbases_(other.predbases_),
    prednames_(other.prednames_), jobs_(other.jobs_),
    geovars_(other.geovars_), hdiags_(other.hdiags_) {
  oops::Log::trace() << "ObsBias::copy ctor starting." << std::endl;

  // Initialize the biascoeffs
  biascoeffs_ = Eigen::MatrixXf::Zero(prednames_.size(), jobs_.size());

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
    predbases_  = rhs.predbases_;
    prednames_  = rhs.prednames_;
    jobs_       = rhs.jobs_;
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

    // Read all channels from the file into std vector
    ioda::Variable channelsvar = obsgroup.vars.open("channels");
    std::vector<int> channels;
    channelsvar.read<int>(channels);

    // Read all predictors from the file into std vector
    ioda::Variable predictorsvar = obsgroup.vars.open("predictors");
    std::vector<std::string> predictors;
    predictorsvar.read<std::string>(predictors);

    // Find indices of predictors and channels that we need in the data read from the file
    std::vector<int> pred_idx = getAllIndices(predictors, prednames_);
    std::vector<int> chan_idx = getAllIndices(channels, jobs_);

    // Filter predictors and channels that we need
    // FIXME: may be possible by indexing allbiascoeffs(pred_idx, chan_idx) when Eigen 3.4
    // is available
    for (size_t jpred = 0; jpred < pred_idx.size(); ++jpred) {
      for (size_t jchan = 0; jchan < chan_idx.size(); ++jchan) {
         biascoeffs_(jpred, jchan) = allbiascoeffs(pred_idx[jpred], chan_idx[jchan]);
      }
    }
  } else {
    oops::Log::warning() << "ObsBias::prior file is NOT available, starting from ZERO"
                         << std::endl;
  }

  oops::Log::trace() << "ObsBias::read and initilization done " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBias::write(const eckit::Configuration & conf) const {
  oops::Log::trace() << "ObsBias::write to file not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBias::computeObsBias(ioda::ObsVector & ybias,
                             ObsDiagnostics & ydiags,
                             const std::vector<ioda::ObsVector> & predData) const {
  oops::Log::trace() << "ObsBias::computeObsBias starting" << std::endl;

  if (this->size() > 0) {
    const std::size_t nlocs  = ybias.nlocs();
    const std::size_t npreds = prednames_.size();
    const std::size_t njobs  = jobs_.size();

    ASSERT(biascoeffs_.size() == npreds*njobs);
    ASSERT(predData.size() == npreds);
    ASSERT(ybias.nvars() == njobs);

    /* predData memory layout (npreds X nlocs X njobs)
     *       Loc     0      1      2       3
     *             --------------------------
     * pred1 Chan1 | 0      3      6       9
     *       Chan2 | 1      4      7      10
     *       ....  | 2      5      8      11 
     *
     * pred2 Chan1 |12     15     18      21
     *       Chan2 |13     16     19      22
     *       ....  |14     17     20      23
     */

    ybias.zero();

    /* ybias memory layout (nlocs X njobs)
     *     ch1    ch2    ch3     ch4
     * Loc --------------------------
     *  0 | 0      1      2       3
     *  1 | 4      5      6       7
     *  2 | 8      9     10      11 
     *  3 |12     13     14      15
     *  4 |16     17     18      19
     * ...|
     */

    /* map bias coeff to eigen matrix npreds X njobs (read only)
     * bias coeff memory layout (npreds X njobs)
     *        ch1    ch2    ch3     ch4
     *       --------------------------
     * pred1 | 0      1      2       3
     * pred2 | 4      5      6       7
     * pred3 | 8      9     10      11 
     * ....  |
     */
    std::vector<double> biasTerm(nlocs, 0.0);
    //  For each channel: ( nlocs X 1 ) =  ( nlocs X npreds ) * (  npreds X 1 )
    for (std::size_t jch = 0; jch < njobs; ++jch) {
      for (std::size_t jp = 0; jp < npreds; ++jp) {
        // axpy
        const double beta = biascoeffs_(jp, jch);
        for (std::size_t jl = 0; jl < nlocs; ++jl) {
          biasTerm[jl] = predData[jp][jl*njobs+jch] * beta;
          ybias[jl*njobs+jch] += biasTerm[jl];
        }
        // Save ObsBiasTerms (bias_coeff * predictor) for QC
        const std::string varname = predbases_[jp]->name() + "_" + std::to_string(jobs_[jch]);
        if (ydiags.has(varname)) {
          ydiags.save(biasTerm, varname, 1);
        } else {
          oops::Log::error() << varname << " is not reserved in ydiags !" << std::endl;
          ABORT("ObsBiasTerm variable is not reserved in ydiags");
        }
      }
    }
  }

  oops::Log::trace() << "ObsBias::computeObsBias done." << std::endl;
}

// -----------------------------------------------------------------------------
std::vector<ioda::ObsVector> ObsBias::computePredictors(const GeoVaLs & geovals,
                                                        const ObsDiagnostics & ydiags) const {
  const std::size_t npreds = predbases_.size();

  std::vector<ioda::ObsVector> predData(npreds, ioda::ObsVector(odb_));

  for (std::size_t p = 0; p < npreds; ++p) {
    predbases_[p]->compute(odb_, geovals, ydiags, predData[p]);
    predData[p].save(predbases_[p]->name() + "Predictor");
  }

  oops::Log::trace() << "ObsBias::computePredictors done." << std::endl;
  return predData;
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
    auto njobs = jobs_.size();
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
