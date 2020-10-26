/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsBias.h"

#include <Eigen/Core>

#include <algorithm>
#include <fstream>
#include <functional>
#include <iomanip>
#include <set>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

#include "ufo/obsbias_io/ObsBiasIO.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBias::ObsBias(ioda::ObsSpace & odb, const eckit::Configuration & conf)
  : predbases_(0), jobs_(0), odb_(odb), conf_(conf), sensor_() {
  oops::Log::trace() << "ObsBias::create starting." << std::endl;

  if (conf_.has("obs bias.sensor")) {
    sensor_ = conf_.getString("obs bias.sensor");
  }

  // Get the jobs(channels)
  if (conf_.has("obs bias.jobs")) {
    const std::set<int> jobs = oops::parseIntSet(conf_.getString("obs bias.jobs"));
    jobs_.assign(jobs.begin(), jobs.end());
  }

  // Predictor factory
  if (conf_.has("obs bias.predictors")) {
    auto confs = conf_.getSubConfigurations("obs bias.predictors");
    predbases_.reserve(confs.size());
    prednames_.reserve(confs.size());
    for (std::size_t j = 0; j < confs.size(); ++j) {
      predbases_.emplace_back(PredictorFactory::create(confs[j], jobs_, sensor_, odb_.comm()));
      prednames_.emplace_back(predbases_.back()->name());
      geovars_ += predbases_.back()->requiredGeovars();
      hdiags_ += predbases_.back()->requiredHdiagnostics();

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
    biascoeffs_.resize(prednames_.size() * jobs_.size());
    std::fill(biascoeffs_.begin(), biascoeffs_.end(), 0.0);

    // Read or initialize bias coefficients
    this->read(conf);
  }

  oops::Log::trace() << "ObsBias::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBias::ObsBias(const ObsBias & other, const bool copy)
  : odb_(other.odb_), conf_(other.conf_), predbases_(other.predbases_),
    prednames_(other.prednames_), jobs_(other.jobs_),
    geovars_(other.geovars_), hdiags_(other.hdiags_),
    biascoeffs_(other.biascoeffs_), sensor_(other.sensor_) {
  oops::Log::trace() << "ObsBias::copy ctor starting" << std::endl;
  oops::Log::trace() << "ObsBias::copy ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsBias & ObsBias::operator+=(const ObsBiasIncrement & dx) {
  oops::Log::trace() << "ObsBias::operator+= starting" << std::endl;

  std::transform(biascoeffs_.cbegin(), biascoeffs_.cend(),
                 dx.cbegin(),
                 biascoeffs_.begin(),
                 std::plus<double>());

  oops::Log::trace() << "ObsBias::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

ObsBias & ObsBias::operator=(const ObsBias & rhs) {
  if (rhs.size() > 0 && this->size() == rhs.size()) {
    conf_       = rhs.conf_;
    biascoeffs_ = rhs.biascoeffs_;
    predbases_  = rhs.predbases_;
    prednames_  = rhs.prednames_;
    jobs_       = rhs.jobs_;
    sensor_     = rhs.sensor_;
    geovars_    = rhs.geovars_;
    hdiags_     = rhs.hdiags_;
  }
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBias::read(const eckit::Configuration & conf) {
  oops::Log::trace() << "ObsBias::read and initialize from file, starting "<< std::endl;

  const auto fileName = conf.getString("obs bias.prior.datain", "");

  const auto root = 0;

  if (!fileName.empty()) {
    if (odb_.comm().rank() == root) {
      ObsBiasIO< Record > biasIO(fileName, std::ios::in);

      for (std::size_t i = 0; i< jobs_.size(); ++i) {
        const auto results =
          biasIO.readByPredictors(sensor_, jobs_[i], prednames_);
        std::copy(results.begin(), results.end(),
                  biascoeffs_.begin() + i*prednames_.size());
      }
    }

    if (odb_.isDistributed()) {
      odb_.comm().broadcast(biascoeffs_, root);
    }
  }

  oops::Log::trace() << "ObsBias::read and initilization done " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBias::write(const eckit::Configuration & conf) const {
  oops::Log::trace() << "ObsBias::write to file starting" << std::endl;

  const auto fileName = conf.getString("obs bias.analysis.dataout", "");

  if (!fileName.empty() && odb_.comm().rank() == 0) {
    ObsBiasIO< Record > biasIO(fileName, std::ios::out);

    std::vector< double > data(prednames_.size());
    for (std::size_t i = 0; i< jobs_.size(); ++i) {
      std::copy(biascoeffs_.begin() + i*prednames_.size(),
                biascoeffs_.begin() + (i + 1)*prednames_.size(),
                data.begin());
      biasIO.addByPredictors(sensor_, jobs_[i], prednames_, data);
    }

    for (auto & pred : predbases_) {
      pred->write(conf, biasIO);
    }

    biasIO.commit();
  }

  oops::Log::trace() << "ObsBias::write to file done" << std::endl;
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
    Eigen::Map<const Eigen::MatrixXd> coeffs(biascoeffs_.data(), npreds, njobs);

    std::vector<double> biasTerm(nlocs, 0.0);
    //  For each channel: ( nlocs X 1 ) =  ( nlocs X npreds ) * (  npreds X 1 )
    for (std::size_t jch = 0; jch < njobs; ++jch) {
      for (std::size_t jp = 0; jp < npreds; ++jp) {
        // axpy
        const double beta = coeffs(jp, jch);
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
  oops::Log::trace() << "ObsBias::norm starting" << std::endl;
  auto zz = std::sqrt(std::inner_product(biascoeffs_.begin(), biascoeffs_.end(),
                                         biascoeffs_.begin(), 0.0));
  oops::Log::trace() << "ObsBias::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

void ObsBias::print(std::ostream & os) const {
  if (this->size() > 0) {
    // map bias coeffs to eigen matrix (writable)
    Eigen::Map<const Eigen::MatrixXd>
      coeffs(biascoeffs_.data(), prednames_.size(), jobs_.size());

    os << std::endl;
    os << "ObsBias::print " << std::endl;
    os << "---------------------------------------------------------------" << std::endl;
    auto njobs = jobs_.size();
    for (std::size_t p = 0; p < prednames_.size(); ++p) {
      os << std::fixed << std::setw(20) << prednames_[p]
         << ":  Min= " << std::setw(15) << std::setprecision(8)
         << coeffs.row(p).minCoeff()
         << ",  Max= " << std::setw(15) << std::setprecision(8)
         << coeffs.row(p).maxCoeff()
         << ",  Norm= " << std::setw(15) << std::setprecision(8)
         << coeffs.row(p).norm()
         << std::endl;
    }
    os << "---------------------------------------------------------------" << std::endl;
    os << "ObsBias::print done" << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
