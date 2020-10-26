/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsBiasIncrement.h"

#include <algorithm>
#include <iomanip>
#include <set>

#include "eckit/mpi/Comm.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ioda::ObsSpace & odb, const eckit::Configuration & conf)
  : predbases_(0), jobs_(0), odb_(odb), conf_(conf), sensor_() {
  oops::Log::trace() << "ObsBiasIncrement::create starting." << std::endl;

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
    }
  }

  // initialize bias coefficient perturbations
  biascoeffsinc_.resize(prednames_.size()*jobs_.size());
  this->zero();

  oops::Log::trace() << "ObsBiasIncrement::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other, const bool copy)
  : odb_(other.odb_), conf_(other.conf_), predbases_(other.predbases_),
    prednames_(other.prednames_), jobs_(other.jobs_),
    biascoeffsinc_(other.biascoeffsinc_), sensor_(other.sensor_) {
  oops::Log::trace() << "ObsBiasIncrement::copy ctor starting" << std::endl;
  oops::Log::trace() << "ObsBiasIncrement::copy ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other,
                                   const eckit::Configuration & conf)
  : odb_(other.odb_), conf_(conf), predbases_(), prednames_(), jobs_(), sensor_() {
  oops::Log::trace() << "ObsBiasIncrement::copy ctor reconf starting" << std::endl;
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
    }
  }

  // initialize bias coefficient perturbations
  biascoeffsinc_.resize(prednames_.size()*jobs_.size());
  this->zero();

  // Copy the data
  if (biascoeffsinc_.size() > 0) *this = other;

  oops::Log::trace() << "ObsBiasIncrement::copy ctor reconf done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::diff(const ObsBias & b1, const ObsBias & b2) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj]= b1[jj] - b2[jj];
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::zero() {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj]= 0.0;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator=(const ObsBiasIncrement & rhs) {
  if (rhs) {
    predbases_.clear();
    jobs_.clear();

    predbases_     = rhs.predbases_;
    jobs_          = rhs.jobs_;
    sensor_        = rhs.sensor_;
    biascoeffsinc_ = rhs.biascoeffsinc_;
  }
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator+=(const ObsBiasIncrement & rhs) {
  oops::Log::trace() << "ObsBiasIncrement::operator+= startng" << std::endl;

  std::transform(biascoeffsinc_.cbegin(), biascoeffsinc_.cend(),
                 rhs.cbegin(),
                 biascoeffsinc_.begin(),
                 std::plus<double>());

  oops::Log::trace() << "ObsBiasIncrement::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator-=(const ObsBiasIncrement & rhs) {
  oops::Log::trace() << "ObsBiasIncrement::operator-= starting" << std::endl;

  std::transform(biascoeffsinc_.cbegin(), biascoeffsinc_.cend(),
                 rhs.cbegin(),
                 biascoeffsinc_.begin(),
                 std::minus<double>());

  oops::Log::trace() << "ObsBiasIncrement::operator-= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator*=(const double fact) {
  oops::Log::trace() << "ObsBiasIncrement::operator*= starting" << std::endl;

  std::transform(biascoeffsinc_.cbegin(), biascoeffsinc_.cend(),
                 biascoeffsinc_.begin(),
                 [fact](double x){return (fact*x);});

  oops::Log::trace() << "ObsBiasIncrement::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::axpy(const double fact, const ObsBiasIncrement & rhs) {
  oops::Log::trace() << "ObsBiasIncrement::axpy starting" << std::endl;

  std::transform(biascoeffsinc_.cbegin(), biascoeffsinc_.cend(),
                 rhs.cbegin(),
                 biascoeffsinc_.begin(),
                 [fact](double x, double y){return (x += fact * y);});

  oops::Log::trace() << "ObsBiasIncrement::axpy done" << std::endl;
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::dot_product_with(const ObsBiasIncrement & rhs) const {
  oops::Log::trace() << "ObsBiasIncrement::dot_product_with starting" << std::endl;

  auto zz = std::inner_product(biascoeffsinc_.cbegin(), biascoeffsinc_.cend(),
                               rhs.cbegin(), 0.0);

  oops::Log::trace() << "ObsBiasIncrement::dot_product_with done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::norm() const {
  oops::Log::trace() << "ObsBiasIncrement::norm starting" << std::endl;

  auto zz = std::sqrt(dot_product_with(*this));

  oops::Log::trace() << "ObsBiasIncrement::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::computeObsBiasTL(const GeoVaLs & geovals,
                                        const std::vector<ioda::ObsVector> & predData,
                                        ioda::ObsVector & ybiasinc) const {
  oops::Log::trace() << "ObsBiasIncrement::computeObsBiasTL starts." << std::endl;

  if (serialSize() > 0) {
    const std::size_t nlocs  = ybiasinc.nlocs();
    const std::size_t npreds = prednames_.size();
    const std::size_t njobs  = jobs_.size();

    ASSERT(biascoeffsinc_.size() == npreds*njobs);
    ASSERT(predData.size() == npreds);
    ASSERT(ybiasinc.nvars() == njobs);

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

    /* ybiasinc memory layout (nlocs X njobs)
     *     ch1    ch2    ch3     ch4
     * Loc --------------------------
     *  0 | 0      1      2       3
     *  1 | 4      5      6       7
     *  2 | 8      9     10      11 
     *  3 |12     13     14      15
     *  4 |16     17     18      19
     * ...|
    */

    ybiasinc.zero();

    // map bias coeffs to eigen matrix (read only)
    Eigen::Map<const Eigen::MatrixXd> coeffs(biascoeffsinc_.data(), npreds, njobs);

    // For each channel: ( nlocs X 1 ) = ( nlocs X npreds ) * ( npreds X 1 )
    for (std::size_t jch = 0; jch < njobs; ++jch) {
      for (std::size_t jp = 0; jp < npreds; ++jp) {
        const double beta = coeffs(jp, jch);
        /// axpy
        for (std::size_t jl = 0; jl < nlocs; ++jl) {
          ybiasinc[jl*njobs+jch] += predData[jp][jl*njobs+jch] * beta;
        }
      }
    }
  }

  oops::Log::trace() << "ObsBiasIncrement::computeObsBiasTL done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::computeObsBiasAD(GeoVaLs & geovals,
                                        const std::vector<ioda::ObsVector> & predData,
                                        const ioda::ObsVector & ybiasinc) {
  oops::Log::trace() << "ObsBiasIncrement::computeObsBiasAD starts." << std::endl;

  if (this->serialSize() > 0) {
    const std::size_t nlocs  = ybiasinc.nlocs();
    const std::size_t npreds = prednames_.size();
    const std::size_t njobs  = jobs_.size();

    ASSERT(biascoeffsinc_.size() == npreds*njobs);
    ASSERT(predData.size() == npreds);
    ASSERT(ybiasinc.nvars() == njobs);

    // map bias coeffs to eigen matrix (writable)
    Eigen::Map<Eigen::MatrixXd> coeffs(biascoeffsinc_.data(), npreds, njobs);

    std::size_t indx;
    Eigen::VectorXd tmp = Eigen::VectorXd::Zero(nlocs, 1);
    for (std::size_t jch = 0; jch < njobs; ++jch) {
      for (std::size_t jl = 0; jl < nlocs; ++jl) {
        indx = jl*njobs+jch;
        if (ybiasinc[indx] != util::missingValue(ybiasinc[indx])) {
          tmp(jl) = ybiasinc[indx];
        }
      }
      // For each channel: ( npreds X 1 ) = ( npreds X nlocs ) * ( nlocs X 1 )
      for (std::size_t jp = 0; jp < npreds; ++jp) {
        for (std::size_t jl = 0; jl < nlocs; ++jl) {
          coeffs(jp, jch) += predData[jp][jl*njobs+jch] * tmp(jl);
        }
      }

      // zero out for next job
      tmp.setConstant(0.0);
    }

    // Sum across the processros
    odb_.comm().allReduceInPlace(biascoeffsinc_.begin(), biascoeffsinc_.end(), eckit::mpi::sum());
  }

  oops::Log::trace() << "ObsBiasIncrement::computeAD done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::print(std::ostream & os) const {
  if (this->serialSize() > 0) {
    // map bias coeffs to eigen matrix (writable)
    Eigen::Map<const Eigen::MatrixXd>
      coeffs(biascoeffsinc_.data(), prednames_.size(), jobs_.size());

    os << std::endl;
    os << "ObsBiasIncrement::print " << std::endl;
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
    os << "ObsBiasIncrement::print done" << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
