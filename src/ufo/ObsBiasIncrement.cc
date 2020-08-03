/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsBiasIncrement.h"

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
  : predbases_(0), jobs_(0), odb_(odb), conf_(conf) {
  oops::Log::trace() << "ObsBiasIncrement::create starting." << std::endl;

  /// Get the jobs(channels)
  if (conf_.has("obs bias.jobs")) {
    const std::set<int> jobs = oops::parseIntSet(conf_.getString("obs bias.jobs"));
    jobs_.assign(jobs.begin(), jobs.end());
  }

  /// Predictor factory
  if (conf_.has("obs bias.predictors")) {
    std::vector<eckit::LocalConfiguration> confs;
    conf_.get("obs bias.predictors", confs);
    typedef std::unique_ptr<PredictorBase> predictor;
    for (std::size_t j = 0; j < confs.size(); ++j) {
      predbases_.push_back(predictor(PredictorFactory::create(confs[j], jobs_)));
      prednames_.push_back(predbases_[j]->name());
    }
  }

  /// initialize bias coefficient perturbations
  biascoeffsinc_.resize(prednames_.size()*jobs_.size(), 0.0);

  oops::Log::trace() << "ObsBiasIncrement::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other, const bool copy)
  : odb_(other.odb_), conf_(other.conf_), predbases_(other.predbases_),
    prednames_(other.prednames_), jobs_(other.jobs_) {
  oops::Log::trace() << "ObsBiasIncrement::copy ctor starting" << std::endl;

  /// initialize bias coefficient perturbations
  biascoeffsinc_.resize(prednames_.size()*jobs_.size(), 0.0);

  /// Copy the bias model coeff data
  if (biascoeffsinc_.size() > 0) *this = other;

  oops::Log::trace() << "ObsBiasIncrement::copy ctor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other,
                                   const eckit::Configuration & conf)
  : odb_(other.odb_), conf_(conf), predbases_(), prednames_(), jobs_() {
  oops::Log::trace() << "ObsBiasIncrement::copy ctor starting." << std::endl;
  /// Get the jobs(channels)
  if (conf_.has("obs bias.jobs")) {
    const std::set<int> jobs = oops::parseIntSet(conf_.getString("obs bias.jobs"));
    jobs_.assign(jobs.begin(), jobs.end());
  }

  /// Predictor factory
  if (conf_.has("obs bias.predictors")) {
    std::vector<eckit::LocalConfiguration> confs;
    conf_.get("obs bias.predictors", confs);
    typedef std::unique_ptr<PredictorBase> predictor;
    for (std::size_t j = 0; j < confs.size(); ++j) {
      predbases_.push_back(predictor(PredictorFactory::create(confs[j], jobs_)));
      prednames_.push_back(predbases_[j]->name());
    }
  }

  /// initialize bias coefficient perturbations
  biascoeffsinc_.resize(prednames_.size()*jobs_.size(), 0.0);

  /// Copy the data
  if (biascoeffsinc_.size() > 0) *this = other;

  oops::Log::trace() << "ObsBiasIncrement::copy ctor done." << std::endl;
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
    biascoeffsinc_ = rhs.biascoeffsinc_;
  }
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator+=(const ObsBiasIncrement & rhs) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] += rhs[jj];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator-=(const ObsBiasIncrement & rhs) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] -= rhs[jj];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator*=(const double fact) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] *= fact;
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::axpy(const double fact, const ObsBiasIncrement & rhs) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] += fact * rhs[jj];
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::dot_product_with(const ObsBiasIncrement & rhs) const {
  double zz = 0.0;
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj) {
    zz += biascoeffsinc_[jj] * rhs[jj];
  }
  return zz;
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::norm() const {
  double zz = 0.0;
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj) {
    zz += biascoeffsinc_[jj] * biascoeffsinc_[jj];
  }
  if (biascoeffsinc_.size() > 0) zz = std::sqrt(zz/biascoeffsinc_.size());
  return zz;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::computeObsBiasTL(const GeoVaLs & geovals,
                                        const Eigen::MatrixXd & predData,
                                        ioda::ObsVector & ybiasinc) const {
  oops::Log::trace() << "ObsBiasIncrement::computeObsBiasTL starts." << std::endl;

  if (serialSize() > 0) {
    const std::size_t nlocs  = ybiasinc.nlocs();
    const std::size_t npreds = prednames_.size();
    const std::size_t njobs  = jobs_.size();

    ASSERT(biascoeffsinc_.size() == npreds*njobs);
    ASSERT(predData.rows() == npreds*njobs);
    ASSERT(predData.cols() == nlocs);
    ASSERT(ybiasinc.nvars() == njobs);

    /* predData memory layout (njobs*npreds X nlocs)
     *     Loc     0      1      2       3
     *           --------------------------
     * ch1 pred1 | 0      1      2       4
     *     pred2 | 5      6      7       8
     *     pred3 | 9     10     11      12 
     * ch2 pred1 |13     14     15      16
     *     pred2 |17     18     19      20
     *     ....  |
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

    Eigen::VectorXd tmp;
    for (std::size_t jch = 0; jch < njobs; ++jch) {
      // ( nlocs X 1 ) = ( nlocs X npreds ) * ( npreds X 1 )
      tmp = predData.transpose().block(0, jch*npreds , nlocs, npreds) *
            coeffs.block(0, jch, npreds, 1);
      for (std::size_t jrow = 0; jrow < nlocs; ++jrow) {
        ybiasinc[jrow*njobs+jch] = tmp(jrow);
      }
    }
  }

  oops::Log::trace() << "ObsBiasIncrement::computeObsBiasTL done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::computeObsBiasAD(GeoVaLs & geovals,
                                        const Eigen::MatrixXd & predData,
                                        const ioda::ObsVector & ybiasinc) {
  oops::Log::trace() << "ObsBiasIncrement::computeObsBiasAD starts." << std::endl;

  if (this->serialSize() > 0) {
    const std::size_t nlocs  = ybiasinc.nlocs();
    const std::size_t npreds = prednames_.size();
    const std::size_t njobs  = jobs_.size();

    ASSERT(biascoeffsinc_.size() == npreds*njobs);
    ASSERT(predData.rows() == npreds*njobs);
    ASSERT(predData.cols() == nlocs);
    ASSERT(ybiasinc.nvars() == njobs);

    // map bias coeffs to eigen matrix (writable)
    Eigen::Map<Eigen::MatrixXd> coeffs(biascoeffsinc_.data(), npreds, njobs);

    std::size_t indx;
    Eigen::VectorXd tmp = Eigen::VectorXd::Zero(nlocs, 1);
    for (std::size_t jch = 0; jch < njobs; ++jch) {
      for (std::size_t jrow = 0; jrow < nlocs; ++jrow) {
        indx = jrow*njobs+jch;
        if (ybiasinc[indx] != util::missingValue(ybiasinc[indx])) {
          tmp(jrow) = ybiasinc[indx];
        }
      }
      // ( npreds X 1 ) = ( npreds X nlocs ) * ( nlocs X 1 )
      coeffs.col(jch) += predData.block(jch*npreds, 0, npreds, nlocs) * tmp;

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
