/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "ufo/obsbias/ObsBiasLinearCombinationTLAD.h"

#include "eckit/mpi/Comm.h"

#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ioda/ObsVector.h"

namespace ufo {

static LinearObsBiasMaker<ObsBiasLinearCombinationTLAD>
       makerBiasLinearCombinationTLAD_("LinearCombination");

// -----------------------------------------------------------------------------

ObsBiasLinearCombinationTLAD::ObsBiasLinearCombinationTLAD(const ioda::ObsSpace & odb,
                                                           const eckit::Configuration & conf,
                                                           const std::vector<std::string> & preds,
                                                           const std::vector<int> & jobs)
  : LinearObsBiasBase(odb, conf), odb_(odb), biascoeffsinc_(0), prednames_(preds), jobs_(jobs) {
  oops::Log::trace() << "ObsBiasLinearCombinationTLAD::create starting." << std::endl;

  biascoeffsinc_.resize(prednames_.size()*jobs_.size(), 0.0);

  oops::Log::trace() << "ObsBiasLinearCombinationTLAD::create done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasLinearCombinationTLAD::read(const eckit::Configuration & conf) const {
  oops::Log::trace() << "ObsBiasLinearCombinationTLAD::read to file not implmented " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasLinearCombinationTLAD::write(const eckit::Configuration & conf) const {
  oops::Log::trace() << "ObsBiasLinearCombinationTLAD::write to file not implmented " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasLinearCombinationTLAD::computeObsBiasTL(const GeoVaLs & geovals,
                                                    const Eigen::MatrixXd & predData,
                                                    ioda::ObsVector & ybiasinc) const {
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

  oops::Log::trace() << "ObsBiasLinearCombinationTLAD::computeTL done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasLinearCombinationTLAD::computeObsBiasAD(GeoVaLs & geovals,
                                                    const Eigen::MatrixXd & predData,
                                                    const ioda::ObsVector & ybiasinc) {
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
  mpi_comm().allReduceInPlace(biascoeffsinc_.begin(), biascoeffsinc_.end(), eckit::mpi::sum());

  oops::Log::trace() << "ObsBiasLinearCombinationTLAD::computeAD done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasLinearCombinationTLAD::diff(const ObsBias & b1, const ObsBias & b2) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj]= b1[jj] - b2[jj];
}

// -----------------------------------------------------------------------------

void ObsBiasLinearCombinationTLAD::zero() {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj]= 0.0;
}

// -----------------------------------------------------------------------------

ObsBiasLinearCombinationTLAD & ObsBiasLinearCombinationTLAD::operator=
                               (const ObsBiasIncrement & rhs) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] = rhs[jj];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasLinearCombinationTLAD & ObsBiasLinearCombinationTLAD::operator+=
                               (const ObsBiasIncrement & rhs) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] += rhs[jj];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasLinearCombinationTLAD & ObsBiasLinearCombinationTLAD::operator-=
                               (const ObsBiasIncrement & rhs) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] -= rhs[jj];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasLinearCombinationTLAD & ObsBiasLinearCombinationTLAD::operator*=(const double fact) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] *= fact;
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBiasLinearCombinationTLAD::axpy(const double fact, const ObsBiasIncrement & rhs) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] += fact * rhs[jj];
}

// -----------------------------------------------------------------------------

double ObsBiasLinearCombinationTLAD::dot_product_with(const ObsBiasIncrement & rhs) const {
  double zz = 0.0;
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj) {
    zz += biascoeffsinc_[jj] * rhs[jj];
  }
  return zz;
}

// -----------------------------------------------------------------------------

double ObsBiasLinearCombinationTLAD::norm() const {
  double zz = 0.0;
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj) {
    zz += biascoeffsinc_[jj] * biascoeffsinc_[jj];
  }
  if (biascoeffsinc_.size() > 0) zz = std::sqrt(zz/biascoeffsinc_.size());
  return zz;
}

// -----------------------------------------------------------------------------

void ObsBiasLinearCombinationTLAD::print(std::ostream & os) const {
  os << "ObsBiasLinearCombinationTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
