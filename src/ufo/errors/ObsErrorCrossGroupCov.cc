/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/errors/ObsErrorCrossGroupCov.h"

#include <math.h>
#include <vector>

#include "ioda/Engines/EngineUtils.h"
#include "ioda/Engines/HH.h"
#include "ioda/Layout.h"
#include "ioda/ObsGroup.h"

#include "oops/generic/gc99.h"

#include "ufo/utils/IodaGroupIndices.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsErrorCrossGroupCov::ObsErrorCrossGroupCov(const Parameters_ & params,
                                             ioda::ObsSpace & obspace,
                                             const eckit::mpi::Comm &timeComm)
  : ObsErrorBase(timeComm), obspace_(obspace),
    stddev_(obspace, "ObsError"), vars_(obspace.assimvariables())
{
  varcorrelations_.reserve(obspace.nrecs());
  std::vector<float> vertcoord(obspace.nlocs());
  obspace.get_db("MetaData", params.var, vertcoord);
  size_t recnum = 0;
  for (auto irec = obspace.recidx_begin(); irec != obspace.recidx_end(); ++irec, recnum++) {
    std::vector<size_t> rec_idx = obspace.recidx_vector(irec);
    size_t rec_nobs = rec_idx.size();
    Eigen::MatrixXd corr = Eigen::MatrixXd::Identity(rec_nobs, rec_nobs);
    // Only lower triangle is needed
    for (size_t iloc = 0; iloc < rec_nobs; ++iloc) {
      for (size_t jloc = iloc+1; jloc < rec_nobs; ++jloc) {
         corr(jloc, iloc) = oops::gc99(std::abs(vertcoord[rec_idx[iloc]]-vertcoord[rec_idx[jloc]]) /
                                       params.lscale.value());
      }
    }
    varcorrelations_.push_back(corr);
  }
}

// -----------------------------------------------------------------------------

void ObsErrorCrossGroupCov::update(const ioda::ObsVector & obserr) {
  stddev_ = obserr;
}

// -----------------------------------------------------------------------------

void ObsErrorCrossGroupCov::multiply(ioda::ObsVector & dy) const {
  // R * dy = D^{1/2} * C * D^{1/2} * dy
  // where D^{1/2} - diagonal matrix with stddev_ on the diagonal
  //       C - correlations

  // D^{1/2] * dy
  dy *= stddev_;

  // C * D^{1/2} * dy
  const size_t nvars = dy.nvars();
  const double missing = util::missingValue(double());

  size_t recnum = 0;
  for (auto irec = obspace_.recidx_begin();
       irec != obspace_.recidx_end(); ++irec, recnum++) {
    std::vector<size_t> rec_idx = obspace_.recidx_vector(irec);
    size_t rec_nobs = rec_idx.size();
    // preallocate containers
    std::vector<int> usedobs_indices(rec_nobs);
    Eigen::VectorXd dy_at_rec(rec_nobs);
    Eigen::MatrixXd corr = Eigen::MatrixXd::Identity(rec_nobs, rec_nobs);
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      // find values to be used (the ones that passed QC), and create subset of used values
      // for this record (dy_at_rec) and submatrix of correlations for the used values (corr,
      // only lower triangle is filled)
      size_t nused = 0;
      for (size_t iloc = 0; iloc < rec_nobs; ++iloc) {
        if (dy[rec_idx[iloc]*nvars + jvar] != missing) usedobs_indices[nused++] = iloc;
      }
      for (size_t iloc = 0; iloc < nused; ++iloc) {
        int ind = usedobs_indices[iloc];
        dy_at_rec(iloc) = dy[rec_idx[ind]*nvars + jvar];
        for (size_t jloc = iloc+1; jloc < nused; ++jloc) {
          int ind2 = usedobs_indices[jloc];
          corr(jloc, iloc) = varcorrelations_[recnum](ind2, ind);
        }
      }

      // multiply by C
      dy_at_rec.head(nused) = corr.block(0, 0, nused, nused).selfadjointView<Eigen::Lower>() *
                              dy_at_rec.head(nused);
      // save results in dy
      for (size_t iloc = 0; iloc < nused; ++iloc) {
        dy[rec_idx[usedobs_indices[iloc]]*nvars + jvar] = dy_at_rec[iloc];
      }
    }
  }

  // D^{1/2} * C * D^{1/2} * dy
  dy *= stddev_;
}

// -----------------------------------------------------------------------------

void ObsErrorCrossGroupCov::inverseMultiply(ioda::ObsVector & dy) const {
  // R^{-1} * dy = D^{-1/2} * C^{-1} * D^{-1/2} * dy
  // where D^{1/2} - diagonal matrix with stddev_ on the diagonal
  //       C - correlations

  // D^{-1/2] * dy
  dy /= stddev_;


  // C^{-1} * D^{-1/2} * dy
  const size_t nvars = dy.nvars();
  const double missing = util::missingValue(double());

  size_t recnum = 0;
  for (auto irec = obspace_.recidx_begin();
       irec != obspace_.recidx_end(); ++irec, recnum++) {
    std::vector<size_t> rec_idx = obspace_.recidx_vector(irec);
    size_t rec_nobs = rec_idx.size();
    // preallocate containers
    std::vector<int> usedobs_indices(rec_nobs);
    Eigen::VectorXd dy_at_rec(rec_nobs);
    Eigen::MatrixXd corr = Eigen::MatrixXd::Identity(rec_nobs, rec_nobs);
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      // find values to be used (the ones that passed QC), and create subset of used values
      // for this record (dy_at_rec) and submatrix of correlations for the used values (corr,
      // only lower triangle is filled)
      size_t nused = 0;
      for (size_t iloc = 0; iloc < rec_nobs; ++iloc) {
        if (dy[rec_idx[iloc]*nvars + jvar] != missing) usedobs_indices[nused++] = iloc;
      }
      for (size_t iloc = 0; iloc < nused; ++iloc) {
        int ind = usedobs_indices[iloc];
        dy_at_rec(iloc) = dy[rec_idx[ind]*nvars + jvar];
        for (size_t jloc = iloc+1; jloc < nused; ++jloc) {
          int ind2 = usedobs_indices[jloc];
          // only need the lower triangle for llt() below; not filling upper triangle
          corr(jloc, iloc) = varcorrelations_[recnum](ind2, ind);
        }
      }
      // Multiply by inverse of C, using standard Cholesky decomposition from Eigen library
      // https://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html
      corr.topLeftCorner(nused, nused).llt().solveInPlace(dy_at_rec.head(nused));
      // save results in dy
      for (size_t iloc = 0; iloc < nused; ++iloc) {
        dy[rec_idx[usedobs_indices[iloc]]*nvars + jvar] = dy_at_rec[iloc];
      }
    }
  }

  // D^{-1/2} * C^{-1} * D^{-1/2} * dy
  dy /= stddev_;
}

// -----------------------------------------------------------------------------

void ObsErrorCrossGroupCov::randomize(ioda::ObsVector & dy) const {
  dy.random();
  multiply(dy);
}

// -----------------------------------------------------------------------------

void ObsErrorCrossGroupCov::save(const std::string & name) const {
  stddev_.save(name);
}

// -----------------------------------------------------------------------------

std::unique_ptr<ioda::ObsVector> ObsErrorCrossGroupCov::getObsErrors() const {
  return std::make_unique<ioda::ObsVector>(stddev_);
}

// -----------------------------------------------------------------------------

std::unique_ptr<ioda::ObsVector> ObsErrorCrossGroupCov::getInverseVariance() const {
  std::unique_ptr<ioda::ObsVector> inverseVariance = std::make_unique<ioda::ObsVector>(stddev_);
  *inverseVariance *= stddev_;
  inverseVariance->invert();
  return inverseVariance;
}

// -----------------------------------------------------------------------------
void ObsErrorCrossGroupCov::print(std::ostream & os) const {
  os << "Observation error covariance with correlations within group." << std::endl;
  os << " Obs error stddev: " << stddev_ << std::endl;
  os << " Cross-variable correlations for the first record (lower triangle): " << std::endl;
  os << varcorrelations_[0] << std::endl;
}

// -----------------------------------------------------------------------------


}  // namespace ufo
