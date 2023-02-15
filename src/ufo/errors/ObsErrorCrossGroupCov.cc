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
  : ObsErrorBase(timeComm),
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
    for (size_t iloc = 0; iloc < rec_nobs; ++iloc) {
      for (size_t jloc = iloc+1; jloc < rec_nobs; ++jloc) {
         corr(iloc, jloc) = oops::gc99(std::abs(vertcoord[rec_idx[iloc]]-vertcoord[rec_idx[jloc]]) /
                                       params.lscale.value());
      }
    }
    oops::Log::info() << corr << std::endl;
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

/*
  // C * D^{1/2} * dy
  const size_t nlocs = dy.nlocs();
  const size_t nvars = dy.nvars();
  const double missing = util::missingValue(double());
  // preallocate data
  std::vector<int> usedobs_indices(nvars);
  Eigen::VectorXd dy_at_loc(nvars);
  Eigen::MatrixXd corr = Eigen::MatrixXd::Identity(nvars, nvars);

  // loop over all observations locations
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    // find values to be used (the ones that passed QC), and create subset of used values
    // at this location (dy_at_loc) and submatrix of correlations for the used values (corr)
    size_t nused = 0;
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      if (dy[jloc*nvars + jvar] != missing) usedobs_indices[nused++] = jvar;
    }
    for (size_t jvar = 0; jvar < nused; ++jvar) {
      int ind = usedobs_indices[jvar];
      dy_at_loc(jvar) = dy[jloc*nvars + ind];
      for (size_t jvar2 = jvar+1; jvar2 < nused; ++jvar2) {
        int ind2 = usedobs_indices[jvar2];
        corr(jvar, jvar2) = varcorrelations_(ind, ind2);
        corr(jvar2, jvar) = varcorrelations_(ind2, ind);
      }
    }
    // multiply by C
    dy_at_loc.head(nused) = corr.block(0, 0, nused, nused) * dy_at_loc.head(nused);
    // save results in dy
    for (size_t jvar = 0; jvar < nused; ++jvar) {
      dy[jloc*nvars + usedobs_indices[jvar]] = dy_at_loc[jvar];
    }
  }
*/
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

/*
  // C^{-1} * D^{-1/2} * dy
  const size_t nlocs = dy.nlocs();
  const size_t nvars = dy.nvars();
  const double missing = util::missingValue(double());
  // preallocate data
  std::vector<int> usedobs_indices(nvars);
  Eigen::VectorXd dy_at_loc(nvars);
  Eigen::MatrixXd corr = Eigen::MatrixXd::Identity(nvars, nvars);
  // loop over all observations locations
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    // find values to be used (the ones that passed QC), and create subset of used values
    // at this location (dy_at_loc) and submatrix of correlations for the used values (corr)
    size_t nused = 0;
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      if (dy[jloc*nvars + jvar] != missing) usedobs_indices[nused++] = jvar;
    }
    for (size_t jvar = 0; jvar < nused; ++jvar) {
      int ind = usedobs_indices[jvar];
      dy_at_loc(jvar) = dy[jloc*nvars + ind];
      for (size_t jvar2 = jvar+1; jvar2 < nused; ++jvar2) {
        int ind2 = usedobs_indices[jvar2];
        // only need the lower triangle for llt() below; not filling upper triangle
        corr(jvar2, jvar) = varcorrelations_(ind2, ind);
      }
    }
    // Multiply by inverse of C, using standard Cholesky decomposition from Eigen library
    // https://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html
    corr.topLeftCorner(nused, nused).llt().solveInPlace(dy_at_loc.head(nused));
    // save results in dy
    for (size_t jvar = 0; jvar < nused; ++jvar) {
      dy[jloc*nvars + usedobs_indices[jvar]] = dy_at_loc[jvar];
    }
  }
*/
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
/*  os << " Cross-variable correlations: " << std::endl;
  if (varcorrelations_.rows() < 8) {
    os << varcorrelations_;
  } else {
     Eigen::MatrixXf::Index maxRow, maxCol;
     float max = varcorrelations_.maxCoeff(&maxRow, &maxCol);
     Eigen::MatrixXf::Index minRow, minCol;
     float min = varcorrelations_.minCoeff(&minRow, &minCol);
     os << "  Maximum correlation: " << max <<  ", between: " <<
            vars_[maxRow] << " and " << vars_[maxCol] << std::endl;
     os << "  Minimum correlation: " << min << ", between: " <<
            vars_[minRow] << " and " << vars_[minCol];
  }*/
}

// -----------------------------------------------------------------------------


}  // namespace ufo
