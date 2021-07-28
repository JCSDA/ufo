/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/errors/ObsErrorCrossVarCov.h"

#include <vector>

#include "ioda/Engines/Factory.h"
#include "ioda/Engines/HH.h"
#include "ioda/Layout.h"
#include "ioda/ObsGroup.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsErrorCrossVarCov::ObsErrorCrossVarCov(const eckit::Configuration & conf,
                                         ioda::ObsSpace & obspace,
                                         const eckit::mpi::Comm &timeComm)
  : ObsErrorBase(timeComm),
    stddev_(obspace, "ObsError"),
    varcorrelations_(Eigen::MatrixXd::Identity(stddev_.nvars(), stddev_.nvars()))
{
  ObsErrorCrossVarCovParameters options;
  options.deserialize(conf);

  // Open and read error correlations from the hdf5 file
  ioda::Engines::BackendNames  backendName = ioda::Engines::BackendNames::Hdf5File;
  ioda::Engines::BackendCreationParameters backendParams;
  backendParams.fileName = options.inputFile;
  backendParams.action   = ioda::Engines::BackendFileActions::Open;
  backendParams.openMode = ioda::Engines::BackendOpenModes::Read_Only;

  ioda::Group backend = constructBackend(backendName, backendParams);
  ioda::ObsGroup obsgroup = ioda::ObsGroup(backend,
                         ioda::detail::DataLayoutPolicy::generate(
                         ioda::detail::DataLayoutPolicy::Policies::None));

  ioda::Variable corrvar = obsgroup.vars["obserror_correlations"];
  corrvar.readWithEigenRegular(varcorrelations_);
  // Check that the sizes are correct
  const size_t nvars = stddev_.nvars();
  if ((varcorrelations_.rows() != nvars) || (varcorrelations_.cols() != nvars)) {
    std::string errormsg = std::string("Correlation matrix for R, specified in ") +
                options.inputFile.value() + std::string(" should be size ") +
                std::to_string(nvars) + std::string(" by ") + std::to_string(nvars);
    throw eckit::UserError(errormsg, Here());
  }
}

// -----------------------------------------------------------------------------

void ObsErrorCrossVarCov::update(const ioda::ObsVector & obserr) {
  stddev_ = obserr;
}

// -----------------------------------------------------------------------------

void ObsErrorCrossVarCov::multiply(ioda::ObsVector & dy) const {
  // R * dy = D^{1/2} * C * D^{1/2} * dy
  // where D^{1/2} - diagonal matrix with stddev_ on the diagonal
  //       C - correlations

  // D^{1/2] * dy
  dy *= stddev_;

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

  // D^{1/2} * C * D^{1/2} * dy
  dy *= stddev_;
}

// -----------------------------------------------------------------------------

void ObsErrorCrossVarCov::inverseMultiply(ioda::ObsVector & dy) const {
  // R^{-1} * dy = D^{-1/2} * C^{-1} * D^{-1/2} * dy
  // where D^{1/2} - diagonal matrix with stddev_ on the diagonal
  //       C - correlations

  // D^{-1/2] * dy
  dy /= stddev_;

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

  // D^{-1/2} * C^{-1} * D^{-1/2} * dy
  dy /= stddev_;
}

// -----------------------------------------------------------------------------

void ObsErrorCrossVarCov::randomize(ioda::ObsVector & dy) const {
  dy.random();
  multiply(dy);
}

// -----------------------------------------------------------------------------

void ObsErrorCrossVarCov::save(const std::string & name) const {
  stddev_.save(name);
}

// -----------------------------------------------------------------------------

std::unique_ptr<ioda::ObsVector> ObsErrorCrossVarCov::getObsErrors() const {
  return std::make_unique<ioda::ObsVector>(stddev_);
}

// -----------------------------------------------------------------------------

std::unique_ptr<ioda::ObsVector> ObsErrorCrossVarCov::getInverseVariance() const {
  std::unique_ptr<ioda::ObsVector> inverseVariance = std::make_unique<ioda::ObsVector>(stddev_);
  *inverseVariance *= stddev_;
  inverseVariance->invert();
  return inverseVariance;
}

// -----------------------------------------------------------------------------
void ObsErrorCrossVarCov::print(std::ostream & os) const {
  os << "Observation error covariance with cross-variable correlations." << std::endl;
  os << " Obs error stddev: " << stddev_ << std::endl;
  os << " Cross-variable correlations: " << std::endl << varcorrelations_;
}

// -----------------------------------------------------------------------------


}  // namespace ufo
