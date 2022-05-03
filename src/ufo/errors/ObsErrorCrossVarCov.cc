/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/errors/ObsErrorCrossVarCov.h"

#include <math.h>
#include <vector>

#include "ioda/Engines/Factory.h"
#include "ioda/Engines/HH.h"
#include "ioda/Layout.h"
#include "ioda/ObsGroup.h"

#include "ufo/utils/IodaGroupIndices.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsErrorCrossVarCov::ObsErrorCrossVarCov(const Parameters_ & options,
                                         ioda::ObsSpace & obspace,
                                         const eckit::mpi::Comm &timeComm)
  : ObsErrorBase(timeComm),
    stddev_(obspace, "ObsError"), vars_(obspace.assimvariables()),
    varcorrelations_(Eigen::MatrixXd::Identity(stddev_.nvars(), stddev_.nvars()))
{
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

  Eigen::MatrixXf allvarcorrelations;
  if (obsgroup.vars.exists("obserror_correlations") &&
      (!obsgroup.vars.exists("obserror_covariances"))) {
    ioda::Variable corrvar = obsgroup.vars["obserror_correlations"];
    corrvar.readWithEigenRegular(allvarcorrelations);
    // Issue an error if it's not a correlation matrix
    if (!allvarcorrelations.diagonal().isOnes()) {
      throw eckit::BadParameter("ObsErrorCrossVarCov: obserror_correlations matrix read from "
                                "file is not a correlation matrix (values different from one "
                                "found on a diagonal).");
    }
  } else if (obsgroup.vars.exists("obserror_covariances") &&
             (!obsgroup.vars.exists("obserror_correlations"))) {
    ioda::Variable covvar = obsgroup.vars["obserror_covariances"];
    Eigen::MatrixXf allvarcovariances;
    covvar.readWithEigenRegular(allvarcovariances);
    // Convert to correlations
    Eigen::MatrixXf stddevinv(allvarcovariances.diagonal().array().rsqrt().matrix().asDiagonal());
    allvarcorrelations = stddevinv * allvarcovariances * stddevinv;
  } else {
    oops::Log::error() << "One of obserror_correlations or obserror_covariances has to "
                       << "be specified in the input file " << options.inputFile.value()
                       << std::endl;
    throw eckit::BadParameter("One of obserror_correlations or obserror_covariances has "
                              "to be specified in the input file.");
  }

  // Get channel/variables indices in the allvarcorrelations; index == -1 means that correlations
  // don't exist in the file for this channel/variable, and we'll keep them at 0.
  const std::vector<int> var_idx =
         getRequiredVarOrChannelIndices(obsgroup, vars_, false);

  for (size_t ivar = 0; ivar < var_idx.size(); ++ivar) {
    if (var_idx[ivar] < 0) {
      oops::Log::warning() << "ObsErrorCrossVarCov: Obs error correlations not provided for "
                           << "variable " << obspace.obsvariables()[ivar] << " in "
                           << options.inputFile.value() << ", correlations set to zero."
                           << std::endl;
    } else {
      for (size_t jvar = 0; jvar < var_idx.size(); ++jvar) {
        if (var_idx[jvar] >= 0) {
          varcorrelations_(ivar, jvar) = allvarcorrelations(var_idx[ivar], var_idx[jvar]);
        }
      }
    }
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
  os << " Cross-variable correlations: " << std::endl;
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
  }
}

// -----------------------------------------------------------------------------


}  // namespace ufo
