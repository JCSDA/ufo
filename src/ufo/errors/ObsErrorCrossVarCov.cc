/*
 * (C) Copyright 2021 UCAR.
 * (C) Crown copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/errors/ObsErrorCrossVarCov.h"

#include <vector>

#include "ioda/Engines/EngineUtils.h"
#include "ioda/Engines/HH.h"
#include "ioda/Layout.h"
#include "ioda/ObsGroup.h"

#include "ufo/errors/ObsErrorReconditioner.h"
#include "ufo/utils/IodaGroupIndices.h"

namespace ufo {

// -----------------------------------------------------------------------------

static ObsErrorMaker<ObsErrorCrossVarCov> makerCrossVarCov_("cross variable covariances");

// -----------------------------------------------------------------------------

ObsErrorCrossVarCov::ObsErrorCrossVarCov(const Parameters_ & params,
                                         ioda::ObsSpace & obspace,
                                         const eckit::mpi::Comm &timeComm)
  : ObsErrorBase(timeComm), params_(params),
    stddev_(obspace, "ObsError"), vars_(obspace.assimvariables()),
    varcorrelations_(Eigen::MatrixXd::Identity(stddev_.nvars(), stddev_.nvars())),
    reconditioner_(nullptr)
{
  oops::Log::trace() << "ObsErrorCrossVarCov::ObsErrorCrossVarCov starting" << std::endl;
  // Create reconditioner
  reconditioner_.reset(new ObsErrorReconditioner(params_.reconditioning.value()));
  // Open and read error correlations from the hdf5 file
  ioda::Engines::BackendNames  backendName = ioda::Engines::BackendNames::Hdf5File;
  ioda::Engines::BackendCreationParameters backendParams;
  backendParams.fileName = params.inputFile;
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
                       << "be specified in the input file " << params.inputFile.value()
                       << std::endl;
    throw eckit::BadParameter("One of obserror_correlations or obserror_covariances has "
                              "to be specified in the input file.");
  }

  // Get channel/variables indices in the allvarcorrelations; index == -1 means that correlations
  // don't exist in the file for this channel/variable, and we'll keep them at 0.
  const std::vector<int> var_idx =
         getRequiredVarOrChannelIndices(obsgroup, vars_, false);

  // Report the number of channels/variables without error correlation provided.
  std::size_t count_noErrCorrect = 0;
  std::size_t count_ttl_vars = var_idx.size();
  for (size_t ivar = 0; ivar < var_idx.size(); ++ivar) {
    if (var_idx[ivar] < 0) {
      ++count_noErrCorrect;
      // Print into to the trace log
      oops::Log::trace() << "ObsErrorCrossVarCov: Obs error correlations not provided for "
                         << "variable " << vars_[ivar] << " in "
                         << params.inputFile.value() << ", correlations set to zero.\n";
    } else {
      for (size_t jvar = 0; jvar < var_idx.size(); ++jvar) {
        if (var_idx[jvar] >= 0) {
          varcorrelations_(ivar, jvar) = allvarcorrelations(var_idx[ivar], var_idx[jvar]);
        }
      }
    }
  }

  // Report how many channels/variables are missing error correlation in the log file.
  if (count_noErrCorrect > 0) {
    oops::Log::warning() << "ObsErrorCrossVarCov: Obs error correlations not provided for "
                         << count_noErrCorrect << " of " << count_ttl_vars
                         << " channels/variables in file: " << params.inputFile.value()
                         << ". To see which channels, turn on OOPS_TRACE\n" << std::endl;
  }
  oops::Log::trace() << "ObsErrorCrossVarCov::ObsErrorCrossVarCov finished" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorCrossVarCov::update(const ioda::ObsVector & obserr) {
  stddev_ = obserr;
  if (params_.reconditioning.value().ReconMethod.value() !=
      ufo::ObsErrorReconditionerMethod::NORECONDITIONING) {
    this->recondition(obserr);
  }
}

// -----------------------------------------------------------------------------

void ObsErrorCrossVarCov::recondition(const ioda::ObsVector & mask) {
  const size_t nlocs = mask.nlocs();
  const size_t nvars = mask.nvars();
  const double missing = util::missingValue<double>();
  // preallocate data
  Eigen::MatrixXd avgcorr = Eigen::MatrixXd::Zero(nvars, nvars);
  size_t nused_locs = 0;
  const double dnlocs = static_cast<double>(nlocs);

  // Masking and packing R matrix at each location
  // loop over all observations locations
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    std::vector<int> usedobs_indices(nvars);
    // Calculating variables that are used from mask
    // These are the variables that pass QC.
    size_t nused = 0;
    for (size_t jvar = 0; jvar < nvars; ++jvar)
      if (mask[jloc * nvars + jvar] != missing) usedobs_indices[nused++] = jvar;

    // Initialising correlation matrix for the given location.
    // This will be a fraction of the correlation matrix after
    // the reconditioning has happened.
    Eigen::MatrixXd corr_at_loc(varcorrelations_ / dnlocs);
    if (nused <= 1) {
      oops::Log::trace() << "nused = " << nused
                         << "at jloc = " << jloc
                         << ", skipping reconditioning.\n";
      continue;
    }
    oops::Log::trace() << "\nReconditioning R matrix at jloc = " << jloc << std::endl;
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(nused, nused);

    // loop over all used variables at location
    for (size_t jvar = 0; jvar < nused; ++jvar) {
      const size_t ivar = usedobs_indices[jvar];
      const size_t ind = jloc * nvars + ivar;
      R(jvar, jvar) = stddev_[ind]*stddev_[ind];
      for (size_t jvar2 = jvar + 1; jvar2 < nused; ++jvar2) {
        const size_t ivar2 = usedobs_indices[jvar2];
        const size_t ind2 = jloc * nvars + ivar2;
        R(jvar, jvar2) = varcorrelations_(ivar, ivar2)
                       * stddev_[ind]
                       * stddev_[ind2];
        R(jvar2, jvar) = varcorrelations_(ivar2, ivar)
                       * stddev_[ind2]
                       * stddev_[ind];
      }
    }

    // Recondition the R matrix
    this->reconditioner_->recondition(R);

    // Unpacking the reconditioned matrix
    // into stddev_ and varcorrelations_ members
    // loop over all used variables at location
    for (size_t jvar = 0; jvar < nused; ++jvar) {
      const size_t ivar = usedobs_indices[jvar];
      const size_t ind = jloc * nvars + ivar;
      // Updating stddev_
      const double stddev_ind = std::sqrt(R(jvar, jvar));
      stddev_[ind] = stddev_ind;
      for (size_t jvar2 = jvar + 1; jvar2 < nused; ++jvar2) {
        // Ensuring location independent correlations
        // by using average of reconditioned correlations at each location
        const double stddev_ind2 = std::sqrt(R(jvar2, jvar2));
        const size_t ivar2 = usedobs_indices[jvar2];
        corr_at_loc(ivar, ivar2) = R(jvar, jvar2)
                                 / (stddev_ind
                                 * stddev_ind2
                                 * dnlocs);
        corr_at_loc(ivar2, ivar) = R(jvar2, jvar)
                                 / (stddev_ind2
                                 * stddev_ind
                                 * dnlocs);
      }
    }
    nused_locs++;
    avgcorr += corr_at_loc;
  }

  // Reassigning the correlations to the renormalised
  // average correlations, if any locations are used
  if (nused_locs > 0) {
    varcorrelations_ = avgcorr
                     * dnlocs
                     / static_cast<double>(nused_locs);
  }
}  // recondition

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
  const double missing = util::missingValue<double>();
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
  const double missing = util::missingValue<double>();
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

void ObsErrorCrossVarCov::localize(ioda::ObsVector & locvector) const {
  oops::Log::trace() << "ufo::ObsErrorCrossVarCov::localize start" << std::endl;

  const double missing = util::missingValue<double>();

  ASSERT_MSG(locvector.size() == stddev_.size(),
             "ObsErrorCrossVarCov::localize: "
             "Localization vector dimension must match dimension of obs space.");
  std::vector<double> localstddev;
  const size_t nlocs = locvector.nlocs();
  const size_t nvars = locvector.nvars();
  local_jvars_.clear();
  local_nobs_.clear();
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    size_t nused = 0;
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      size_t jj = jloc*nvars + jvar;
      if (locvector[jj] != missing && locvector[jj] <= 0) {
        throw eckit::BadValue("Localization weights must be positive. Use "
                              "oops::util::missingValue<double>() to indicate "
                              "an observation with a weight of zero.");
      }
      if (locvector[jj] != missing && stddev_[jj] != missing) {
        localstddev.push_back(stddev_[jj] * std::pow(locvector[jj], -0.5));
        local_jvars_.push_back(jvar);
        ++nused;
      }
    }
    if (nused > 0) local_nobs_.push_back(nused);
  }
  local_stddev_ = Eigen::Map<Eigen::VectorXd>(localstddev.data(), localstddev.size());
}

int ObsErrorCrossVarCov::localDim() const {
  return local_stddev_.size();
}

Eigen::MatrixXf ObsErrorCrossVarCov::localInverseMultiply(const Eigen::MatrixXf & zz) const {
  ASSERT_MSG(zz.cols() == localDim(),
             "ObsErrorCrossVarCov::localInverseMultiply: "
             "Input vector dimension must match local R matrix dimension");
  Eigen::MatrixXf zzRinv(zz.rows(), zz.cols());
  size_t j_obs = 0;
  // loop over obs locations in local volume
  for (int pt_nobs : local_nobs_) {
    // construct pt std dev and correlations for obs at this point
    Eigen::VectorXf pt_stddev = local_stddev_.segment(j_obs, pt_nobs).cast<float>();
    Eigen::MatrixXf pt_corr = Eigen::MatrixXf::Identity(pt_nobs, pt_nobs);
    for (size_t ipt = 0; ipt < pt_nobs; ++ipt) {
      for (size_t jpt = 0; jpt < ipt; ++jpt) {
        // only fill lower triangle as that's all we need for LLT
        pt_corr(ipt, jpt) = static_cast<float>(varcorrelations_(local_jvars_[j_obs+ipt],
                                                                local_jvars_[j_obs+jpt]));
      }
    }
    // Cholesky decomposition to allow inverse multiplication
    Eigen::LLT<Eigen::MatrixXf> LLT(pt_corr);
    for (size_t ii = 0; ii < zz.rows(); ++ii) {
      Eigen::VectorXf row = zz(ii, Eigen::seq(j_obs, j_obs+pt_nobs-1)).transpose();
      // D^{-1/2} * row
      row.array() /= pt_stddev.array();
      // C^{-1} * D^{-1/2} * row
      LLT.solveInPlace(row);
      // D^{-1/2} * C^{-1} * D^{-1/2} * row
      row.array() /= pt_stddev.array();
      // store in output
      zzRinv(ii, Eigen::seq(j_obs, j_obs+pt_nobs-1)) = row.transpose();
    }
    j_obs += pt_nobs;
  }
  return zzRinv;
}

Eigen::VectorXd ObsErrorCrossVarCov::local_invVarR() const {
  throw eckit::BadParameter("ObsErrorCrossVarCov::local_invVarR not implemented.");
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
  os << "Observation error covariance with cross-variable correlations.\n";
  os << " Obs error stddev: " << stddev_ << std::endl;
  os << " Cross-variable correlations: \n";
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
