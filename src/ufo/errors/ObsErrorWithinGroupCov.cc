/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/errors/ObsErrorWithinGroupCov.h"

#include <math.h>
#include <vector>

#include "ioda/Engines/EngineUtils.h"
#include "ioda/Engines/HH.h"
#include "ioda/Layout.h"
#include "ioda/ObsGroup.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/utils/GeometryCalculations.h"
#include "ufo/utils/IodaGroupIndices.h"

namespace ufo {

// -----------------------------------------------------------------------------

static ObsErrorMaker<ObsErrorWithinGroupCov> makerWithinGroupCov_("within group covariances");

// -----------------------------------------------------------------------------

namespace {

// -----------------------------------------------------------------------------
// Correlation functions
// -----------------------------------------------------------------------------
double gc99(const double & distnorm) {
  // computes Gaspari-Cohn 99 localization
  // distnorm - normalized distance
  double gc99value = 0.0;
  if (distnorm < 0.5) {
    gc99value = -8.0*pow(distnorm, 5.0)+8.0*pow(distnorm, 4.0)+5.0*pow(distnorm, 3.0)-
                20.0/3.0*pow(distnorm, 2.0)+1.0;
  } else if (distnorm < 1.0) {
    gc99value = 8.0/3.0*pow(distnorm, 5.0)-8.0*pow(distnorm, 4.0)+5.0*pow(distnorm, 3.0)+
                20.0/3.0*pow(distnorm, 2.0)-10.0*distnorm+4.0-1.0/(3.0*distnorm);
  }
  return gc99value;
}

double markov(const double & distnorm, const double & maxnormdist) {
  // computes Markov localization
  // distnorm - normalized distance
  // maxdist - maximum distance for localization
  // Returns 0.0 for distances larger than maxdist
  double markovvalue = 0.0;
  if (distnorm < maxnormdist) {
    markovvalue = std::exp(-distnorm);
  }
  return markovvalue;
}

// -----------------------------------------------------------------------------
// Distance functions - base class and linear distance
// -----------------------------------------------------------------------------
double distance_linear(double loc1, double loc2) {
  return std::abs(loc1 - loc2);
}

double distance_haversine(double lat1, double lon1, double lat2, double lon2) {
  return ufo::haversine(lat1, lon1, lat2, lon2);  // Distance in meters
}

// -----------------------------------------------------------------------------
// Coordinate constructor function dependent on distance function
// -----------------------------------------------------------------------------
ioda::ObsDataVector<float> coord_constructor(const ObsErrorWithinGroupCovParameters & params,
                                             ioda::ObsSpace & obspace) {
  // Check the distance function will work and create the data vector for coordinates
  std::vector<std::string> vars;
  switch (params.distanceFunction.value()) {
    case DistanceFunctions::LINEAR:
      if (!params.var.value().has_value()) {
        throw eckit::BadParameter("ObsErrorWithinGroupCov: 'var' parameter must be set for "
                                  "linear distance function", Here());
      }
      vars = params.var.value().get();
      assert(vars.size() == 1);
      break;
    case DistanceFunctions::HAVERSINE:
      vars = {"latitude", "longitude"};
      break;
    default:
      throw eckit::BadParameter("ObsErrorWithinGroupCov: unimplemented distance function", Here());
  }
  return ioda::ObsDataVector<float>(obspace, oops::ObsVariables(vars), "MetaData");
}

}  // anonymous namespace

ObsErrorWithinGroupCov::ObsErrorWithinGroupCov(const Parameters_ & params,
                                               ioda::ObsSpace & obspace,
                                               const eckit::mpi::Comm &timeComm)
  : ObsErrorBase(timeComm),
    params_(params),
    obspace_(obspace),
    coord_(coord_constructor(params_, obspace)),
    stddev_(obspace, "ObsError")
{
  correlations_.reserve(obspace.nrecs());
  double missing_double = util::missingValue<double>();

  for (auto irec = obspace.recidx_begin(); irec != obspace.recidx_end(); ++irec) {
    std::vector<size_t> rec_idx = obspace.recidx_vector(irec);
    size_t rec_nobs = rec_idx.size();
    Eigen::MatrixXd corr = Eigen::MatrixXd::Identity(rec_nobs, rec_nobs);
    // Only lower triangle is needed
    for (size_t iloc = 0; iloc < rec_nobs; ++iloc) {
      for (size_t jloc = iloc+1; jloc < rec_nobs; ++jloc) {
         // Compute the distance between the two locations
        double distance = missing_double;
        switch (params_.distanceFunction.value()) {
          case DistanceFunctions::LINEAR: {
            double val1 = coord_[0][rec_idx[iloc]];
            double val2 = coord_[0][rec_idx[jloc]];
            distance = distance_linear(val1, val2);
            break;
          }
          case DistanceFunctions::HAVERSINE: {
            double lat1 = coord_["latitude"][rec_idx[iloc]];
            double lon1 = coord_["longitude"][rec_idx[iloc]];
            double lat2 = coord_["latitude"][rec_idx[jloc]];
            double lon2 = coord_["longitude"][rec_idx[jloc]];
            distance = distance_haversine(lat1, lon1, lat2, lon2);  // in meters
            break;
          }
        }
        // Compute correlation value
        switch (params_.correlationFunction.value()) {
          case CorrelationFunctions::GC99:
            corr(jloc, iloc) = gc99(distance / params_.lscale.value());
            break;
          case CorrelationFunctions::MARKOV:
            corr(jloc, iloc) = markov(distance / params_.lscale.value(),
                                      params_.markovLengthscaleFactor.value());
            break;
        }
      }
    }
    // The Markov correlation method does not guarantee a positive definite matrix therefore a
    // check is needed. If the matrix is not positive definite the diagonal is slightly inflated
    // to make it positive definite which has the effect of slightly inflating the observation
    // errors.
    if (params_.correlationFunction.value() == CorrelationFunctions::MARKOV) {
      Eigen::SelfAdjointView<Eigen::MatrixXd, Eigen::Lower> corr_view(corr);
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(corr_view);
      double min_eigenvalue = es.eigenvalues().minCoeff();
      if (min_eigenvalue < 0.0) {
        oops::Log::trace() << "ObsErrorWithinGroupCov::ObsErrorWithinGroupCov "
                          << "Inflating diagonal of correlation matrix to make it positive "
                          << "definite. Minimum eigenvalue was " << min_eigenvalue << std::endl;
        // Make positive definite by adding a small value to the diagonal
        corr += 1.1 * std::abs(min_eigenvalue) * Eigen::MatrixXd::Identity(rec_nobs, rec_nobs);
      }
    }
    correlations_.push_back(corr);
  }
}

// -----------------------------------------------------------------------------

void ObsErrorWithinGroupCov::update(const ioda::ObsVector & obserr) {
  stddev_ = obserr;
}

// -----------------------------------------------------------------------------

void ObsErrorWithinGroupCov::multiply(ioda::ObsVector & dy) const {
  // R * dy = D^{1/2} * C * D^{1/2} * dy
  // where D^{1/2} - diagonal matrix with stddev_ on the diagonal
  //       C - correlations

  // D^{1/2] * dy
  dy *= stddev_;

  // C * D^{1/2} * dy
  multiplyCorrelations(dy);

  // D^{1/2} * C * D^{1/2} * dy
  dy *= stddev_;
}

// -----------------------------------------------------------------------------

void ObsErrorWithinGroupCov::multiplyCorrelations(ioda::ObsVector & dy) const {
  const size_t nvars = dy.nvars();
  const double missing = util::missingValue<double>();

  size_t recnumLocal = 0;
  for (auto irec = obspace_.recidx_begin(); irec != obspace_.recidx_end();
       ++irec, ++recnumLocal) {
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
        for (size_t jloc = iloc; jloc < nused; ++jloc) {
          int ind2 = usedobs_indices[jloc];
          corr(jloc, iloc) = correlations_[recnumLocal](ind2, ind);
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
}

// -----------------------------------------------------------------------------

void ObsErrorWithinGroupCov::inverseMultiply(ioda::ObsVector & dy) const {
  // R^{-1} * dy = D^{-1/2} * C^{-1} * D^{-1/2} * dy
  // where D^{1/2} - diagonal matrix with stddev_ on the diagonal
  //       C - correlations

  // D^{-1/2] * dy
  dy /= stddev_;


  // C^{-1} * D^{-1/2} * dy
  const size_t nvars = dy.nvars();
  const double missing = util::missingValue<double>();

  size_t recnumLocal = 0;
  for (auto irec = obspace_.recidx_begin(); irec != obspace_.recidx_end();
       ++irec, ++recnumLocal) {
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
        for (size_t jloc = iloc; jloc < nused; ++jloc) {
          int ind2 = usedobs_indices[jloc];
          // only need the lower triangle for llt() below; not filling upper triangle
          corr(jloc, iloc) = correlations_[recnumLocal](ind2, ind);
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

void ObsErrorWithinGroupCov::randomize(ioda::ObsVector & dy) const {
  dy.random();
  multiply(dy);
}

// -----------------------------------------------------------------------------

void ObsErrorWithinGroupCov::save(const std::string & name) const {
  stddev_.save(name);
}

// -----------------------------------------------------------------------------

void ObsErrorWithinGroupCov::saveCorrelations(const std::string & filename,
                             size_t recnum, ioda::ObsVector & randomVec) const {
  if (obspace_.recidx_has(recnum)) {
    // find the index of the correlation matrix for this record on the local task
    const std::vector<size_t> recnums = obspace_.recidx_all_recnums();
    const size_t recnumLocal = std::find(recnums.begin(), recnums.end(), recnum) -
                               recnums.begin();
    // Create a file, overwrite if exists
    ioda::Group group = ioda::Engines::HH::createFile(filename,
                        ioda::Engines::BackendCreateModes::Truncate_If_Exists);
    std::vector<size_t> rec_idx = obspace_.recidx_vector(recnum);
    size_t rec_nlocs = rec_idx.size();

    ioda::NewDimensionScales_t dims {ioda::NewDimensionScale<int>("nlocs", rec_nlocs)};
    ioda::ObsGroup ogrp = ioda::ObsGroup::generate(group, dims);

    // save the coordinates used for correlation computation
    for (size_t ivar=0; ivar < coord_.nvars(); ++ivar) {
      ioda::Variable coordVar = ogrp.vars.createWithScales<float>(
                                coord_.varnames()[ivar], {ogrp.vars["nlocs"]});
      std::vector<float> recCoord(rec_nlocs);
      for (size_t jloc = 0; jloc < rec_nlocs; ++jloc) {
        recCoord[jloc] = coord_[ivar][rec_idx[jloc]];
      }
      coordVar.write(recCoord);
    }

    // Set up the creation parameters for the correlation matrix
    ioda::VariableCreationParameters float_params;
    float_params.chunk = true;               // allow chunking
    float_params.compressWithGZIP();         // compress using gzip
    const float missing_value = util::missingValue<float>();
    float_params.setFillValue<float>(missing_value);

    // Create a variable for correlations, save the values to the variable
    ioda::Variable corrVar = ogrp.vars.createWithScales<float>("correlations",
                         {ogrp.vars["nlocs"], ogrp.vars["nlocs"]}, float_params);
    corrVar.writeWithEigenRegular(
            Eigen::MatrixXd(correlations_[recnumLocal].selfadjointView<Eigen::Lower>()));

    // For diagnostics, output the random vector, and the random vector multiplied
    // by the covariance
    const size_t nvars = randomVec.nvars();
    ioda::Variable randVar = ogrp.vars.createWithScales<float>(
                              "randomVector", {ogrp.vars["nlocs"]});
    std::vector<float> randVect(rec_nlocs);
    for (size_t jloc = 0; jloc < rec_nlocs; ++jloc) {
      randVect[jloc] = randomVec[rec_idx[jloc]*nvars];
    }
    randVar.write(randVect);

    multiplyCorrelations(randomVec);
    ioda::Variable randMultVar = ogrp.vars.createWithScales<float>(
                              "randomVectorMultipliedByCov", {ogrp.vars["nlocs"]});
    std::vector<float> randMultVect(rec_nlocs);
    for (size_t jloc = 0; jloc < rec_nlocs; ++jloc) {
      randMultVect[jloc] = randomVec[rec_idx[jloc]*nvars];
    }
    randMultVar.write(randMultVect);
  }
}

// -----------------------------------------------------------------------------

std::unique_ptr<ioda::ObsVector> ObsErrorWithinGroupCov::getObsErrors() const {
  return std::make_unique<ioda::ObsVector>(stddev_);
}

// -----------------------------------------------------------------------------

std::unique_ptr<ioda::ObsVector> ObsErrorWithinGroupCov::getInverseVariance() const {
  std::unique_ptr<ioda::ObsVector> inverseVariance = std::make_unique<ioda::ObsVector>(stddev_);
  *inverseVariance *= stddev_;
  inverseVariance->invert();
  return inverseVariance;
}

// -----------------------------------------------------------------------------
void ObsErrorWithinGroupCov::print(std::ostream & os) const {
  os << "Observation error covariance with correlations within group." << std::endl;
  os << " Obs error stddev: " << stddev_ << std::endl;
  os << " Cross-variable correlations for the first record (lower triangle): " << std::endl;
  os << correlations_[0] << std::endl;
}

// -----------------------------------------------------------------------------


}  // namespace ufo
