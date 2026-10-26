/*
 * (C) Copyright 2017-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/columnretrieval/ObsColumnRetrievalTLAD.h"

#include <algorithm>
#include <ostream>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variable.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/columnretrieval/SatColumnAlgorithm.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsColumnRetrievalTLAD> makerColumnRetrievalTL_("ColumnRetrieval");
// -----------------------------------------------------------------------------

ObsColumnRetrievalTLAD::ObsColumnRetrievalTLAD(const ioda::ObsSpace & odb,
                                               const Parameters_ & params)
  : LinearObsOperatorBase(odb),
    odb_(odb),
    tracerVar_(oops::Variable(params.tracerVariables.value()[0])),
    prsiVar_(oops::Variable("air_pressure_levels")),
    nlayers_(params.nlayers_retrieval.value()),
    isApriori_(params.isApriori.value()),
    isAveragingKernel_(params.isAveragingKernel.value()),
    totalNoVertice_(params.totalNoVertice.value()),
    stretch_(params.stretchVertices.value()),
    modelUnitsCoeff_(params.modelUnitsCoeff.value()),
    tropopause_(params.tropopause.value()),
    nlocs_(0), nval_(0)
{
  varin_ = oops::Variables({tracerVar_});
  oops::Log::trace() << "ObsColumnRetrievalTLAD constructor done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsColumnRetrievalTLAD::~ObsColumnRetrievalTLAD() {
  oops::Log::trace() << "ObsColumnRetrievalTLAD destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsColumnRetrievalTLAD::setTrajectory(const GeoVaLs & gv, ObsDiagnostics &,
                                           const QCFlags_t &) {
  oops::Log::trace() << "ObsColumnRetrievalTLAD::setTrajectory start" << std::endl;

  nlocs_ = odb_.nlocs();
  nval_  = gv.nlevs(prsiVar_);

  // Store model pressure interfaces for all locations
  prsiModel_.resize(nval_ * nlocs_);
  std::vector<double> prsiLoc(nval_);
  for (size_t iobs = 0; iobs < nlocs_; ++iobs) {
    gv.getAtLocation(prsiLoc, prsiVar_, iobs);
    std::copy(prsiLoc.begin(), prsiLoc.end(), prsiModel_.begin() + iobs * nval_);
  }

  // Store averaging kernel (read as double so missing values convert correctly)
  avgKernelObs_.assign(nlayers_ * nlocs_, 1.0);
  if (isAveragingKernel_) {
    odb_.get_db("RetrievalAncillaryData", "averagingKernel", avgKernelObs_);
    checkRetrievalDataSize("averagingKernel", avgKernelObs_.size(), nlayers_ * nlocs_);
  }

  // Store obs pressure vertices
  const size_t nvert = nlayers_ + 1;
  prsiObs_.assign(nvert * nlocs_, 0.0);
  if (totalNoVertice_) {
    // A single layer spanning the model top to the model bottom
    for (size_t iobs = 0; iobs < nlocs_; ++iobs) {
      prsiObs_[iobs * nvert]            = prsiModel_[iobs * nval_];
      prsiObs_[iobs * nvert + nlayers_] = prsiModel_[iobs * nval_ + nval_ - 1];
    }
  } else {
    odb_.get_db("RetrievalAncillaryData", "pressureVertice", prsiObs_);
    checkRetrievalDataSize("pressureVertice", prsiObs_.size(), nvert * nlocs_);
  }

  // The a priori term does not enter the TL/AD, but locations where it is missing
  // have a missing H(x) and are skipped here too.
  std::vector<double> aprioriTerm(nlocs_, 0.0);
  if (isApriori_) {
    odb_.get_db("RetrievalAncillaryData", "aprioriTerm", aprioriTerm);
    checkRetrievalDataSize("aprioriTerm", aprioriTerm.size(), nlocs_);
  }

  // Flag locations with missing obs data
  std::vector<double> prsiModel(nval_);
  std::vector<double> avgKernel(nlayers_);
  std::vector<double> prsiObs(nvert);
  validLoc_.assign(nlocs_, true);
  for (size_t iobs = 0; iobs < nlocs_; ++iobs) {
    getTrajectoryAtLocation(iobs, prsiModel, avgKernel, prsiObs);
    validLoc_[iobs] = !hasMissingRetrievalData(avgKernel, prsiObs, aprioriTerm[iobs]);
  }

  oops::Log::trace() << "ObsColumnRetrievalTLAD::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsColumnRetrievalTLAD::simulateObsTL(const GeoVaLs & gv,
                                           ioda::ObsVector & ovec) const {
  oops::Log::trace() << "ObsColumnRetrievalTLAD::simulateObsTL start" << std::endl;

  const size_t nlocs = ovec.nlocs();
  const size_t nvars = ovec.nvars();
  const size_t nmod  = nval_ - 1;  // number of model layers
  checkSingleSimulatedVariable(nvars);

  std::vector<double> tracerTL(nmod);
  std::vector<double> prsiModel(nval_);
  std::vector<double> avgKernel(nlayers_);
  std::vector<double> prsiObs(nlayers_ + 1);

  const double missingDouble = util::missingValue<double>();

  for (size_t iobs = 0; iobs < nlocs; ++iobs) {
    if (!validLoc_[iobs]) {
      ovec[iobs * nvars] = missingDouble;
      continue;
    }

    gv.getAtLocation(tracerTL, tracerVar_, iobs);
    for (double & tracer : tracerTL) tracer *= modelUnitsCoeff_;

    getTrajectoryAtLocation(iobs, prsiModel, avgKernel, prsiObs);

    ovec[iobs * nvars] = simulateColumnObTL(avgKernel, prsiObs, prsiModel,
                                            tracerTL, stretch_, tropopause_);
  }

  oops::Log::trace() << "ObsColumnRetrievalTLAD::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsColumnRetrievalTLAD::simulateObsAD(GeoVaLs & gv,
                                           const ioda::ObsVector & ovec) const {
  oops::Log::trace() << "ObsColumnRetrievalTLAD::simulateObsAD start" << std::endl;

  const size_t nlocs = ovec.nlocs();
  const size_t nvars = ovec.nvars();
  const size_t nmod  = nval_ - 1;
  checkSingleSimulatedVariable(nvars);

  std::vector<double> tracerAD(nmod);
  std::vector<double> contributionAD(nmod);
  std::vector<double> prsiModel(nval_);
  std::vector<double> avgKernel(nlayers_);
  std::vector<double> prsiObs(nlayers_ + 1);

  const double missingDouble = util::missingValue<double>();

  for (size_t iobs = 0; iobs < nlocs; ++iobs) {
    if (!validLoc_[iobs] || ovec[iobs * nvars] == missingDouble) continue;

    const double hofxAD = ovec[iobs * nvars];

    // Read existing adjoint GeoVaLs to accumulate into
    gv.getAtLocation(tracerAD, tracerVar_, iobs);

    getTrajectoryAtLocation(iobs, prsiModel, avgKernel, prsiObs);

    // Accumulate this location's contribution into a zeroed buffer rather than
    // straight into the GeoVaLs, so that the unit conversion below scales only
    // our own contribution and not any adjoint already deposited by another
    // operator sharing these GeoVaLs.
    contributionAD.assign(nmod, 0.0);
    simulateColumnObAD(avgKernel, prsiObs, prsiModel,
                       hofxAD, contributionAD, stretch_, tropopause_);

    // The TL scales its input perturbation by modelUnitsCoeff_, so the adjoint
    // scales our contribution by the same factor before accumulating it.
    for (size_t j = 0; j < nmod; ++j) tracerAD[j] += contributionAD[j] * modelUnitsCoeff_;

    gv.putAtLocation(tracerAD, tracerVar_, iobs);
  }

  oops::Log::trace() << "ObsColumnRetrievalTLAD::simulateObsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsColumnRetrievalTLAD::getTrajectoryAtLocation(size_t iobs,
                                                     std::vector<double> & prsiModel,
                                                     std::vector<double> & avgKernel,
                                                     std::vector<double> & prsiObs) const {
  std::copy_n(prsiModel_.begin() + iobs * nval_, nval_, prsiModel.begin());
  std::copy_n(avgKernelObs_.begin() + iobs * nlayers_, nlayers_, avgKernel.begin());
  std::copy_n(prsiObs_.begin() + iobs * (nlayers_ + 1), nlayers_ + 1, prsiObs.begin());
}

// -----------------------------------------------------------------------------

void ObsColumnRetrievalTLAD::print(std::ostream & os) const {
  os << "ObsColumnRetrievalTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
