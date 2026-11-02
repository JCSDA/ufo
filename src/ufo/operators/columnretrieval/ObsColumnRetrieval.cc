/*
 * (C) Copyright 2017-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/columnretrieval/ObsColumnRetrieval.h"

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
static ObsOperatorMaker<ObsColumnRetrieval> makerColumnRetrieval_("ColumnRetrieval");
// -----------------------------------------------------------------------------

ObsColumnRetrieval::ObsColumnRetrieval(const ioda::ObsSpace & odb,
                                       const Parameters_ & params)
  : ObsOperatorBase(odb),
    odb_(odb),
    tracerVar_(oops::Variable(params.tracerVariables.value()[0])),
    prsiVar_(oops::Variable("air_pressure_levels")),
    nlayers_(params.nlayers_retrieval.value()),
    isApriori_(params.isApriori.value()),
    isAveragingKernel_(params.isAveragingKernel.value()),
    totalNoVertice_(params.totalNoVertice.value()),
    stretch_(params.stretchVertices.value()),
    modelUnitsCoeff_(params.modelUnitsCoeff.value()),
    tropopause_(params.tropopause.value())
{
  varin_ = oops::Variables({tracerVar_, prsiVar_});
  oops::Log::trace() << "ObsColumnRetrieval constructor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsColumnRetrieval::~ObsColumnRetrieval() {
  oops::Log::trace() << "ObsColumnRetrieval destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsColumnRetrieval::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                     ObsDiagnostics &, const QCFlags_t &) const {
  oops::Log::trace() << "ObsColumnRetrieval::simulateObs start" << std::endl;

  const size_t nlocs = ovec.nlocs();
  const size_t nvars = ovec.nvars();
  const size_t nval  = gv.nlevs(prsiVar_);  // number of model interface levels
  const size_t nvert = nlayers_ + 1;        // number of obs pressure vertices
  checkSingleSimulatedVariable(nvars);

  // --- Read obs data (2D fields flattened, row-major: location x level) ---
  // Read as double so that missing values are converted to the double missing value.

  // Averaging kernel [nlocs x nlayers_]
  std::vector<double> avgKernelFlat(nlayers_ * nlocs, 1.0);
  if (isAveragingKernel_) {
    odb_.get_db("RetrievalAncillaryData", "averagingKernel", avgKernelFlat);
    checkRetrievalDataSize("averagingKernel", avgKernelFlat.size(), nlayers_ * nlocs);
  }

  // Obs pressure vertices [nlocs x nvert]
  std::vector<double> prsiObsFlat(nvert * nlocs, 0.0);
  if (totalNoVertice_) {
    // A single layer spanning the model top to the model bottom
    std::vector<double> prsiTop(nlocs), prsiBottom(nlocs);
    gv.getAtLevel(prsiTop, prsiVar_, 0);
    gv.getAtLevel(prsiBottom, prsiVar_, nval - 1);
    for (size_t iobs = 0; iobs < nlocs; ++iobs) {
      prsiObsFlat[iobs * nvert]            = prsiTop[iobs];
      prsiObsFlat[iobs * nvert + nlayers_] = prsiBottom[iobs];
    }
  } else {
    odb_.get_db("RetrievalAncillaryData", "pressureVertice", prsiObsFlat);
    checkRetrievalDataSize("pressureVertice", prsiObsFlat.size(), nvert * nlocs);
  }

  // A priori term [nlocs]
  std::vector<double> aprioriTerm(nlocs, 0.0);
  if (isApriori_) {
    odb_.get_db("RetrievalAncillaryData", "aprioriTerm", aprioriTerm);
    checkRetrievalDataSize("aprioriTerm", aprioriTerm.size(), nlocs);
  }

  // Per-location work arrays
  std::vector<double> prsiModel(nval);
  std::vector<double> tracerProfile(nval - 1);
  std::vector<double> avgKernel(nlayers_);
  std::vector<double> prsiObs(nvert);

  const double missingDouble = util::missingValue<double>();

  for (size_t iobs = 0; iobs < nlocs; ++iobs) {
    // Gather obs pressure vertices and averaging kernel at this location
    std::copy_n(prsiObsFlat.begin() + iobs * nvert, nvert, prsiObs.begin());
    std::copy_n(avgKernelFlat.begin() + iobs * nlayers_, nlayers_, avgKernel.begin());

    // Skip locations with any missing obs data
    if (hasMissingRetrievalData(avgKernel, prsiObs, aprioriTerm[iobs])) {
      ovec[iobs * nvars] = missingDouble;
      continue;
    }

    // Gather model profiles at this location
    gv.getAtLocation(prsiModel,     prsiVar_,   iobs);
    gv.getAtLocation(tracerProfile, tracerVar_, iobs);

    // Apply unit conversion
    for (double & tracer : tracerProfile) tracer *= modelUnitsCoeff_;

    const double hofx = simulateColumnOb(avgKernel, prsiObs, prsiModel,
                                         tracerProfile, stretch_, tropopause_);

    ovec[iobs * nvars] = hofx + aprioriTerm[iobs];
  }

  oops::Log::trace() << "ObsColumnRetrieval::simulateObs done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsColumnRetrieval::print(std::ostream & os) const {
  os << "ObsColumnRetrieval::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
