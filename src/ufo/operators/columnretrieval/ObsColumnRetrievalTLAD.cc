/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/columnretrieval/ObsColumnRetrievalTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsColumnRetrievalTLAD> makerColumnRetrievalTL_("ColumnRetrieval");
// -----------------------------------------------------------------------------

ObsColumnRetrievalTLAD::ObsColumnRetrievalTLAD(const ioda::ObsSpace & odb,
                               const Parameters_ & parameters)
  : LinearObsOperatorBase(odb), keyOperColumnRetrieval_(0), varin_()
{
  ufo_columnretrieval_tlad_setup_f90(keyOperColumnRetrieval_, parameters.toConfiguration(),
                               odb.assimvariables(), varin_);
  oops::Log::trace() << "ObsColumnRetrievalTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsColumnRetrievalTLAD::~ObsColumnRetrievalTLAD() {
  ufo_columnretrieval_tlad_delete_f90(keyOperColumnRetrieval_);
  oops::Log::trace() << "ObsColumnRetrievalTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsColumnRetrievalTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                           const QCFlags_t & qc_flags) {
  oops::Log::trace() << "ObsColumnRetrievalTLAD::setTrajectory entering" << std::endl;

  ufo_columnretrieval_tlad_settraj_f90(keyOperColumnRetrieval_, geovals.toFortran(), obsspace());

  oops::Log::trace() << "ObsColumnRetrievalTLAD::setTrajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsColumnRetrievalTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                           const QCFlags_t & qc_flags) const {
  ufo_columnretrieval_simobs_tl_f90(keyOperColumnRetrieval_, geovals.toFortran(), obsspace(),
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsColumnRetrievalTLAD::simulateObsTL exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsColumnRetrievalTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                           const QCFlags_t & qc_flags) const {
  ufo_columnretrieval_simobs_ad_f90(keyOperColumnRetrieval_, geovals.toFortran(), obsspace(),
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsColumnRetrievalTLAD::simulateObsAD exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsColumnRetrievalTLAD::print(std::ostream & os) const {
  os << "ObsColumnRetrievalTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
