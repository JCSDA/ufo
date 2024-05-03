/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/marine/insitutemperature/ObsInsituTemperatureTLAD.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsInsituTemperatureTLAD>
   makerInsituTemperatureTL_("InsituTemperature");
// -----------------------------------------------------------------------------

ObsInsituTemperatureTLAD::ObsInsituTemperatureTLAD(const ioda::ObsSpace & odb,
                                                   const ObsInsituTemperatureParameters &
                                                   params)
  : LinearObsOperatorBase(odb), keyOper_(0), varin_()
{
  // get optional list of variables to operate on (to be consistent with other
  // operators that work under "composite") eventhough we assume on the Fortran
  // side that this operator will operate on just a single variable (sea_water_temperature)
  std::vector<int> operatorVarIndices;
  getOperatorVariables(params.variables.value(), odb.assimvariables(),
    operatorVars_, operatorVarIndices);

  // sanity check to make sure waterTemperature is the ONLY variable
  ASSERT_MSG(
    operatorVars_.size() == 1 && operatorVars_[0] == "waterTemperature",
    "InsituTemperature can only work on variable \"waterTemperature\"");

  ufo_insitutemperature_tlad_setup_f90(keyOper_, params.toConfiguration(),
    operatorVars_, operatorVarIndices.data(), operatorVarIndices.size(), varin_);

  oops::Log::trace() << "ObsInsituTemperatureTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsInsituTemperatureTLAD::~ObsInsituTemperatureTLAD() {
  ufo_insitutemperature_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsInsituTemperatureTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperatureTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                             const QCFlags_t & qc_flags) {
  ufo_insitutemperature_tlad_settraj_f90(keyOper_, geovals.toFortran(), obsspace());
  oops::Log::trace() << "ObsInsituTemperatureTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperatureTLAD::simulateObsTL(const GeoVaLs & geovals,
                                             ioda::ObsVector & ovec,
                                             const QCFlags_t & qc_flags) const {
  ufo_insitutemperature_simobs_tl_f90(keyOper_, geovals.toFortran(), obsspace(),
                                      ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsInsituTemperatureTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperatureTLAD::simulateObsAD(GeoVaLs & geovals,
                                             const ioda::ObsVector & ovec,
                                             const QCFlags_t & qc_flags) const {
  ufo_insitutemperature_simobs_ad_f90(keyOper_, geovals.toFortran(), obsspace(),
                                      ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsInsituTemperatureTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperatureTLAD::print(std::ostream & os) const {
  os << "ObsInsituTemperatureTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
