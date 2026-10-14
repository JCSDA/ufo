/*
 * (C) Copyright 2017-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/crtm/ObsRadianceCRTMTLAD.h"

#include <set>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/util/dateFunctions.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/TimeWindow.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/crtm/ObsRadianceCRTMTLAD.interface.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsRadianceCRTMTLAD> makerCRTMTL_("CRTM");
// -----------------------------------------------------------------------------

void ObsRadianceCRTMTLAD::nullifyJacParamSetup(
                        const CRTMLinearObsOperatorParameters & linearParams,
                        eckit::LocalConfiguration & conf) {
  if (linearParams.NullifyJacobianParams.value() == boost::none) {
    oops::Log::info() << "No nullifyJacobianParams specified in the config. "
                    << "Jacobian will not be nullified for any variable or channel."
                    << std::endl;
    return;
  }
  const std::vector<NullifyJacobianParameters> & nullifyJacobianParams =
        linearParams.NullifyJacobianParams.value().value();

  oops::Log::info() << "Nullifying Jacobian for size:" << nullifyJacobianParams.size()
                  << std::endl;
  std::vector<eckit::LocalConfiguration> nullifyJacobianConfs(nullifyJacobianParams.size());
  for (size_t jj = 0; jj < nullifyJacobianParams.size(); ++jj) {
    const std::set<int> channels =
                        oops::parseIntSet(nullifyJacobianParams[jj].channelList.value());
    nullifyJacobianConfs[jj].set("variable", nullifyJacobianParams[jj].variable.value());
    nullifyJacobianConfs[jj].set("channels", std::vector<int>(channels.begin(), channels.end()));
  }

  eckit::LocalConfiguration linearConf = conf.getSubConfiguration("linear obs operator");
  linearConf.set("nullify jacobian params", nullifyJacobianConfs);
  conf.set("linear obs operator", linearConf);
}

// -----------------------------------------------------------------------------

ObsRadianceCRTMTLAD::ObsRadianceCRTMTLAD(const ioda::ObsSpace & odb,
                                         const Parameters_ & params)
  : LinearObsOperatorBase(odb), keyOperRadianceCRTM_(0), varin_()
{
  // parse channels from the config and create variable names
  const oops::ObsVariables & observed = odb.assimvariables();
  std::vector<int> channels_list = observed.channels();

  // get a single central of middle time from observation space
  const util::DateTime midPoint = odb.timeWindow().midpoint();
  std::string year, month, day, hour, minute, second;
  midPoint.toYYYYMMDDhhmmss(year, month, day, hour,  minute,  second);
  // Julian Day Number since noon Universal Time (UT) on January 1, 4713 BCE
  uint64_t midPointJulday = util::datefunctions::dateToJulian(std::stoi(year),
                                                              std::stoi(month),
                                                              std::stoi(day));

  eckit::LocalConfiguration conf = params.toConfiguration();
  if (params.LinearObsOperator.value() != boost::none) {
    nullifyJacParamSetup(params.LinearObsOperator.value().value(), conf);
  }

  // call Fortran setup routine
  ufo_radiancecrtm_tlad_setup_f90(keyOperRadianceCRTM_, conf,
                                  channels_list.size(), channels_list[0], midPointJulday,
                                  varin_, odb.comm());

  oops::Log::trace() << "ObsRadianceCRTMTLAD constructor done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceCRTMTLAD::~ObsRadianceCRTMTLAD() {
  ufo_radiancecrtm_tlad_delete_f90(keyOperRadianceCRTM_);
  oops::Log::trace() << "ObsRadianceCRTMTLAD destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTMTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics & ydiags,
                                        const QCFlags_t & qc_flags) {
  ufo_radiancecrtm_tlad_settraj_f90(keyOperRadianceCRTM_, geovals.toFortran(), obsspace(),
                                    ydiags.toFortran(), qc_flags);
  oops::Log::trace() << "ObsRadianceCRTMTLAD::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTMTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_radiancecrtm_simobs_tl_f90(keyOperRadianceCRTM_, geovals.toFortran(), obsspace(),
                                 ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsRadianceCRTMTLAD::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTMTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_radiancecrtm_simobs_ad_f90(keyOperRadianceCRTM_, geovals.toFortran(), obsspace(),
                                 ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsRadianceCRTMTLAD::simulateObsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTMTLAD::print(std::ostream & os) const {
  os << "ObsRadianceCRTMTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
