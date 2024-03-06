/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <algorithm>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/gnssro/BendMetOffice/ObsGnssroBendMetOffice.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGnssroBendMetOffice> makerGnssroBendMetOffice_("GnssroBendMetOffice");
// -----------------------------------------------------------------------------

ObsGnssroBendMetOffice::ObsGnssroBendMetOffice(const ioda::ObsSpace & odb,
                                       const Parameters_ & parameters)
  : ObsOperatorBase(odb), keyOperGnssroBendMetOffice_(0), odb_(odb), varin_()
{
  oops::Log::trace() << "Constructing obs operator" << std::endl;
  oops::Log::debug() << "Number of channels " << odb.assimvariables().channels() << std::endl;

  const std::vector<std::string> vv{"air_pressure_levels", "specific_humidity",
                                    "geopotential_height", "geopotential_height_levels"};
  varin_.reset(new oops::Variables(vv));

  oops::Log::debug() << "nlocs " << odb.nlocs() << std::endl;
  oops::Log::debug() << "nchans " << odb.nchans() << std::endl;
  oops::Log::debug() << "nrecs " << odb.nrecs() << std::endl;
  oops::Log::debug() << "nvars " << odb.nvars() << std::endl;

  std::set<int> channelset = oops::parseIntSet(parameters.channelList);
  std::vector<int> channels;
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels));

  ufo_gnssro_bendmetoffice_setup_f90(keyOperGnssroBendMetOffice_,
                                     parameters.vertInterpOPS,
                                     parameters.pseudoLevels,
                                     parameters.minTempGrad,
                                     channels.size(),
                                     channels[0],
                                     parameters.noSuperCheck);

  oops::Log::trace() << "ObsGnssroBendMetOffice created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroBendMetOffice::~ObsGnssroBendMetOffice() {
  ufo_gnssro_bendmetoffice_delete_f90(keyOperGnssroBendMetOffice_);
  oops::Log::trace() << "ObsGnssroBendMetOffice destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBendMetOffice::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                       ObsDiagnostics & ydiags, const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "Starting simulateObs" << std::endl;
  oops::Log::debug() << "ObsVector: nvars = " << ovec.nvars() << "  nlocs = "
                     << ovec.nlocs() << "  size = " << ovec.size() << std::endl;

  ufo_gnssro_bendmetoffice_simobs_f90(keyOperGnssroBendMetOffice_, gom.toFortran(), odb_,
                                  ovec.nvars(), ovec.nlocs(), ovec.toFortran(), ydiags.toFortran());
  oops::Log::trace() << "Finishing simulateObs" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroBendMetOffice::print(std::ostream & os) const {
  os << "ObsGnssroBendMetOffice::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
