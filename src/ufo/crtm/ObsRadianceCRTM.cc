/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/crtm/ObsRadianceCRTM.h"

#include <algorithm>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadianceCRTM> makerCRTM_("CRTM");

// -----------------------------------------------------------------------------

ObsRadianceCRTM::ObsRadianceCRTM(const ioda::ObsSpace & odb,
                                 const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOperRadianceCRTM_(0),
    odb_(odb), varin_()
{
  int c_name_size = 3000;
  char *buffin = new char[c_name_size];

  // parse channels from the config and create variable names
  const oops::Variables & observed = odb.obsvariables();
  std::vector<int> channels_list = observed.channels();

  // call Fortran setup routine
  const eckit::Configuration * configc = &config;
  ufo_radiancecrtm_setup_f90(keyOperRadianceCRTM_, &configc, channels_list.size(),
                             channels_list[0], buffin, c_name_size);

  std::string vstr_in(buffin);
  std::vector<std::string> vvin;
  boost::split(vvin, vstr_in, boost::is_any_of("\t"));
  varin_.reset(new oops::Variables(vvin));

  oops::Log::info() << "ObsRadianceCRTM channels: " << channels_list << std::endl;
  oops::Log::trace() << "ObsRadianceCRTM created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceCRTM::~ObsRadianceCRTM() {
  ufo_radiancecrtm_delete_f90(keyOperRadianceCRTM_);
  oops::Log::trace() << "ObsRadianceCRTM destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTM::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec) const {
  ufo_radiancecrtm_simobs_f90(keyOperRadianceCRTM_, gom.toFortran(), odb_,
                          ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsRadianceCRTM simulateObs done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTM::print(std::ostream & os) const {
  os << "ObsRadianceCRTM::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
