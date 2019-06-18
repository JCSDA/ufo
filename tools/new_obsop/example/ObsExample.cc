/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "tools/new_obsop/example/ObsExample.h"

#include <ostream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsExample> makerExample_("Example");
// -----------------------------------------------------------------------------

ObsExample::ObsExample(const ioda::ObsSpace & odb,
                       const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), keyOper_(0), odb_(odb), varin_()
{
  int c_name_size = 800;
  char *buffin = new char[c_name_size];
  const eckit::Configuration * configc = &config;

  const oops::Variables & observed = odb.obsvariables();
  const eckit::Configuration * varconfig = &observed.toFortran();
  ufo_example_setup_f90(keyOper_, &configc, &varconfig, buffin, c_name_size);

  std::string vstr_in(buffin);
  std::vector<std::string> vvin;
  boost::split(vvin, vstr_in, boost::is_any_of("\t"));
  varin_.reset(new oops::Variables(vvin));

  oops::Log::trace() << "ObsExample created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsExample::~ObsExample() {
  ufo_example_delete_f90(keyOper_);
  oops::Log::trace() << "ObsExample destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExample::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec) const {
  ufo_example_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                         ovec.toFortran());
  oops::Log::trace() << "ObsExample: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsExample::print(std::ostream & os) const {
  os << "ObsExample::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
