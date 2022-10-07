/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/OceanConversions/OceanDensity.h"

#include <cmath>

#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static TransformMaker<OceanDensity>
    makerOceanDensity_("OceanDensity");

OceanDensity::OceanDensity(
        const Parameters_ &options,
        const ObsFilterData &data,
        const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
        const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr),
      salinityvariable_(options.SalinityVariable),
      salinitygroup_(options.SalinityGroup),
      temperaturevariable_(options.TemperatureVariable),
      temperaturegroup_(options.TemperatureGroup),
      pressurevariable_(options.PressureVariable),
      pressuregroup_(options.PressureGroup),
      densityvariable_(options.DensityVariable)
{}

// -----------------------------------------------------------------------------

void OceanDensity::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << "Retrieve Ocean Density from Salinity, Temperature, Pressure" << std::endl;

  // dimension
  const size_t nlocs = obsdb_.nlocs();

  // Get all required data
  std::vector<float> pressure;
  std::vector<float> temp;
  std::vector<float> sal;
  std::vector<float> density(nlocs, missingValueFloat);
  std::vector<float> densityError(nlocs, 1.0);
  std::vector<float> densitypge;
  std::vector<int> densityflags;
  getObservation(salinitygroup_,
                 salinityvariable_,
                 sal, true);
  getObservation(temperaturegroup_,
                 temperaturevariable_,
                 temp, true);
  getObservation(pressuregroup_,
                 pressurevariable_,
                 pressure, true);
  getObservation("GrossErrorProbability", densityvariable_,
                 densitypge);
  getObservation("QCFlags", densityvariable_,
                 densityflags);

  if (densitypge.empty()) {
    densitypge.assign(nlocs, missingValueFloat);
  }
  if (densityflags.empty()) {
    densityflags.assign(nlocs, 0);
  }

  for (size_t loc = 0; loc < nlocs; ++loc) {
    if (!apply[loc]) continue;
    if (sal[loc] != missingValueFloat &&
        temp[loc] != missingValueFloat &&
        pressure[loc] != missingValueFloat) {
      density[loc] = gsw_rho_t_exact_f90(sal[loc],
                                           temp[loc],
                                           pressure[loc]);
    }
  }
  obsdb_.put_db("DerivedObsValue", densityvariable_, density);
  const size_t iv = obserr_.varnames().find(densityvariable_);
  for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
    if (!apply[jobs])
      continue;
    obserr_[iv][jobs] = densityError[jobs];
  }

  // copy PGEFinal and QCflags to new density
  obsdb_.put_db("GrossErrorProbability", densityvariable_, densitypge);
  obsdb_.put_db("QCFlags", densityvariable_, densityflags);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
