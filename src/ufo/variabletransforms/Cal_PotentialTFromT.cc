/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/Cal_PotentialTFromT.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/metoffice/MetOfficeQCFlags.h"


namespace ufo {

/************************************************************************************/
//  Cal_PotentialTFromT
/************************************************************************************/

static TransformMaker<Cal_PotentialTFromT>
    makerCal_PotentialTFromT_("PotentialTFromT");

Cal_PotentialTFromT::Cal_PotentialTFromT(
        const Parameters_ &options,
        const ObsFilterData &data,
        const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
        const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr),
      pressurevariable_(options.PressureVariable),
      pressuregroup_(options.PressureGroup),
      temperaturevariable_(options.TemperatureVariable),
      potentialtempvariable_(options.PotentialTempVariable)
{}

/************************************************************************************/
/*
 The potential temperature is the temperature
 an unsaturated parcel would have if lowered (or raised) to a level of 100 kPa.

Method: -
 potential_temperature = (reference_pressure/pressure)^kappa  * temperature

 where kappa=gas constant dry air / specific heat capacity dry air
*/
void Cal_PotentialTFromT::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " Retrieve Potential Temperature from Temperature" << std::endl;

  SetUseValidDataOnly(true);

  const size_t nlocs = obsdb_.nlocs();
  const float kappa = Constants::rd_over_cp;

  // Get all required data
  std::vector<float> pressure;
  std::vector<float> temperature;
  std::vector<float> tempError;
  std::vector<float> potTemp(nlocs, missingValueFloat);
  std::vector<float> potTempErr(nlocs, missingValueFloat);
  std::vector<float> temppge;
  std::vector<int> tempflags;
  getObservation(pressuregroup_, pressurevariable_,
                 pressure, true);
  getObservation("ObsValue", temperaturevariable_,
                 temperature, true);
  data_.get(Variable(std::string("ObsErrorData/") + temperaturevariable_), tempError);
  getObservation("GrossErrorProbability", temperaturevariable_,
                 temppge);
  getObservation("QCFlags", temperaturevariable_,
                 tempflags);

  if (temppge.empty()) {
    temppge.assign(nlocs, missingValueFloat);
  }
  if (tempflags.empty()) {
    tempflags.assign(nlocs, 0);
  }

  if (!oops::allVectorsSameNonZeroSize(pressure, temperature, tempError)) {
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(pressure,
                                                    temperature,
                                                    tempError)
                         << std::endl;
    throw eckit::BadValue("At least one vector is the wrong size or empty out of "
                          "pressure, temperature and temperature error", Here());
  }

  for (ioda::ObsSpace::RecIdxIter irec = obsdb_.recidx_begin();
       irec != obsdb_.recidx_end(); ++irec) {
    const std::vector<std::size_t> &rSort = obsdb_.recidx_vector(irec);

    // Loop over each record
    for (size_t iloc : rSort) {
      if (!apply[iloc]) continue;
      const float temp = temperature[iloc];
      const float pres = pressure[iloc];
      const float terr = tempError[iloc];

      if (pres == missingValueFloat ||
          pres <= 0.0f ||
          temp == missingValueFloat ||
          terr == missingValueFloat) continue;

      const float conversion = std::pow((Constants::pref / pres), kappa);
      potTemp[iloc] = conversion * temp;
      potTempErr[iloc] = conversion * terr;
    }
  }

  obsdb_.put_db("DerivedObsValue", potentialtempvariable_, potTemp);
  const size_t iv = obserr_.varnames().find(potentialtempvariable_);
  for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
    if (!apply[jobs])
      continue;
    obserr_[iv][jobs] = potTempErr[jobs];
  }

  // copy temperature's PGEFinal and QCflags to new potentialTemperature
  obsdb_.put_db("GrossErrorProbability", potentialtempvariable_, temppge);
  obsdb_.put_db("QCFlags", potentialtempvariable_, tempflags);
}
}  // namespace ufo
