/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/SatwindInversionCorrection.h"

#include <limits>

#include "ioda/distribution/Accumulator.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/utils/metoffice/MetOfficeQCFlags.h"

namespace ufo {

// -----------------------------------------------------------------------------

SatwindInversionCorrection::SatwindInversionCorrection(ioda::ObsSpace & obsdb,
                               const Parameters_ & parameters,
                               std::shared_ptr<ioda::ObsDataVector<int> > flags,
                               std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "SatwindInversion contructor starting" << std::endl;
  // Get parameters from options
  allvars_ += parameters_.obs_pressure;
  // Include list of required data from GeoVaLs
  allvars_ += Variable("GeoVaLs/air_temperature");
  allvars_ += Variable("GeoVaLs/relative_humidity");
  allvars_ += Variable("GeoVaLs/air_pressure");
}

// -----------------------------------------------------------------------------

SatwindInversionCorrection::~SatwindInversionCorrection() {
  oops::Log::trace() << "SatwindInversion destructed" << std::endl;
}

// -----------------------------------------------------------------------------
/*! \brief A filter that modifies the assigned pressure of AMV observations if a
 *  temperature inversion is detected in the model profile and defined criteria
 *  are met.
 *
 * \details The model profile is searched for the presence of a temperature
 * inversion. Where there are multiple temperature inversions, only the lowest one is found.
 * This is intended only for use on low level AMVs, typically below 700 hPa height.
 *
 * The pressure of the AMV is corrected downwards in height if the following conditions are true:
 *  a) Originally assigned pressure is greater than or equal to min_pressure (Pa).
 *  b) AMV derived from IR and visible channels only.
 *  c) Temperature inversion is present in the model profile for pressures less than or equal to
 *     max_pressure (Pa).
 *  d) In order to be considered significant, the temperature difference across the top and base of
 *     the inversion must be greater than or equal to the inversion_temperature (K) value.
 *  e) Relative humidity at the top of the inversion is less than the rh_threshold value.
 *  f) AMV has been assigned above the height of the inversion base.
 *
 * The AMV is then re-assigned to the base of the inversion.
 *
 * Reference for initial implementation:
 * Cotton, J., Forsythe, M., Warrick, F., (2016). Towards improved height assignment and
 * quality control of AMVs in Met Office NWP. Proceedings for the 13th International Winds
 * Workshop 27 June - 1 July 2016, Monterey, California, USA.
 *
 * Example:
 * \code{.unparsed}
 *  obs filters:
 *  - filter: Satwind Inversion Correction
 *    observation pressure:
 *      name: MetaData/pressure
 *    RH threshold: 50
 *    maximum pressure: 96000
 * \endcode
 *
 * \author J.Cotton (Met Office)
 *
 * \date 07/05/2021: Created
 */

void SatwindInversionCorrection::applyFilter(const std::vector<bool> & apply,
                                             const Variables & filtervars,
                                             std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "SatwindInversionCorrection priorFilter" << std::endl;
  print(oops::Log::trace());

  const float missing = util::missingValue<float>();
  const size_t nlocs = obsdb_.nlocs();

// Get parameters from options.
  const float rh_threshold = parameters_.rh_threshold.value();
  const float min_pressure = parameters_.min_pressure.value();
  const float max_pressure = parameters_.max_pressure.value();
  const float inversion_temperature = parameters_.inversion_temperature.value();
// get names of GeoVal variables
  const oops::Variable model_temp_name{"air_temperature"};
  const oops::Variable model_rh_name{"relative_humidity"};
  const oops::Variable model_vcoord_name{"air_pressure"};

// Get variables from ObsSpace
// Get the observation pressure
  std::vector<float> obs_pressure;
  data_.get(parameters_.obs_pressure, obs_pressure);
// Get wind computation method
  std::vector<int> comp_method(obsdb_.nlocs());
  obsdb_.get_db("MetaData", "windComputationMethod", comp_method);
// Get flags
  std::vector<int> u_flags(obsdb_.nlocs());
  std::vector<int> v_flags(obsdb_.nlocs());
  if (obsdb_.has("QCFlags", "windEastward") && obsdb_.has("QCFlags", "windNorthward")) {
    obsdb_.get_db("QCFlags", "windEastward", u_flags);
    obsdb_.get_db("QCFlags", "windNorthward", v_flags);
  } else {
    throw eckit::Exception("QCFlags/windEastward or QCFlags/windNorthward not initialised",
                           Here());
  }
// Get GeoVaLs
  const ufo::GeoVaLs * gvals = data_.getGeoVaLs();
// Get number of vertical levels in GeoVaLs
  const size_t nlevs = gvals->nlevs(model_vcoord_name);
// Vectors storing GeoVaL column for current location.
  std::vector <double> model_temp_profile(gvals->nlevs(model_temp_name), 0.0);
  std::vector <double> model_rh_profile(gvals->nlevs(model_rh_name), 0.0);
  std::vector <double> model_vcoord_profile(gvals->nlevs(model_vcoord_name), 0.0);

// Vector to store original pressure
  std::vector<float> original_pressure(obs_pressure);

// diagnostic variables to be summed over all processors at the end of the routine
  std::unique_ptr<ioda::Accumulator<size_t>> countAccumulator =
      obsdb_.distribution()->createAccumulator<size_t>();
  std::unique_ptr<ioda::Accumulator<double>> pdiffAccumulator =
      obsdb_.distribution()->createAccumulator<double>();

// Loop through locations
  for (size_t iloc=0; iloc < nlocs; ++iloc) {
    if (apply[iloc]) {
      //  only consider low level AMVs from infrared or visible channels
      if (obs_pressure[iloc] >= min_pressure &&
          (comp_method[iloc] == CloudMotionMethod::infrared ||
           comp_method[iloc] == CloudMotionMethod::visible)) {
        // Get GeoVaLs at the specified location.
        gvals->getAtLocation(model_temp_profile, model_temp_name, iloc);
        gvals->getAtLocation(model_rh_profile, model_rh_name, iloc);
        gvals->getAtLocation(model_vcoord_profile, model_vcoord_name, iloc);
        // Check GeoVaLs are in correct vertical order
        if (model_vcoord_profile.front() > model_vcoord_profile.back()) {
          throw eckit::BadValue("GeoVaLs are not ordered from model top to bottom", Here());
        }
        // ---------------------------------------------------------------------------
        //  Search for inversion and if present find T and P of base and top
        // ---------------------------------------------------------------------------
        bool inversion = false;
        bool firsttime = true;
        float inversion_base = std::numeric_limits<float>::max();
        float inversion_top = std::numeric_limits<float>::max();
        float temp_inversion_base = std::numeric_limits<float>::max();
        float temp_inversion_top = std::numeric_limits<float>::max();
        //  loop over levels starting from highest pressure (bottom to top)
        for (int ilev  = nlevs-1; ilev >= 1; ilev--) {
          //  if haven't found inversion and pressure is less than min_pressure Pa then exit
          if (inversion == false && model_vcoord_profile[ilev] < min_pressure) {
            break;
          }
          // if model pressure is greater than max_pressure Pa then skip to next level
          if (model_vcoord_profile[ilev] > max_pressure) {
            continue;
          }
          //  if haven't found inversion, check for increase in temperature
          if (!inversion && model_temp_profile[ilev-1] > model_temp_profile[ilev]) {
            //  T of level above is greater so take base at current level
            inversion_base =  model_vcoord_profile[ilev];
            temp_inversion_base = model_temp_profile[ilev];
            inversion = true;
          }
          // if inversion found, then detect level the temperature starts to decrease again above
          // the inversion base
          if (inversion &&
              model_temp_profile[ilev-1] < model_temp_profile[ilev] &&
              model_vcoord_profile[ilev] < inversion_base &&
              firsttime) {
            //  Check humidity of inversion top
            if (model_rh_profile[ilev] < rh_threshold) {
              inversion_top = model_vcoord_profile[ilev];
              temp_inversion_top = model_temp_profile[ilev];
              firsttime = false;
            } else {
              // otherwise discard and keep looking for a new base/top higher up
              inversion = false;
            }
          }
        }  // level loop
        //
        //---------------------------------------------------------------------------
        // Correct pressure if inversion found and conditions are met
        //---------------------------------------------------------------------------
        if (inversion &&
            (temp_inversion_top - temp_inversion_base) >= inversion_temperature &&
            obs_pressure[iloc] < inversion_base) {
          // re-assign to base of inversion
          countAccumulator->addTerm(iloc, 1);
          pdiffAccumulator->addTerm(iloc, inversion_base - obs_pressure[iloc]);
          obs_pressure[iloc] = inversion_base;
          // set flag
          u_flags[iloc] |= ufo::MetOfficeQCFlags::SatWind::SatwindInversionFlag;
          v_flags[iloc] |= ufo::MetOfficeQCFlags::SatWind::SatwindInversionFlag;
        }
      }
    }  // apply
  }  // location loop
  //  write back corrected pressure, updated flags and original pressure
  obsdb_.put_db(parameters_.obs_pressure.value().group(),
                parameters_.obs_pressure.value().variable(), obs_pressure);
  obsdb_.put_db("QCFlags", "windEastward", u_flags);
  obsdb_.put_db("QCFlags", "windNorthward", v_flags);
  obsdb_.put_db(parameters_.obs_pressure.value().group(),
                parameters_.obs_pressure.value().variable() + std::string("_original"),
                original_pressure);

  // sum number corrected and pressure differences
  const std::size_t count = countAccumulator->computeResult();
  const double pdiff = pdiffAccumulator->computeResult();
  if (count) {
    oops::Log::info() << "Satwind Inversion: "<< count
                                       << " observations with modified pressure" << std::endl;
    oops::Log::info() << "Satwind Inversion: "<< pdiff / count
                                       << " Pa mean pressure difference" << std::endl;
  }
}

// -----------------------------------------------------------------------------

void SatwindInversionCorrection::print(std::ostream & os) const {
  os << "SatwindInversionCorrection: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
