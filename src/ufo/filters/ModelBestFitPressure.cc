/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ModelBestFitPressure.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <vector>

#include "ioda/distribution/Accumulator.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/utils/metoffice/MetOfficeQCFlags.h"

namespace ufo {

// -----------------------------------------------------------------------------

ModelBestFitPressure::ModelBestFitPressure(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr),
    parameters_(parameters)
{
  oops::Log::trace() << "ModelBestFitPressure contructor starting" << std::endl;
  // Get parameters from options
  allvars_ += parameters_.obs_pressure;
  // Include list of required data from GeoVals
  allvars_ += Variable("air_pressure_levels@GeoVaLs");
  allvars_ += Variable("eastward_wind@GeoVaLs");
  allvars_ += Variable("northward_wind@GeoVaLs");
}

// -----------------------------------------------------------------------------

ModelBestFitPressure::~ModelBestFitPressure() {
  oops::Log::trace() << "ModelBestFitPressure destructed" << std::endl;
}

// -----------------------------------------------------------------------------
/*! \brief A filter that calculates the pressure at which the AMV wind vector (u,v)
 *  is a best match to the model wind profile. The best-fit pressure is also checked
 *  to see if it is well-constrained (True/False)
 *
 *  Example:
 *  \code{.unparsed}
 *  obs filter:
 *  - filter: Model Best Fit Pressure
 *    observation pressure:
 *      name: air_pressure_levels@MetaData
 *    top pressure: 10000
 *    pressure band half-width: 10000
 *    upper vector diff: 4
 *    lower vector diff: 2
 *    tolerance vector diff: 1.0e-8
 *    tolerance pressure: 0.01
 *    calculate bestfit winds: true
 *  \endcode
 *
 *  \author A.Martins (Met Office)
 *
 *  \date 18/05/2021: Created
 */
void ModelBestFitPressure::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "ModelBestFitPressure applyFilter" << std::endl;

  const float missing = util::missingValue(missing);
  const size_t nlocs = obsdb_.nlocs();

  // Get parameters from options.
  const float top_pressure = parameters_.top_pressure.value();
  const float pressure_band_half_width = parameters_.pressure_band_half_width.value();
  const float upper_vector_diff = parameters_.upper_vector_diff.value();
  const float lower_vector_diff = parameters_.lower_vector_diff.value();
  const float tolerance_vector_diff = parameters_.tolerance_vector_diff.value();
  const float tolerance_pressure = parameters_.tolerance_pressure.value();
  const bool calculate_best_fit_winds = parameters_.calculate_best_fit_winds.value();
  // get names of GeoVal variables
  const std::string model_pressure_name = parameters_.obs_pressure.value().variable();
  const std::string model_eastvec_name = "eastward_wind";
  const std::string model_northvec_name = "northward_wind";

  // Get GeoVaLs
  const ufo::GeoVaLs * gvals = data_.getGeoVaLs();
  // Get number of vertical levels in GeoVaLs
  const size_t num_level = data_.nlevs(Variable(model_eastvec_name + "@GeoVaLs"));

  std::vector<float> satwind_best_fit_press(nlocs, missing);
  std::vector<float> satwind_best_fit_eastward_wind;
  std::vector<float> satwind_best_fit_northward_wind;
  if (calculate_best_fit_winds) {
    satwind_best_fit_eastward_wind.resize(nlocs, missing);
    satwind_best_fit_northward_wind.resize(nlocs, missing);
  }

  // Get the observation pressure
  std::vector<float> obs_pressure(nlocs);
  data_.get(parameters_.obs_pressure, obs_pressure);

  // wind vector obs
  std::vector<float> obs_eastward(nlocs);
  obsdb_.get_db("ObsValue", model_eastvec_name, obs_eastward);
  std::vector<float> obs_northward(nlocs);
  obsdb_.get_db("ObsValue", model_northvec_name, obs_northward);

  // Get flags
  std::vector<int> u_flags(obsdb_.nlocs());
  std::vector<int> v_flags(obsdb_.nlocs());
  if (obsdb_.has("QCFlags", model_eastvec_name) &&
      obsdb_.has("QCFlags", model_northvec_name)) {
    obsdb_.get_db("QCFlags", model_eastvec_name, u_flags);
    obsdb_.get_db("QCFlags", model_northvec_name, v_flags);
  } else {
    throw eckit::Exception("eastward_wind@QCFlags or northward_wind@QCFlags not initialised");
  }

  // Vectors storing GeoVaL column for each location.
  std::vector <float> model_pressure_profile(gvals->nlevs(model_pressure_name), 0.0);
  std::vector <float> model_eastvec_profile(gvals->nlevs(model_eastvec_name), 0.0);
  std::vector <float> model_northvec_profile(gvals->nlevs(model_northvec_name), 0.0);
  std::vector <float> vec_diff(num_level);
  // diagnostic variable to be summed over all processors at the end of the routine
  std::unique_ptr<ioda::Accumulator<size_t>> countAccumulator =
      obsdb_.distribution()->createAccumulator<size_t>();

  for (size_t idata = 0; idata < nlocs; ++idata) {
    if (apply[idata]) {
      // Get GeoVaLs at the specified location.
      gvals->getAtLocation(model_pressure_profile, model_pressure_name, idata);
      gvals->getAtLocation(model_eastvec_profile, model_eastvec_name, idata);
      gvals->getAtLocation(model_northvec_profile, model_northvec_name, idata);

      // check ascending or decending (this code assume pressures are in descending order)
      // if the model pressures are ascending, reverse order
      if (std::is_sorted(model_pressure_profile.begin(),
                         model_pressure_profile.end()) ||
          model_pressure_profile.front() < model_pressure_profile.back()) {
        std::reverse(model_pressure_profile.begin(), model_pressure_profile.end());
        std::reverse(model_eastvec_profile.begin(), model_eastvec_profile.end());
        std::reverse(model_northvec_profile.begin(), model_northvec_profile.end());
      }

      // check gvals are descending
      ASSERT(!std::is_sorted(model_pressure_profile.begin(),
                             model_pressure_profile.end()));

      float min_vector_diff = std::numeric_limits<float>::max();
      const size_t UNINITIALIZED = std::numeric_limits<size_t>::max();
      size_t imin = UNINITIALIZED;

      // 1) Calculate vector difference between observed and background at all levels.
      //    Calculate best-fit pressure using vector difference.
      //    Find model level best-fit pressure (minimum vector difference).
      //    Use parabolic fit to find best-fit pressure.
      for (size_t ilev = 0; ilev < num_level - 1; ++ilev) {
        vec_diff[ilev] = std::hypot(obs_eastward[idata] - model_eastvec_profile[ilev],
                                    obs_northward[idata] - model_northvec_profile[ilev]);
        if (model_pressure_profile[ilev] < top_pressure) continue;
        if (vec_diff[ilev] < min_vector_diff) {
          min_vector_diff = vec_diff[ilev];
          imin = ilev;
        }
      }
      // check if imin set
      if (imin == UNINITIALIZED) {
        throw eckit::Exception("No model level pressure above top_pressure",
                               Here());
      }

      // use parabolic fit to find best-fit pressure
      float pressure_1;
      const float pressure_2 = model_pressure_profile[imin];
      float pressure_3;
      float vec_diff_1;
      const float vec_diff_2 = vec_diff[imin];
      float vec_diff_3;

      // if bottom model level
      if (imin == 0) {
        satwind_best_fit_press[idata] = pressure_2;
      } else {
        pressure_1 = model_pressure_profile[imin - 1];
        pressure_3 = model_pressure_profile[imin + 1];
        vec_diff_1 = vec_diff[imin - 1];
        vec_diff_3 = vec_diff[imin + 1];

        // if top of allowed region
        if (pressure_3 < top_pressure) {
          satwind_best_fit_press[idata] = pressure_2;
        // if vec_diff_2 /= vec_diff_3
        } else if (!(std::fabs(vec_diff_2 - vec_diff_3) <= tolerance_vector_diff)) {
          const float top = (((pressure_2 - pressure_1) *
                              (pressure_2 - pressure_1) *
                              (vec_diff_2 - vec_diff_3)) -
                             ((pressure_2 - pressure_3) *
                              (pressure_2 - pressure_3) *
                              (vec_diff_2 - vec_diff_1)));
          const float bottom = (((pressure_2 - pressure_1) *
                                 (vec_diff_2 - vec_diff_3)) -
                                ((pressure_2 - pressure_3) *
                                 (vec_diff_2 - vec_diff_1)));
          satwind_best_fit_press[idata] = pressure_2 -
              (0.5f * (top / bottom));
        // if vec_diff_2 = vec_diff_3 set best fit pressure = pressure_2
        } else {
          satwind_best_fit_press[idata] = pressure_2;
        }
      }

      // 2) Find bestfit eastward and northwards winds by linear interpolation,
      //    if calculate_best_fit_winds set to true (default: false)
      if (calculate_best_fit_winds) {
        size_t lev_below;
        size_t lev_above;
        float prop;
        if (std::fabs(pressure_2 - satwind_best_fit_press[idata]) <=
            tolerance_pressure) {
          satwind_best_fit_eastward_wind[idata] = model_eastvec_profile[imin];
          satwind_best_fit_northward_wind[idata] = model_northvec_profile[imin];
        } else {
          if (pressure_2 < satwind_best_fit_press[idata]) {
            lev_below = imin - 1;
            lev_above = imin;
            prop = (satwind_best_fit_press[idata] - pressure_1) /
                (pressure_2 - pressure_1);
          } else {
            lev_below = imin;
            lev_above = imin + 1;
            prop = (satwind_best_fit_press[idata] - pressure_2) /
                (pressure_3 - pressure_2);
          }
          ASSERT(prop >= 0 && prop <= 1);
          satwind_best_fit_eastward_wind[idata] = model_eastvec_profile[lev_below] *
              (1.0f - prop) + model_eastvec_profile[lev_above] * prop;
          satwind_best_fit_northward_wind[idata] = model_northvec_profile[lev_below] *
              (1.0f - prop) + model_northvec_profile[lev_above] * prop;
        }
      }

      // 3) Check if best-fit pressure is well constrained then set flag SatwindPoorConstraint.
      if (min_vector_diff <= upper_vector_diff) {
        for (size_t ilev = 0; ilev < num_level - 1; ++ilev) {
          if (model_pressure_profile[ilev] < top_pressure) continue;
          if ((model_pressure_profile[ilev] <
               satwind_best_fit_press[idata] - pressure_band_half_width ||
               model_pressure_profile[ilev] >
               satwind_best_fit_press[idata] + pressure_band_half_width) &&
              vec_diff[ilev] <= min_vector_diff + lower_vector_diff) {
            countAccumulator->addTerm(idata, 1);
            u_flags[idata] |= ufo::MetOfficeQCFlags::SatWind::SatwindPoorConstraint;
            v_flags[idata] |= ufo::MetOfficeQCFlags::SatWind::SatwindPoorConstraint;
            break;
          }
        }
      } else {
        countAccumulator->addTerm(idata, 1);
        u_flags[idata] |= ufo::MetOfficeQCFlags::SatWind::SatwindPoorConstraint;
        v_flags[idata] |= ufo::MetOfficeQCFlags::SatWind::SatwindPoorConstraint;
      }
    }  // apply
  }  // location loop

  const std::size_t iconstraint = countAccumulator->computeResult();
  if (iconstraint  > 0) {
    oops::Log::info() << "Satwind Poor constraint: "<< iconstraint
                      << " observations with modified pressure" << std::endl;
  }
  // write back flags and best-fit pressure/ winds
  obsdb_.put_db("QCFlags", model_eastvec_name, u_flags);
  obsdb_.put_db("QCFlags", model_northvec_name, v_flags);
  obsdb_.put_db("DerivedValue", "model_bestfit_pressure", satwind_best_fit_press);
  if (calculate_best_fit_winds) {
    obsdb_.put_db("DerivedValue", "model_bestfit_eastward_wind",
                  satwind_best_fit_eastward_wind);
    obsdb_.put_db("DerivedValue", "model_bestfit_northward_wind",
                  satwind_best_fit_northward_wind);
  }
}

// -----------------------------------------------------------------------------

void ModelBestFitPressure::print(std::ostream & os) const {
  os << "ModelBestFitPressure filter" << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
