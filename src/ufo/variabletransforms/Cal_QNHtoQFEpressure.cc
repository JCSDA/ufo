/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/Cal_QNHtoQFEpressure.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/metoffice/MetOfficeQCFlags.h"


namespace ufo {

/************************************************************************************/
//  Cal_QNHtoQFEpressure
/************************************************************************************/

static TransformMaker<Cal_QNHtoQFEpressure>
    makerCal_QNHtoQFEpressure_("QNHtoQFEpressure");

Cal_QNHtoQFEpressure::Cal_QNHtoQFEpressure(
        const GenericVariableTransformParameters &options,
        const ObsFilterData &data,
        const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
        const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
      : TransformBase(options, data, flags, obserr) {}

/************************************************************************************/

void Cal_QNHtoQFEpressure::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " --> Retrieve QFE (pressure reduced to official aerodrome altitude)"
            << std::endl;
  oops::Log::trace() << "      --> method: " << method() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  const size_t nlocs = obsdb_.nlocs();
  const float missing = util::missingValue(missing);

  // Get all required obs, metadata and geovals
  // Obs
  std::vector<float> PStn, Pmsl;
  getObservation("ObsValue", "station_pressure",
                 PStn, true);
  getObservation("ObsValue", "mean_sea_level_pressure",
                 Pmsl, true);

  // MetaData
  std::vector<float> ZStn, PStn_error;
  getObservation("MetaData", "station_elevation",
                 ZStn, true);
  getObservation("ObsError", "station_pressure",
                 PStn_error, true);
  std::vector<int> PStn_flag;
  getObservation("QCFlags", "station_pressure",
                 PStn_flag, true);

  if (!oops::allVectorsSameNonZeroSize(Pmsl, PStn)) {
    oops::Log::warning() << "Vector sizes: " << oops::listOfVectorSizes(Pmsl, PStn)
                         << std::endl;
    throw eckit::BadValue("At least one vector is the wrong size or empty out of "
                          "Station pressure and PMSL", Here());
  }

  const float RoundErr = 2500.0/3.0;  // Error variance caused by rounding to nearest 100 Pa
  // The constants below are used in the calculation of QFE; the stated values are taken from
  // the ICAO (2011) Manual on Automatic Meterological Observing Systems at Aerodromes Doc 9837
  const float ExpQFE1 = Constants::icao_lapse_rate_l*Constants::rd_over_g;  // 0.190263;
  const float ExpQFE2 = Constants::g_over_rd/Constants::icao_lapse_rate_l;  // 5.25588
  const float C1 = Constants::icao_temp_surface/Constants::icao_lapse_rate_l;  // 44330.77
  const float C2 = C1/std::pow(Constants::icao_pressure_surface, ExpQFE1);  // 11880.32

  // Loop over all obs to calculate QFE
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (!apply[jj]) continue;
    if (ZStn[jj] != missing && Pmsl[jj] > 0.0) {
      float Palt = Pmsl[jj];  // Altimeter pressure (QNH)
      if ((!(PStn_flag[jj] & MetOfficeQCFlags::Surface::notRoundedFlag))
           && (std::fmod(Palt, 100.0) == 0.0)) {
        Palt += 50.0;  // To avoid systematic bias add 0.5hPa to values that have been rounded down.
        PStn_error[jj] = std::sqrt(std::pow(PStn_error[jj], 2) + RoundErr);
        PStn_flag[jj] |= ufo::MetOfficeQCFlags::Surface::QNHhPaFlag;
      } else {
        PStn_flag[jj] |= ufo::MetOfficeQCFlags::Surface::QNHinHgFlag;
      }
      const float H = Constants::icao_temp_surface
                *(1.0 - std::pow(Palt/(100.0*Constants::icao_pressure_surface), ExpQFE1))
                / Constants::icao_lapse_rate_l + ZStn[jj];
      PStn[jj] = 100.0 * std::pow((C1 - H)/C2, ExpQFE2);
    }
  }

  putObservation("station_pressure", PStn);
  obsdb_.put_db("DerivedQCFlags", "station_pressure", PStn_flag);
  obsdb_.put_db("DerivedObsError", "station_pressure", PStn_error);
}
}  // namespace ufo

