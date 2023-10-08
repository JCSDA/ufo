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
  const float missing = util::missingValue<float>();

  // Get all required obs, metadata and geovals
  // Obs
  std::vector<float> PStn, Pmsl;
  getObservation("ObsValue", "stationPressure",
                 PStn);
  getObservation("ObsValue", "pressureReducedToMeanSeaLevel",
                 Pmsl);

  // MetaData
  std::vector<float> ZStn, PStn_error;
  if (data_.has(Variable("MetaData/correctedStationAltitude"))) {
    getObservation("MetaData", "correctedStationAltitude",
                     ZStn);
  } else {
    getObservation("MetaData", "stationElevation",
                     ZStn);
  }
  // Error statistics
  data_.get(Variable("ObsErrorData/stationPressure"), PStn_error);

  // Flags
  std::vector<bool> notRounded_flag;
  std::vector<bool> QNHinHg_flag;
  std::vector<bool> QNHhPa_flag;

  // get diagnostic flags from ObsSpace (warn if they have not yet been created)
  if (obsdb_.has("DiagnosticFlags/notRounded", "stationPressure")) {
    obsdb_.get_db("DiagnosticFlags/notRounded", "stationPressure", notRounded_flag);
  } else {
    throw eckit::UserError("Variable 'DiagnosticFlags/notRounded/stationPressure' does not "
                           "exist yet. It needs to be set up with the 'Create Diagnostic "
                           "Flags' filter prior to using the 'set' or 'unset' action.");
  }
  if (obsdb_.has("DiagnosticFlags/QNHinHg", "stationPressure")) {
    obsdb_.get_db("DiagnosticFlags/QNHinHg", "stationPressure", QNHinHg_flag);
  } else {
    throw eckit::UserError("Variable 'DiagnosticFlags/QNHinHg/stationPressure' does not exist yet. "
                           "It needs to be set up with the 'Create Diagnostic Flags' filter "
                           "prior to using the 'set' or 'unset' action.");
  }
  if (obsdb_.has("DiagnosticFlags/QNHhPa", "stationPressure")) {
    obsdb_.get_db("DiagnosticFlags/QNHhPa", "stationPressure", QNHhPa_flag);
  } else {
    throw eckit::UserError("Variable 'DiagnosticFlags/QNHhPa/stationPressure' does not exist yet. "
                           "It needs to be set up with the 'Create Diagnostic Flags' filter "
                           "prior to using the 'set' or 'unset' action.");
  }

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
      if ((!notRounded_flag[jj]) && (std::fmod(Palt, 100.0) == 0.0)) {
        Palt += 50.0;  // To avoid systematic bias add 0.5hPa to values that have been rounded down.
        PStn_error[jj] = std::sqrt(std::pow(PStn_error[jj], 2) + RoundErr);
        QNHhPa_flag[jj] = true;
      } else {
        QNHinHg_flag[jj] = true;
      }
      const float H = Constants::icao_temp_surface
                *(1.0 - std::pow(Palt/(100.0*Constants::icao_pressure_surface), ExpQFE1))
                / Constants::icao_lapse_rate_l + ZStn[jj];
      PStn[jj] = 100.0 * std::pow((C1 - H)/C2, ExpQFE2);
    }
  }

  putObservation("stationPressure", PStn);
  obsdb_.put_db("DiagnosticFlags/QNHhPa", "stationPressure", QNHhPa_flag);
  obsdb_.put_db("DiagnosticFlags/QNHinHg", "stationPressure", QNHinHg_flag);
  obsdb_.put_db("DerivedObsError", "stationPressure", PStn_error);
  const size_t iv = obserr_.varnames().find("stationPressure");
  for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
    if (!apply[jobs])
      continue;
    obserr_[iv][jobs] = PStn_error[jobs];
  }
}
}  // namespace ufo

