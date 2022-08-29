/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <string>

#include "oops/util/missingValues.h"
#include "ufo/utils/metoffice/MetOfficeQCFlags.h"
#include "ufo/variabletransforms/Cal_PStar.h"

namespace ufo {

/************************************************************************************/
//  Cal_PStar
/************************************************************************************/

static TransformMaker<Cal_PStar>
    makerCal_PStar_("PStar");

Cal_PStar::Cal_PStar(
    const GenericVariableTransformParameters &options,
    const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr), gvals_() {
    gvals_ += Variable("surf_param_a@GeoVaLs");
    gvals_ += Variable("surf_param_b@GeoVaLs");
    gvals_ += Variable("surface_altitude@GeoVaLs");
    gvals_ += Variable("surface_pressure@GeoVaLs");
}

/************************************************************************************/

void Cal_PStar::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " --> Calculate pressure at model surface"
            << std::endl;
  oops::Log::trace() << "      --> method: " << method() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  const size_t nlocs = obsdb_.nlocs();
  const float missing = util::missingValue(missing);

  // Get all required obs, metadata and geovals
  // Obs
  std::vector<float> PStn, PStd, Pmsl, ZStn, ZStd;
  getObservation("ObsValue", "stationPressure",
                 PStn, true);
  getObservation("ObsValue", "standardPressure",
                 PStd, true);
  getObservation("ObsValue", "pressureReducedToMeanSeaLevel",
                 Pmsl, true);

  // MetaData
  if (data_.has(Variable("MetaData/correctedStationAltitude"))) {
    getObservation("MetaData", "correctedStationAltitude",
                     ZStn);
  } else {
    getObservation("MetaData", "stationAltitude",
                     ZStn);
  }
  getObservation("MetaData", "standardAltitude",
                 ZStd, true);

  // Errors
  std::vector<float> PStn_error, PStd_error, Pmsl_error;
  getObservation("ObsError", "stationPressure",
                 PStn_error, true);
  getObservation("ObsError", "standardPressure",
                 PStd_error, true);
  getObservation("ObsError", "pressureReducedToMeanSeaLevel",
                 Pmsl_error, true);

  // PGEs
  std::vector<float> PStn_PGE, PStd_PGE, Pmsl_PGE;
  getObservation("GrossErrorProbability", "stationPressure",
                 PStn_PGE, true);
  getObservation("GrossErrorProbability", "standardPressure",
                 PStd_PGE, true);
  getObservation("GrossErrorProbability", "pressureReducedToMeanSeaLevel",
                 Pmsl_PGE, true);

  // Geovals
  std::vector<float> PSurfParamA, PSurfParamB, ModelHeight, BkPStar;
  data_.get(Variable("surf_param_a@GeoVaLs"), PSurfParamA);
  data_.get(Variable("surf_param_b@GeoVaLs"), PSurfParamB);
  data_.get(Variable("surface_altitude@GeoVaLs"), ModelHeight);
  data_.get(Variable("surface_pressure@GeoVaLs"), BkPStar);

  if (!oops::allVectorsSameNonZeroSize(PStn, PStd, Pmsl)) {
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(PStn, PStd, Pmsl)
                         << std::endl;
    throw eckit::BadValue("At least one vector is the wrong size or empty out of "
                          "Station pressure, Standard pressure and PMSL", Here());
  }

  // Assign vectors to write obs/flags/errors/PGEs to
  std::vector<float> PStar(nlocs), PStar_error(nlocs), PStar_PGE;
  PStar.assign(nlocs, missing);
  PStar_error.assign(nlocs, missing);
  PStar_PGE.assign(nlocs, missing);

  std::vector<int> PStar_flag(nlocs);
  PStar_flag.assign(nlocs, 0);

  std::vector<float> BkPStn(nlocs), BkPStd(nlocs), BkPmsl(nlocs);
  BkPStn.assign(nlocs, missing);
  BkPStd.assign(nlocs, missing);
  BkPmsl.assign(nlocs, missing);

  std::vector<bool> PreferredVariable_flag;
  std::vector<bool> PmslUsed_flag(nlocs, false);
  std::vector<bool> PstdUsed_flag(nlocs, false);
  std::vector<bool> PstnUsed_flag(nlocs, false);

  // get diagnostic flags from ObsSpace (warn if they have not yet been created)
  if (obsdb_.has("DiagnosticFlags/PreferredVariable", "stationPressure")) {
    obsdb_.get_db("DiagnosticFlags/PreferredVariable", "stationPressure", PreferredVariable_flag);
  } else {
    throw eckit::UserError("Variable 'DiagnosticFlags/PreferredVariable/stationPressure' does not "
                           "exist yet. It needs to be set up with the 'Create Diagnostic "
                           "Flags' filter prior to using the 'set' or 'unset' action.");
  }

  // Loop over all obs to calculate PStar
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (!apply[jj]) continue;

    // Set indicators 0 = missing, 1 = rejected, 2 = present
    int IvPStn = 0, IvPStd = 0, IvPmsl = 0;
    if (ZStn[jj] != missing && PStn[jj] != missing) {
      IvPStn = 2;
      if ((flags_)["stationPressure"][jj] != QCflags::pass) {IvPStn = 1;}
      BkPStn[jj] = formulas::BackgroundPressure(PSurfParamA[jj], PSurfParamB[jj], ZStn[jj]);
    }

    if (ZStd[jj] != missing && PStd[jj] != missing) {
      IvPStd = 2;
      if ((flags_)["standardPressure"][jj] != QCflags::pass) {IvPStd = 1;}
      BkPStd[jj] = formulas::BackgroundPressure(PSurfParamA[jj], PSurfParamB[jj], ZStd[jj]);
    }

    if (Pmsl[jj] != missing) {
      IvPmsl = 2;
      if ((flags_)["pressureReducedToMeanSeaLevel"][jj] != QCflags::pass) {IvPmsl = 1;}
      BkPmsl[jj] = formulas::BackgroundPressure(PSurfParamA[jj], PSurfParamB[jj], 0.0);
    }

    int IvPmax = std::max(IvPStn, std::max(IvPStd, IvPmsl));

    if (IvPmax == 0) continue;
    if (IvPStn == IvPmax && BkPStn[jj] != missing && PreferredVariable_flag[jj]) {
      PStar[jj] = (BkPStar[jj]*PStn[jj])/BkPStn[jj];
      PstnUsed_flag[jj] = true;
      PStar_error[jj] = PStn_error[jj];
      PStar_PGE[jj] = PStn_PGE[jj];
    } else if (BkPStd[jj] != missing && BkPmsl[jj] == missing) {
      PStar[jj] = (BkPStar[jj]*PStd[jj])/BkPStd[jj];
      PstdUsed_flag[jj] = true;
      PStar_error[jj] = PStd_error[jj];
      PStar_PGE[jj] = PStd_PGE[jj];
    } else if (BkPmsl[jj] != missing) {
      PStar[jj] = BkPStar[jj]*Pmsl[jj]/BkPmsl[jj];
      PmslUsed_flag[jj] = true;
      PStar_error[jj] = Pmsl_error[jj];
      PStar_PGE[jj] = Pmsl_PGE[jj];
    } else {
        continue;
    }
  }

  putObservation("surface_pressure", PStar);
  obsdb_.put_db("DerivedObsError", "surface_pressure", PStar_error);
  obsdb_.put_db("GrossErrorProbability", "surface_pressure", PStar_PGE);
  obsdb_.put_db("DiagnosticFlags/PmslUsed", "surface_pressure", PmslUsed_flag);
  obsdb_.put_db("DiagnosticFlags/PstdUsed", "surface_pressure", PstdUsed_flag);
  obsdb_.put_db("DiagnosticFlags/PstnUsed", "surface_pressure", PstnUsed_flag);
}
}  // namespace ufo

