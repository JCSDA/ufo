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
  getObservation("ObsValue", "station_pressure",
                 PStn, true);
  getObservation("ObsValue", "air_pressure",
                 PStd, true);
  getObservation("ObsValue", "mean_sea_level_pressure",
                 Pmsl, true);

  // MetaData
  getObservation("MetaData", "station_elevation",
                 ZStn, true);
  getObservation("MetaData", "standard_elevation",
                 ZStd, true);

  std::vector<int> PStn_flag, PStd_flag, Pmsl_flag;
  getObservation("QCFlags", "station_pressure",
                 PStn_flag, true);
  getObservation("QCFlags", "air_pressure",
                 PStd_flag, true);
  getObservation("QCFlags", "mean_sea_level_pressure",
                 Pmsl_flag, true);

  std::vector<float> PStn_error, PStd_error, Pmsl_error;
  getObservation("ObsError", "station_pressure",
                 PStn_error, true);
  getObservation("ObsError", "air_pressure",
                 PStd_error, true);
  getObservation("ObsError", "mean_sea_level_pressure",
                 Pmsl_error, true);

  std::vector<float> PStn_PGE, PStd_PGE, Pmsl_PGE;
  getObservation("GrossErrorProbability", "station_pressure",
                 PStn_PGE, true);
  getObservation("GrossErrorProbability", "air_pressure",
                 PStd_PGE, true);
  getObservation("GrossErrorProbability", "mean_sea_level_pressure",
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

  // Loop over all obs to calculate PStar
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (!apply[jj]) continue;

    // Set indicators 0 = missing, 1 = rejected, 2 = present
    int IvPStn = 0, IvPStd = 0, IvPmsl = 0;
    if (ZStn[jj] != missing && PStn[jj] != missing) {
      IvPStn = 2;
      if (PStn_flag[jj] & MetOfficeQCFlags::Elem::PermRejectFlag) {IvPStn = 1;}
      BkPStn[jj] = formulas::BackgroundPressure(PSurfParamA[jj], PSurfParamB[jj], ZStn[jj]);
    }

    if (ZStd[jj] != missing && PStd[jj] != missing) {
      IvPStd = 2;
      if (PStd_flag[jj] & MetOfficeQCFlags::Elem::PermRejectFlag) {IvPStd = 1;}
      BkPStd[jj] = formulas::BackgroundPressure(PSurfParamA[jj], PSurfParamB[jj], ZStd[jj]);
    }

    if (Pmsl[jj] != missing) {
      IvPmsl = 2;
      if (Pmsl_flag[jj] & MetOfficeQCFlags::Elem::PermRejectFlag) {IvPmsl = 1;}
      BkPmsl[jj] = formulas::BackgroundPressure(PSurfParamA[jj], PSurfParamB[jj], 0.0);
    }

    int IvPmax = std::max(IvPStn, std::max(IvPStd, IvPmsl));

    if (IvPmax == 0) continue;
    if (IvPStn == IvPmax && BkPStn[jj] != missing && PStn_flag[jj]
        & MetOfficeQCFlags::Surface::PstnPrefFlag) {
      PStar[jj] = (BkPStar[jj]*PStn[jj])/BkPStn[jj];
      PStar_flag[jj] = PStn_flag[jj];
      PStar_flag[jj] |= ufo::MetOfficeQCFlags::Surface::PstnUsedFlag;
      PStar_error[jj] = PStn_error[jj];
      PStar_PGE[jj] = PStn_PGE[jj];
    } else if (BkPStd[jj] != missing && BkPmsl[jj] == missing) {
      PStar[jj] = (BkPStar[jj]*PStd[jj])/BkPStd[jj];
      PStar_flag[jj] = PStd_flag[jj];
      PStar_flag[jj] |= ufo::MetOfficeQCFlags::Surface::PstdUsedFlag;
      PStar_error[jj] = PStd_error[jj];
      PStar_PGE[jj] = PStd_PGE[jj];
    } else if (BkPmsl[jj] != missing) {
      PStar[jj] = BkPStar[jj]*Pmsl[jj]/BkPmsl[jj];
      PStar_flag[jj] = Pmsl_flag[jj];
      PStar_flag[jj] |= ufo::MetOfficeQCFlags::Surface::PmslUsedFlag;
      PStar_error[jj] = Pmsl_error[jj];
      PStar_PGE[jj] = Pmsl_PGE[jj];
    } else {
        continue;
    }
  }

  putObservation("PStar", PStar);
  obsdb_.put_db("DerivedQCFlags", "PStar", PStar_flag);
  obsdb_.put_db("DerivedObsError", "PStar", PStar_error);
  obsdb_.put_db("DerivedGrossErrorProbability", "PStar", PStar_PGE);
}
}  // namespace ufo

