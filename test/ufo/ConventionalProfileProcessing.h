/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_CONVENTIONALPROFILEPROCESSING_H_
#define TEST_UFO_CONVENTIONALPROFILEPROCESSING_H_

#include <iomanip>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FloatCompare.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/ConventionalProfileProcessing.h"
#include "ufo/filters/ConventionalProfileProcessingParameters.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variables.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

#include "ufo/profile/EntireSampleDataHandler.h"
#include "ufo/profile/ModelHeightCalculator.h"
#include "ufo/profile/ProfileCheckBackgroundGeopotentialHeight.h"
#include "ufo/profile/ProfileCheckBackgroundRelativeHumidity.h"
#include "ufo/profile/ProfileCheckBackgroundTemperature.h"
#include "ufo/profile/ProfileCheckBackgroundWindSpeed.h"
#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckTime.h"
#include "ufo/profile/ProfileCheckUInterp.h"

#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileDataHolder.h"
#include "ufo/profile/ProfileVerticalAveraging.h"
#include "ufo/profile/ProfileVerticalInterpolation.h"
#include "ufo/profile/VariableNames.h"

#include "ufo/utils/metoffice/MetOfficeQCFlags.h"

namespace ufo {
namespace test {

void testConventionalProfileProcessing(const eckit::LocalConfiguration &conf) {
  util::DateTime bgn(conf.getString("window begin"));
  util::DateTime end(conf.getString("window end"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "obs space");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::world(), bgn, end, oops::mpi::myself());

  const Variables filtervars = Variables(obsspace.obsvariables());

  ObsFilterData filterdata(obsspace);

  ioda::ObsVector hofx(obsspace);

  const eckit::LocalConfiguration obsdiagconf(conf, "obs diagnostics");
  std::vector<eckit::LocalConfiguration> varconfs;
  obsdiagconf.get("variables", varconfs);
  const Variables diagvars(varconfs);
  const ObsDiagnostics obsdiags(obsdiagconf, obsspace, diagvars.toOopsVariables());

  std::shared_ptr<ioda::ObsDataVector<float>> obserr(new ioda::ObsDataVector<float>(
      obsspace, obsspace.obsvariables(), "ObsError"));

  std::shared_ptr<ioda::ObsDataVector<int>> qcflags(new ioda::ObsDataVector<int>(
      obsspace, obsspace.obsvariables()));

  const eckit::LocalConfiguration filterConf(conf, "Conventional Profile Processing");
  ufo::ConventionalProfileProcessingParameters filterParameters;
  filterParameters.validateAndDeserialize(filterConf);

  // Determine whether an exception is expected to be thrown.
  // Exceptions can be thrown in the following places:
  // - on instantiation of the filter,
  // - during the operation of the filter,
  bool expectThrowOnInstantiation = conf.getBool("ExpectThrowOnInstantiation", false);
  bool expectThrowDuringOperation = conf.getBool("ExpectThrowDuringOperation", false);

  if (expectThrowOnInstantiation) {
    EXPECT_THROWS(ufo::ConventionalProfileProcessing filterThrow(obsspace, filterParameters,
                                                         qcflags, obserr));
    // Do not proceed further in this case.
    return;
  }

  // Instantiate filter.
  ufo::ConventionalProfileProcessing filter(obsspace, filterParameters, qcflags, obserr);

  // Obtain GeoVaLs.
  const bool ignoreGeoVaLs = conf.getBool("IgnoreGeoVaLs", false);
  const auto geovars = filter.requiredVars();
  std::unique_ptr <GeoVaLs> geovals;
  if (!ignoreGeoVaLs && geovars.size() > 0) {
    const eckit::LocalConfiguration geovalsConf(conf, "geovals");
    geovals.reset(new GeoVaLs(geovalsConf, obsspace, geovars));
  } else {
    geovals.reset(new GeoVaLs(obsspace.distribution(), oops::Variables()));
  }

  filter.preProcess();
  filter.priorFilter(*geovals);
  if (expectThrowDuringOperation)
    EXPECT_THROWS(filter.postFilter(hofx, obsdiags));
  else
    filter.postFilter(hofx, obsdiags);

  // Determine whether the mismatch check should be bypassed or not.
  // It might be necessary to disable the mismatch check in tests which are
  // designed to reach code paths that would normally result in failure.
  const bool bypassMismatchComparison = conf.getBool("BypassMismatchComparison", false);

  // Check there are no mismatches between the values produced by this code and the OPS equivalents
  if (!bypassMismatchComparison) {
    for (auto nMM : filter.getMismatches())
      EXPECT_EQUAL(nMM, 0);
  }

  // === Additional tests of exceptions === //

  // Test whether adding the same check twice throws an exception.
  const bool addDuplicateCheck = conf.getBool("AddDuplicateCheck", false);
  if (addDuplicateCheck) {
    static ProfileCheckMaker<ProfileCheckUInterp> makerDuplicate1_("duplicate");
    EXPECT_THROWS(static ProfileCheckMaker<ProfileCheckUInterp> makerDuplicate2_("duplicate"));
  }

  // Test whether using the get function with the wrong type throws an exception.
  const bool getWrongType = conf.getBool("GetWrongType", false);
  if (getWrongType) {
    ConventionalProfileProcessingParameters options;
    options.deserialize(conf);
    EntireSampleDataHandler entireSampleDataHandler(obsspace,
                                                    options.DHParameters);
    // Load data from obsspace
    entireSampleDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);
    // Attempt to access data with incorrect type
    EXPECT_THROWS(entireSampleDataHandler.get<int>(ufo::VariableNames::obs_air_pressure));
    std::vector<bool> apply(obsspace.nlocs(), true);
    std::vector<std::vector<bool>> flagged;
    ProfileDataHandler profileDataHandler(filterdata,
                                          options.DHParameters,
                                          apply,
                                          filtervars,
                                          flagged);
    // Obtain profile data
    profileDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);
    // Attempt to access data with incorrect type
    EXPECT_THROWS(profileDataHandler.get<int>(ufo::VariableNames::obs_air_pressure));
  }

  // Manually modify Processing flags in order to cover rare code paths.
  const bool ManualFlagModification = conf.getBool("ManualFlagModification", false);
  if (ManualFlagModification) {
    ConventionalProfileProcessingParameters options;
    options.deserialize(conf);
    EntireSampleDataHandler entireSampleDataHandler(obsspace,
                                                    options.DHParameters);
    // Load data from obsspace
    entireSampleDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);
    entireSampleDataHandler.get<int>(ufo::VariableNames::qcflags_air_temperature);
    std::vector<bool> apply(obsspace.nlocs(), true);
    std::vector<std::vector<bool>> flagged;
    ProfileDataHandler profileDataHandler(filterdata,
                                          options.DHParameters,
                                          apply,
                                          filtervars,
                                          flagged);
    ProfileCheckValidator profileCheckValidator(options);

    profileDataHandler.initialiseNextProfile();

    // Obtain profile data
    profileDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);

    // Modify flags
    std::vector <int> &ReportFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_observation_report);
    std::vector <int> &tFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_air_temperature);
    std::vector <int> &rhFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_relative_humidity);
    std::vector <int> &uFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_eastward_wind);
    std::vector <int> &zFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_geopotential_height);
    std::vector <int> &timeFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_time);

    ReportFlags[0] |= ufo::MetOfficeQCFlags::WholeObReport::PermRejectReport;
    tFlags[0] |= ufo::MetOfficeQCFlags::Profile::SuperadiabatFlag;
    tFlags[0] |= ufo::MetOfficeQCFlags::Profile::InterpolationFlag;
    tFlags[0] |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;
    rhFlags[0] |= ufo::MetOfficeQCFlags::Elem::PermRejectFlag;
    uFlags[0] |= ufo::MetOfficeQCFlags::Profile::InterpolationFlag;
    zFlags[0] |= ufo::MetOfficeQCFlags::Profile::InterpolationFlag;
    zFlags[0] |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;

    // Create checks
    ProfileCheckTime profileCheckTime(options);
    ProfileCheckBackgroundTemperature profileCheckBackgroundT(options);
    ProfileCheckBackgroundRelativeHumidity profileCheckBackgroundRH(options);
    ProfileCheckBackgroundWindSpeed profileCheckBackgroundUV(options);
    ProfileCheckBackgroundGeopotentialHeight profileCheckBackgroundZ(options);

    // Run time check
    profileCheckTime.runCheck(profileDataHandler);

    // Modify time flag
    timeFlags[0] = true;

    // Run remaining checks
    profileCheckBackgroundT.runCheck(profileDataHandler);
    profileCheckBackgroundRH.runCheck(profileDataHandler);
    profileCheckBackgroundUV.runCheck(profileDataHandler);
    profileCheckBackgroundZ.runCheck(profileDataHandler);
  }

  // Test the profile vertical interpolation
  const bool testProfileVerticalInterpolation =
    conf.getBool("testProfileVerticalInterpolation", false);
  if (testProfileVerticalInterpolation) {
    ConventionalProfileProcessingParameters options;
    options.deserialize(conf);

    std::vector<bool> apply(obsspace.nlocs(), true);
    std::vector<std::vector<bool>> flagged;
    ProfileDataHandler profileDataHandler(filterdata,
                                          options.DHParameters,
                                          apply,
                                          filtervars,
                                          flagged);

    // Get interpolation options for each profile.
    const auto interpMethodNames = conf.getStringVector("interpMethodNames");
    const auto coordOrderNames = conf.getStringVector("coordOrderNames");
    const auto outOfBoundsNames = conf.getStringVector("outOfBoundsNames");

    for (size_t jprof = 0; jprof < obsspace.nrecs(); ++jprof) {
      profileDataHandler.initialiseNextProfile();

      const std::string interpMethodName = interpMethodNames[jprof];
      const std::string coordOrderName = coordOrderNames[jprof];
      const std::string outOfBoundsName = outOfBoundsNames[jprof];

      // Calculate level heights for GeoVaLs.
      std::vector <float> orogGeoVaLs(obsspace.nlocs(), 0.0);
      geovals->getAtLevel(orogGeoVaLs, ufo::VariableNames::geovals_orog, 0);
      std::vector <float> zRhoGeoVaLs;
      std::vector <float> zThetaGeoVaLs;
      ufo::CalculateModelHeight(options.DHParameters.ModParameters,
                                orogGeoVaLs[0],
                                zRhoGeoVaLs,
                                zThetaGeoVaLs);

      // Reverse coordinate order if required.
      if (coordOrderName == "Descending")
        std::reverse(zRhoGeoVaLs.begin(), zRhoGeoVaLs.end());

      // Create column of pressure GeoVaLs.
      std::vector <float> pressureGeoVaLs(obsspace.nlocs(), 0.0);
      const size_t gvnlevs = geovals->nlevs(ufo::VariableNames::geovals_pressure);
      std::vector <float> pressureGeoVaLs_column;
      for (int jlev = 0; jlev < gvnlevs; ++jlev) {
        geovals->getAtLevel(pressureGeoVaLs, ufo::VariableNames::geovals_pressure, jlev);
        pressureGeoVaLs_column.push_back(pressureGeoVaLs[0]);
      }

      // Get observed geopotential height and (empty) pressure vector.
      const auto &zObs = profileDataHandler.get<float>(ufo::VariableNames::obs_geopotential_height);
      auto &pressures = profileDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);

      auto interpMethod = ProfileInterpolation::InterpolationMethod::Linear;
      if (interpMethodName == "LogLinear")
        interpMethod = ProfileInterpolation::InterpolationMethod::LogLinear;
      auto coordOrder = ProfileInterpolation::CoordinateOrder::Ascending;
      if (coordOrderName == "Descending")
        coordOrder = ProfileInterpolation::CoordinateOrder::Descending;
      auto outOfBounds = ProfileInterpolation::OutOfBoundsTreatment::SetToBound;
      if (outOfBoundsName == "SetMissing")
        outOfBounds = ProfileInterpolation::OutOfBoundsTreatment::SetMissing;
      if (outOfBoundsName == "Extrapolate")
        outOfBounds = ProfileInterpolation::OutOfBoundsTreatment::Extrapolate;

      // Interpolate to determine pressure.
      profileVerticalInterpolation(zRhoGeoVaLs,
                                   pressureGeoVaLs_column,
                                   zObs,
                                   pressures,
                                   interpMethod,
                                   coordOrder,
                                   outOfBounds);

      // Compare each value of pressure.
      const std::string expectedPressureName = "OPS_" +
        static_cast<std::string>(ufo::VariableNames::obs_air_pressure);
      const auto &expected_pressures = profileDataHandler.get<float>(expectedPressureName);
      for (size_t jlev = 0; jlev < pressures.size(); ++jlev)
        EXPECT(oops::is_close_relative(pressures[jlev], expected_pressures[jlev], 1e-5f));
    }
  }

  // Test the profile vertical averaging.
  const bool testProfileVerticalAveraging =
    conf.getBool("testProfileVerticalAveraging", false);
  if (testProfileVerticalAveraging) {
    ConventionalProfileProcessingParameters options;
    options.deserialize(conf);

    std::vector<bool> apply(obsspace.nlocs(), true);
    std::vector<std::vector<bool>> flagged;
    ProfileDataHandler profileDataHandler(filterdata,
                                          options.DHParameters,
                                          apply,
                                          filtervars,
                                          flagged);

    for (size_t jprof = 0; jprof < obsspace.nrecs(); ++jprof) {
      profileDataHandler.initialiseNextProfile();

      const auto &flagsIn = profileDataHandler.get<int>(ufo::VariableNames::qcflags_eastward_wind);
      const auto &valuesIn = profileDataHandler.get<float>(ufo::VariableNames::obs_eastward_wind);
      const auto &coordIn = profileDataHandler.get<float>(ufo::VariableNames::LogP_derived);
      const auto &bigGap = profileDataHandler.get<float>(ufo::VariableNames::bigPgaps_derived);
      const auto &coordOut = profileDataHandler.get<float>
        (ufo::VariableNames::modellevels_logPWB_rho_derived);
      const float DZFrac = 0.5;
      const ProfileAveraging::Method method =
        ProfileAveraging::Method::Averaging;

      std::vector <int> flagsOut;
      std::vector <float> valuesOut;
      int numGaps = 0;
      std::vector<float> ZMax;
      std::vector<float> ZMin;
      ufo::calculateVerticalAverage(flagsIn,
                                    valuesIn,
                                    coordIn,
                                    bigGap,
                                    coordOut,
                                    DZFrac,
                                    method,
                                    flagsOut,
                                    valuesOut,
                                    numGaps,
                                    &ZMax,
                                    &ZMin);

      // Compare output values with OPS equivalents.
      // The name of the OPS variables are hardcoded because they are purely used for testing.
      // todo(ctgh): check whether any hardcoded names can be substituted.
      const auto &expected_flagsOut =
        profileDataHandler.get<int>("OPS_eastward_wind@ModelLevelsQCFlags");
      for (size_t jlev = 0; jlev < flagsOut.size(); ++jlev)
        EXPECT(flagsOut[jlev] == expected_flagsOut[jlev]);
      const auto &expected_valuesOut =
        profileDataHandler.get<float>("OPS_eastward_wind@ModelLevelsDerivedValue");
      for (size_t jlev = 0; jlev < flagsOut.size(); ++jlev)
        EXPECT(oops::is_close_relative(valuesOut[jlev], expected_valuesOut[jlev], 1e-4f));
      const auto &expected_ZMin =
        profileDataHandler.get<float>("OPS_LogP_u_Min@ModelLevelsDerivedValue");
      for (size_t jlev = 0; jlev < ZMin.size(); ++jlev)
        EXPECT(oops::is_close_relative(ZMin[jlev], expected_ZMin[jlev], 1e-14f));
      const auto &expected_ZMax =
        profileDataHandler.get<float>("OPS_LogP_u_Max@ModelLevelsDerivedValue");
      for (size_t jlev = 0; jlev < ZMax.size(); ++jlev)
        EXPECT(oops::is_close_relative(ZMax[jlev], expected_ZMax[jlev], 1e-14f));
    }
  }

  // Test the profile data holder.
  const bool testProfileDataHolder =
    conf.getBool("testProfileDataHolder", false);
  if (testProfileDataHolder) {
    ConventionalProfileProcessingParameters options;
    options.deserialize(conf);

    std::vector<bool> apply(obsspace.nlocs(), true);
    std::vector<std::vector<bool>> flagged;
    ProfileDataHandler profileDataHandler(filterdata,
                                          options.DHParameters,
                                          apply,
                                          filtervars,
                                          flagged);

    profileDataHandler.initialiseNextProfile();

    // Create a ProfileDataHolder and request some variables.
    ProfileDataHolder profileDataHolder(profileDataHandler);
    profileDataHolder.fill({ufo::VariableNames::extended_obs_space},
                           {ufo::VariableNames::obs_air_pressure},
                           {ufo::VariableNames::station_ID},
                           {ufo::VariableNames::geovals_pressure},
                           {ufo::VariableNames::bkgerr_air_temperature});

    // Get GeoVaLs
    profileDataHolder.getGeoVaLVector(ufo::VariableNames::geovals_pressure);

    // Attempt to access data with incorrect type.
    EXPECT_THROWS(profileDataHolder.get<int>(ufo::VariableNames::obs_air_pressure));

    // Attempt to access nonexistent data.
    EXPECT_THROWS(profileDataHolder.get<int>("wrong@MetaData"));
    EXPECT_THROWS(profileDataHolder.getGeoVaLVector("wrong@MetaData"));
    EXPECT_THROWS(profileDataHolder.getObsDiagVector("wrong@MetaData"));

    // Check this profile has been marked as being in the correct section of the ObsSpace.
    profileDataHolder.checkObsSpaceSection(ufo::ObsSpaceSection::Original);
    EXPECT_THROWS(profileDataHolder.checkObsSpaceSection(ufo::ObsSpaceSection::Extended));

    // Advance to the profile in the extended section and perform the same check.
    profileDataHandler.initialiseNextProfile();
    ProfileDataHolder profileDataHolderExt(profileDataHandler);
    profileDataHolderExt.fill({ufo::VariableNames::extended_obs_space}, {}, {}, {}, {});
    EXPECT_THROWS(profileDataHolderExt.checkObsSpaceSection(ufo::ObsSpaceSection::Original));
    profileDataHolderExt.checkObsSpaceSection(ufo::ObsSpaceSection::Extended);

    // Exercise the set routine
    std::vector <int> int1 {1};
    profileDataHolder.set<int>("test", std::move(int1));
    std::vector <int> int2 {1};
    profileDataHolder.set<int>("test", std::move(int2));
  }
}

class ConventionalProfileProcessing : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::ConventionalProfileProcessing";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
    for (const std::string & testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(), testCaseName);
      ts.emplace_back(CASE("ufo/ConventionalProfileProcessing/" + testCaseName, testCaseConf)
                      {
                        testConventionalProfileProcessing(testCaseConf);
                      });
    }
  }

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_CONVENTIONALPROFILEPROCESSING_H_
