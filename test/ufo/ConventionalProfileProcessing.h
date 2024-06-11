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
#include "ioda/ObsDataVector.h"
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
#include "ufo/profile/VariableNames.h"

#include "ufo/utils/metoffice/MetOfficeQCFlags.h"

namespace ufo {
namespace test {

void testConventionalProfileProcessing(const eckit::LocalConfiguration &conf) {
  util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "obs space");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::world(), timeWindow, oops::mpi::myself());

  const Variables filtervars = Variables(obsspace.obsvariables());

  ObsFilterData filterdata(obsspace);

  ioda::ObsVector hofx(obsspace, "HofX");

  ioda::ObsVector bias(obsspace);
  bias.zero();

  const eckit::LocalConfiguration obsdiagconf = conf.getSubConfiguration("obs diagnostics");
  std::vector<eckit::LocalConfiguration> varconfs;
  conf.get("obs diagnostics variables", varconfs);
  const Variables diagvars(varconfs);
  const ObsDiagnostics obsdiags(obsdiagconf, obsspace, diagvars.toOopsObsVariables());

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
  // - during the operation of the filter at the prior stage
  // - during the operation of the filter at the post stage
  const bool expectThrowOnInstantiation = conf.getBool("ExpectThrowOnInstantiation", false);
  const bool expectThrowDuringPriorFilter = conf.getBool("ExpectThrowDuringPriorFilter", false);
  const bool expectThrowDuringPostFilter = conf.getBool("ExpectThrowDuringPostFilter", false);

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
  geovals->setDefaultFormat(GeoVaLFormat::REDUCED);

  filter.preProcess();
  if (expectThrowDuringPriorFilter) {
    EXPECT_THROWS(filter.priorFilter(*geovals));
    return;
  }
  filter.priorFilter(*geovals);
  if (expectThrowDuringPostFilter) {
    EXPECT_THROWS(filter.postFilter(*geovals, hofx, bias, obsdiags));
    return;
  }
  filter.postFilter(*geovals, hofx, bias, obsdiags);

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
    EntireSampleDataHandler entireSampleDataHandler(filterdata,
                                                    options.DHParameters);
    // Load data from obsspace
    entireSampleDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);
    // Attempt to access data with incorrect type
    EXPECT_THROWS(entireSampleDataHandler.get<int>(ufo::VariableNames::obs_air_pressure));
    std::vector<bool> apply(obsspace.nlocs(), true);
    std::vector<std::vector<bool>> flagged;
    ProfileDataHandler profileDataHandler(filterdata,
                                          *qcflags,
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
    EntireSampleDataHandler entireSampleDataHandler(filterdata,
                                                    options.DHParameters);
    // Load data from obsspace
    entireSampleDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);
    entireSampleDataHandler.get<int>(ufo::VariableNames::qcflags_air_temperature);
    std::vector<bool> apply(obsspace.nlocs(), true);
    std::vector<std::vector<bool>> flagged;
    ProfileDataHandler profileDataHandler(filterdata,
                                          *qcflags,
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

    // Run remaining checks
    profileCheckBackgroundT.runCheck(profileDataHandler);
    profileCheckBackgroundRH.runCheck(profileDataHandler);
    profileCheckBackgroundUV.runCheck(profileDataHandler);
    profileCheckBackgroundZ.runCheck(profileDataHandler);
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
                                          *qcflags,
                                          options.DHParameters,
                                          apply,
                                          filtervars,
                                          flagged);

    for (size_t jprof = 0; jprof < obsspace.nrecs(); ++jprof) {
      profileDataHandler.initialiseNextProfile();

      const auto &flagsIn = profileDataHandler.get<int>(ufo::VariableNames::qcflags_eastward_wind);
      const auto &valuesIn = profileDataHandler.get<float>(ufo::VariableNames::obs_eastward_wind);
      const auto &coordIn = profileDataHandler.get<float>("logP@DerivedValue");
      const auto &bigGap = profileDataHandler.get<float>("bigPgaps@DerivedValue");
      const auto &coordOut = profileDataHandler.get<float>
        ("LogPWB@ModelRhoLevelsDerivedValue");
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
      const auto &expected_flagsOut =
        profileDataHandler.get<int>("ModelLevelsQCFlags/OPS_windEastward");
      for (size_t jlev = 0; jlev < flagsOut.size(); ++jlev)
        EXPECT(flagsOut[jlev] == expected_flagsOut[jlev]);
      const auto &expected_valuesOut =
        profileDataHandler.get<float>("ModelLevelsDerivedValue/OPS_windEastward");
      for (size_t jlev = 0; jlev < flagsOut.size(); ++jlev)
        EXPECT(oops::is_close_relative(valuesOut[jlev], expected_valuesOut[jlev], 1e-4f));
      const auto &expected_ZMin =
        profileDataHandler.get<float>("ModelLevelsDerivedValue/OPS_LogP_u_Min");
      for (size_t jlev = 0; jlev < ZMin.size(); ++jlev)
        EXPECT(oops::is_close_relative(ZMin[jlev], expected_ZMin[jlev], 1e-14f));
      const auto &expected_ZMax =
        profileDataHandler.get<float>("ModelLevelsDerivedValue/OPS_LogP_u_Max");
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
                                          *qcflags,
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
                           oops::Variables{{oops::Variable{ufo::VariableNames::geovals_pressure}}});

    // Get GeoVaLs
    profileDataHolder.getGeoVaLVector(oops::Variable{ufo::VariableNames::geovals_pressure});

    // Attempt to access data with incorrect type.
    EXPECT_THROWS(profileDataHolder.get<int>(ufo::VariableNames::obs_air_pressure));

    // Attempt to access nonexistent data.
    EXPECT_THROWS(profileDataHolder.get<int>("wrong@MetaData"));
    EXPECT_THROWS(profileDataHolder.getGeoVaLVector(oops::Variable{"wrong@MetaData"}));

    // Check this profile has been marked as being in the correct section of the ObsSpace.
    profileDataHolder.checkObsSpaceSection(ufo::ObsSpaceSection::Original);
    EXPECT_THROWS(profileDataHolder.checkObsSpaceSection(ufo::ObsSpaceSection::Extended));

    // Advance to the profile in the extended section and perform the same check.
    profileDataHandler.initialiseNextProfile();
    ProfileDataHolder profileDataHolderExt(profileDataHandler);
    profileDataHolderExt.fill({ufo::VariableNames::extended_obs_space}, {}, {}, {});
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
