/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_METOFFICERADIANCEERRORMATRICES_H_
#define TEST_UFO_METOFFICERADIANCEERRORMATRICES_H_

#include <Eigen/Dense>
#include <algorithm>
#include <set>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"
#include "ufo/utils/metoffice/MetOfficeBMatrixStatic.h"
#include "ufo/utils/metoffice/MetOfficeRMatrixRadiance.h"

namespace ufo {
namespace test {

void testMetOfficeRadianceErrorMatrices(const eckit::LocalConfiguration &conf) {
  // Get test values from configuration
  std::string chlist = conf.getString("channels");
  const float latitude = conf.getFloat("latitude");
  const size_t nelements = conf.getUnsigned("nelements");
  const float BHT_value = conf.getFloat("BHT_value");
  const float HBHT_value = conf.getFloat("HBHT_value");
  const float HBHT_R_value = conf.getFloat("HBHT_R_value");
  const float tol = conf.getFloat("tol");

  // Setup local variables
  std::vector<int> channels;
  std::set<int> chans = oops::parseIntSet(chlist);
  std::copy(chans.begin(), chans.end(), std::back_inserter(channels));
  const size_t nchans = channels.size();

  // -----------------
  // Bmatrix testing
  // -----------------

  // Construct object
  MetOfficeBMatrixStatic bmatrix(conf);

  // Get an eigen matrix for a specific latitude
  Eigen::MatrixXf Hmatrix = Eigen::MatrixXf::Constant(nchans, nelements, 2.0);
  Eigen::MatrixXf BHT;
  bmatrix.multiply(latitude, Hmatrix.transpose(), BHT);

  // Test print function
  oops::Log::info() << "bmatrix = " << bmatrix << std::endl;

  // -----------------
  // Rmatrix testing
  // -----------------

  // Construct object
  MetOfficeRMatrixRadiance rmatrix(conf);

  // Create a matrix to add to
  Eigen::MatrixXf HBHT;
  Eigen::MatrixXf HBHT_R;
  HBHT = Hmatrix * BHT;
  // HBH' + R
  rmatrix.add(channels, HBHT, HBHT_R);

  // Test print function
  oops::Log::info() << "rmatrix = " << rmatrix << std::endl;

  // -----------------------------
  // Check values during testing
  // ----------------------------

  // Test to check the multiplication appears to be correct
  ASSERT(std::abs(BHT(0, 0) - BHT_value) < tol);
  ASSERT(std::abs(HBHT(0, 0) - HBHT_value) < tol);
  ASSERT(std::abs(HBHT_R(0, 0) - HBHT_R_value) < tol);
}

class MetOfficeRadianceErrorMatrices : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::MetOfficeRadianceErrorMatrices";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
    for (const std::string & testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(), testCaseName);
      ts.emplace_back(CASE("ufo/MetOfficeRadianceErrorMatrices/" + testCaseName, testCaseConf)
                      {
                        testMetOfficeRadianceErrorMatrices(testCaseConf);
                      });
    }
  }

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_METOFFICERADIANCEERRORMATRICES_H_
