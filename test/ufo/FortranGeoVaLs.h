/*
 * (C) Crown copyright 2023, UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_FORTRANGEOVALS_H_
#define TEST_UFO_FORTRANGEOVALS_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------
extern "C" {
  void test_ufo_geovals_setup_with_mismatched_nvars_f90();
  void test_ufo_geovals_setup_with_mismatched_nreduced_vars_f90();
  void test_ufo_geovals_setup_with_sampling_method_set_to_0_f90();
  void test_ufo_geovals_setup_with_sampling_method_set_to_3_f90();
  void test_ufo_geovals_setup_with_method_mislabelled_as_trivial_f90();
  void test_ufo_geovals_partial_setup_with_mismatched_nvars_f90();
  void test_ufo_geovals_partial_setup_with_sampling_method_set_to_0_f90();
  void test_ufo_geovals_partial_setup_with_sampling_method_set_to_3_f90();
  void test_ufo_geovals_setup_sampling_method_0_f90();
  void test_ufo_geovals_setup_sampling_method_3_f90();
  void test_ufo_geovals_setup_sampling_method_with_mismatched_nlocs_f90();
  void test_ufo_geovals_setup_sampling_method_mislabelled_as_trivial_f90();
  void test_ufo_geovals_setup_trivial_sampling_method_0_f90();
  void test_ufo_geovals_setup_trivial_sampling_method_3_f90();
  void test_ufo_geovals_set_invalid_default_format_f90();
  void test_ufo_geovals_get_nonexistent_var_f90();
}

// -----------------------------------------------------------------------------

/// Test error checking in Fortran. Each of these test cases is expected to cause the program
/// to print an error message and abort. The test will pass if the error message matches a regular
/// expression defined in CMakeLists.txt.
void abortingTests() {
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  if (conf.has("aborting test")) {
    const std::string name = conf.getString("aborting test");
    if (name == "ufo_geovals_setup_with_mismatched_nvars")
      test_ufo_geovals_setup_with_mismatched_nvars_f90();
    else if (name == "ufo_geovals_setup_with_mismatched_nreduced_vars")
      test_ufo_geovals_setup_with_mismatched_nreduced_vars_f90();
    else if (name == "ufo_geovals_setup_with_sampling_method_set_to_0")
      test_ufo_geovals_setup_with_sampling_method_set_to_0_f90();
    else if (name == "ufo_geovals_setup_with_sampling_method_set_to_3")
      test_ufo_geovals_setup_with_sampling_method_set_to_3_f90();
    else if (name == "ufo_geovals_setup_with_method_mislabelled_as_trivial")
      test_ufo_geovals_setup_with_method_mislabelled_as_trivial_f90();
    else if (name == "ufo_geovals_partial_setup_with_mismatched_nvars")
      test_ufo_geovals_partial_setup_with_mismatched_nvars_f90();
    else if (name == "ufo_geovals_partial_setup_with_sampling_method_set_to_0")
      test_ufo_geovals_partial_setup_with_sampling_method_set_to_0_f90();
    else if (name == "ufo_geovals_partial_setup_with_sampling_method_set_to_3")
      test_ufo_geovals_partial_setup_with_sampling_method_set_to_3_f90();
    else if (name == "ufo_geovals_setup_sampling_method_0")
      test_ufo_geovals_setup_sampling_method_0_f90();
    else if (name == "ufo_geovals_setup_sampling_method_3")
      test_ufo_geovals_setup_sampling_method_3_f90();
    else if (name == "ufo_geovals_setup_sampling_method_with_mismatched_nlocs")
      test_ufo_geovals_setup_sampling_method_with_mismatched_nlocs_f90();
    else if (name == "ufo_geovals_setup_sampling_method_mislabelled_as_trivial")
      test_ufo_geovals_setup_sampling_method_mislabelled_as_trivial_f90();
    else if (name == "ufo_geovals_setup_trivial_sampling_method_0")
      test_ufo_geovals_setup_trivial_sampling_method_0_f90();
    else if (name == "ufo_geovals_setup_trivial_sampling_method_3")
      test_ufo_geovals_setup_trivial_sampling_method_3_f90();
    else if (name == "ufo_geovals_set_invalid_default_format")
      test_ufo_geovals_set_invalid_default_format_f90();
    else if (name == "ufo_geovals_get_nonexistent_var")
      test_ufo_geovals_get_nonexistent_var_f90();
    else
      throw std::runtime_error("Invalid test name: " + name);
  }
}

// -----------------------------------------------------------------------------

class FortranGeoVaLs : public oops::Test {
 public:
  FortranGeoVaLs() {}
  virtual ~FortranGeoVaLs() = default;

 private:
  std::string testid() const override {return "ufo::test::FortranGeoVaLs";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/FortranGeoVaLs/abortingTests")
      { abortingTests(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_FORTRANGEOVALS_H_
