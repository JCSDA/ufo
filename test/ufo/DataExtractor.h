/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_DATAEXTRACTOR_H_
#define TEST_UFO_DATAEXTRACTOR_H_

#include "ufo/utils/dataextractor/DataExtractor.h"

#include <iomanip>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "eckit/testing/Test.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FloatCompare.h"

namespace ufo {
namespace test {

float missing = util::missingValue(missing);


template <typename T, typename R>
float run_basic(const T obVal0, const R obVal1, const std::vector<T> &varValues0,
                const std::vector<R> &varValues1) {
    const int dimIndex0 = 0;
    const int dimIndex1 = 1;
    const std::string &varName0 = "var0";
    const std::string &varName1 = "var1";
    const std::array<ConstrainedRange, 2> ranges {
        ConstrainedRange(varValues0.size()),
        ConstrainedRange(varValues1.size())};
    assert(varValues0.size() == 5);
    assert(varValues1.size() == 3);

    boost::multi_array<float, 2> interpolatedArray(boost::extents[5][3]);
    std::vector<float> tmp = {1, 2, 3, 4, 5,
                              6, 7, 8, 9, 10,
                              11, 12, 13, 14, 15};
    size_t ind = 0;
    for (int j=0; j < interpolatedArray.shape()[1]; j++) {
      for (int i=0; i < interpolatedArray.shape()[0]; i++) {
        interpolatedArray[i][j] = tmp[ind];
        ind++;
      }
    }
    return bilinearInterpolation(dimIndex0, varName0, varValues0, obVal0,
                                 dimIndex1, varName1, varValues1, obVal1,
                                 ranges, interpolatedArray);
}


// For int/float calls
template <typename T, typename R>
float run_basic(const T obVal0, const R obVal1) {
  const std::vector<T> varValues0 {2, 4, 6, 8, 10};
  const std::vector<R> varValues1 {2, 4, 6};
  return run_basic(obVal0, obVal1, varValues0, varValues1);
}


CASE("ufo/DataExtractor/bilinearinterp/float_linear") {
  // Effectively becomes linear interpolation along dim1.
  const float res = run_basic(4.0, 4.2);
  EXPECT(oops::is_close_absolute(res, 7.5f, 1e-5f, 0,
                                 oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE));
}


CASE("ufo/DataExtractor/bilinearinterp/float_linear_at_lower_boundary_dim0") {
  // Check handling where our point is located on the lower boundary.
  const float res = run_basic(2, 4);
  EXPECT(oops::is_close_absolute(res, 6.0f, 1e-5f, 0,
                                 oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE));
}


CASE("ufo/DataExtractor/bilinearinterp/float_linear_at_lower_boundary_dim1") {
  // Check handling where our point is located on the lower boundary.
  const float res = run_basic(4, 2);
  EXPECT(oops::is_close_absolute(res, 2.0f, 1e-5f, 0,
                                 oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE));
}


CASE("ufo/DataExtractor/bilinearinterp/float_blinear") {
  // Simple bilinear interpolation.
  const float res = run_basic(4.2, 4.2);
  EXPECT(oops::is_close_absolute(res, 7.6f, 1e-5f, 0,
                                 oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE));
}


CASE("ufo/DataExtractor/bilinearinterp/extrapolation_lower_bound_dim0") {
  // Lower bound extrapolation dim0
  const float res = run_basic(0, 4.2);
  EXPECT_EQUAL(res, missing);
}


CASE("ufo/DataExtractor/bilinearinterp/extrapolation_lower_bound_dim1") {
  // Lower bound extrapolation dim1
  const float res = run_basic(4.2, 0.0);
  EXPECT_EQUAL(res, missing);
}


CASE("ufo/DataExtractor/bilinearinterp/extrapolation_upper_bound_dim0") {
  // Upper bound extrapolation dim0
  const float res = run_basic(20.0, 4.2);
  EXPECT_EQUAL(res, missing);
}


CASE("ufo/DataExtractor/bilinearinterp/extrapolation_upper_bound_dim1") {
  // Upper bound extrapolation dim1
  const float res = run_basic(4.2, 20);
  EXPECT_EQUAL(res, missing);
}


CASE("ufo/DataExtractor/bilinearinterp/int_int_dtype") {
  const float res = run_basic(3, 3);
  EXPECT(oops::is_close_absolute(res, 4.0f, 1e-5f, 0,
                                 oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE));
}


CASE("ufo/DataExtractor/bilinearinterp/int_float_dtype") {
  const float res = run_basic(3, 3.0);
  EXPECT(oops::is_close_absolute(res, 4.0f, 1e-5f, 0,
                                 oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE));
}

CASE("ufo/DataExtractor/bilinearinterp/float_int_dtype") {
  const float res = run_basic(3.0, 3);
  EXPECT(oops::is_close_absolute(res, 4.0f, 1e-5f, 0,
                                 oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE));
}


CASE("ufo/DataExtractor/bilinearinterp/string_dtype") {
  const std::vector<std::string> strVarValues0 {"2", "4", "6", "8", "10"};
  const std::vector<std::string> strVarValues1 {"2", "4", "6"};
  const std::vector<float> floatVarValues0 {2, 4, 6, 8, 10};
  const std::vector<float> floatVarValues1 {2, 4, 6};
  const std::string msg = "Bilinear interpolation cannot be performed along coordinate axes "
                          "indexed by string variables such as ";
  EXPECT_THROWS_MSG(run_basic(4.2f, std::string("4.2"), floatVarValues0, strVarValues1),
                    (msg + "var1.").c_str());
  EXPECT_THROWS_MSG(run_basic(std::string("4.2"), 4.2f, strVarValues0, floatVarValues1),
                    (msg + "var0.").c_str());
  EXPECT_THROWS_MSG(run_basic(std::string("4.2"), std::string("4.2"),
                              strVarValues0, strVarValues1),
                    (msg + "var0 or var1.").c_str());
}


float run_missing(const float obVal0, const float obVal1, const std::vector<float> data) {
  const std::vector<float> varValues0 {2, 4, 6, 8, 10};
  const std::vector<float> varValues1 {2, 4, 6};
  const std::array<ConstrainedRange, 2> ranges {
      ConstrainedRange(varValues0.size()),
      ConstrainedRange(varValues1.size())};

  boost::multi_array<float, 2> interpolatedArray(boost::extents[5][3]);
  size_t ind = 0;
  for (int j=0; j < interpolatedArray.shape()[1]; j++) {
    for (int i=0; i < interpolatedArray.shape()[0]; i++) {
      interpolatedArray[i][j] = data[ind];
      ind++;
    }
  }
  return bilinearInterpolation(0, "var0", varValues0, obVal0,
                               1, "var1", varValues1, obVal1,
                               ranges, interpolatedArray);
}


CASE("ufo/DataExtractor/bilinearinterp/one_missing") {
  // If one missing, pick closes of non missing neighbours.
  const std::vector<float> data = {missing, 2, 3, 4, 5,
                                   6, 7, 8, 9, 10,
                                   11, 12, 13, 14, 15};
  // Pick closes non-missing neighbour.
  float res = run_missing(3.1, 2.5, data);
  EXPECT(oops::is_close_absolute(res, 2.0f, 1e-5f, 0,
                                 oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE));

  // Different neighbour is closest.
  res = run_missing(2.5, 3.1, data);
  EXPECT(oops::is_close_absolute(res, 6.0f, 1e-5f, 0,
                                 oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE));

  // Avoid missing closest value.
  res = run_missing(2.5, 2.1, data);
  EXPECT(oops::is_close_absolute(res, 2.0f, 1e-5f, 0,
                                 oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE));
}


CASE("ufo/DataExtractor/bilinearinterp/all_missing") {
  // If one missing, pick closest of non missing neighbours.
  const std::vector<float> data = {missing, missing, 3, 4, 5,
                                   missing, missing, 8, 9, 10,
                                   11, 12, 13, 14, 15};
  const float res = run_missing(3.1, 2.5, data);
  EXPECT_EQUAL(res, missing);
}


CASE("ufo/DataExtractor/bilinearinterp/range_constrain") {
  // Let's constrain the range and change the dimension mapping.
  const float obVal0 = 3.0, obVal1 = 3.0;
  // Coordinates containing values that would change the answer if not excluded
  // via the "range" specified.
  const std::vector<float> varValues0 {2, 4, 2, 4, 6};
  const std::vector<float> varValues1 {3, 2, 4};

  ConstrainedRange con0 = ConstrainedRange(varValues0.size());
  con0.constrain(2, varValues0.size());  // range will ignore 1st two.
  ConstrainedRange con1 = ConstrainedRange(varValues1.size());
  con1.constrain(1, varValues1.size());  // range will ignore 1st.
  const std::array<ConstrainedRange, 2> ranges {con1, con0};

  const std::vector<float> data = {1, 6, 11,
                                   2, 7, 12,
                                   3, 8, 13,
                                   4, 9, 14,
                                   5, 10, 15};
  boost::multi_array<float, 2> interpolatedArray(boost::extents[3][5]);
  size_t ind = 0;
  for (int j=0; j < interpolatedArray.shape()[1]; j++) {
    for (int i=0; i < interpolatedArray.shape()[0]; i++) {
      interpolatedArray[i][j] = data[ind];
      ind++;
    }
  }
  const float res = bilinearInterpolation(1, "var0", varValues0, obVal0,
                                          0, "var1", varValues1, obVal1,
                                          ranges, interpolatedArray);
  EXPECT(oops::is_close_absolute(res, 11.0f, 1e-5f, 0,
                                 oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE));
}


CASE("ufo/DataExtractor/bilinearinterp/range_constrain_extrapolation") {
  // Constraining the range, such that, what would otherwise be within bounds is now out of bounds
  // so returns missing.
  const float obVal0 = 5.0, obVal1 = 3.0;
  // Coordinates containing values that would change the answer if not excluded
  // via the "range" specified.
  const std::vector<float> varValues0 {2, 4, 2, 4, 6};
  const std::vector<float> varValues1 {3, 2, 4};

  ConstrainedRange con0 = ConstrainedRange(varValues0.size());
  con0.constrain(2, varValues0.size()-1);  // range will ignore 1st two and the final element.
  ConstrainedRange con1 = ConstrainedRange(varValues1.size());
  con1.constrain(1, varValues1.size());  // range will ignore 1st.
  const std::array<ConstrainedRange, 2> ranges {con1, con0};

  const std::vector<float> data = {1, 6, 11,
                                   2, 7, 12,
                                   3, 8, 13,
                                   4, 9, 14,
                                   5, 10, 15};
  boost::multi_array<float, 2> interpolatedArray(boost::extents[3][5]);
  size_t ind = 0;
  for (int j=0; j < interpolatedArray.shape()[1]; j++) {
    for (int i=0; i < interpolatedArray.shape()[0]; i++) {
      interpolatedArray[i][j] = data[ind];
      ind++;
    }
  }
  const float res = bilinearInterpolation(1, "var0", varValues0, obVal0,
                                          0, "var1", varValues1, obVal1,
                                          ranges, interpolatedArray);
  EXPECT_EQUAL(res, missing);
}


class DataExtractor : public oops::Test {
 public:
  DataExtractor() {}

 private:
  std::string testid() const override {return "ufo::test::DataExtractor";}

  void register_tests() const override {}

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_DATAEXTRACTOR_H_
