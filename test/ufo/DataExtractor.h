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
    const std::array<ConstrainedRange, 3> ranges {
        ConstrainedRange(varValues0.size()),
        ConstrainedRange(varValues1.size()),
        ConstrainedRange(1)};
    assert(varValues0.size() == 5);
    assert(varValues1.size() == 3);

    boost::multi_array<float, 3> interpolatedArray(boost::extents[5][3][1]);
    std::vector<float> tmp = {1, 2, 3, 4, 5,
                              6, 7, 8, 9, 10,
                              11, 12, 13, 14, 15};
    size_t ind = 0;
    for (int j=0; j < interpolatedArray.shape()[1]; j++) {
      for (int i=0; i < interpolatedArray.shape()[0]; i++) {
        interpolatedArray[i][j][0] = tmp[ind];
        ind++;
      }
    }
    auto array = get2DSlice(interpolatedArray, dimIndex0, dimIndex1, ranges);
    return bilinearInterpolation(varName0, varValues0, obVal0, ranges[dimIndex0],
                                 varName1, varValues1, obVal1, ranges[dimIndex1],
                                 array);
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
  const float res = run_basic(4.0f, 4.2f);
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
  const std::string msg = "No match found for 'bilinear' interpolation of value '0' of the "
                          "variable 'var0'.  Value is out of bounds.  Consider using "
                          "extrapolation.";
  EXPECT_THROWS_MSG(run_basic(0, 4.2), msg.c_str());
}


CASE("ufo/DataExtractor/bilinearinterp/extrapolation_lower_bound_dim1") {
  // Lower bound extrapolation dim1
  const std::string msg = "No match found for 'bilinear' interpolation of value '0' of the "
                          "variable 'var1'.  Value is out of bounds.  Consider using "
                          "extrapolation.";
  EXPECT_THROWS_MSG(run_basic(4.2, 0.0), msg.c_str());
}


CASE("ufo/DataExtractor/bilinearinterp/extrapolation_upper_bound_dim0") {
  // Upper bound extrapolation dim0
  const std::string msg = "No match found for 'bilinear' interpolation of value '20' of the "
                          "variable 'var0'.  Value is out of bounds.  Consider using "
                          "extrapolation.";
  EXPECT_THROWS_MSG(run_basic(20.0, 4.2), msg.c_str());
}


CASE("ufo/DataExtractor/bilinearinterp/extrapolation_upper_bound_dim1") {
  // Upper bound extrapolation dim1
  const std::string msg = "No match found for 'bilinear' interpolation of value '20' of the "
                          "variable 'var1'.  Value is out of bounds.  Consider using "
                          "extrapolation.";
  EXPECT_THROWS_MSG(run_basic(4.2, 20), msg.c_str());
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
  const std::array<ConstrainedRange, 3> ranges {
      ConstrainedRange(varValues0.size()),
      ConstrainedRange(varValues1.size()),
      ConstrainedRange(1)};

  boost::multi_array<float, 3> interpolatedArray(boost::extents[5][3][1]);
  size_t ind = 0;
  for (int j=0; j < interpolatedArray.shape()[1]; j++) {
    for (int i=0; i < interpolatedArray.shape()[0]; i++) {
      interpolatedArray[i][j][0] = data[ind];
      ind++;
    }
  }
  auto array = get2DSlice(interpolatedArray, 0, 1, ranges);
  return bilinearInterpolation("var0", varValues0, obVal0, ranges[0],
                               "var1", varValues1, obVal1, ranges[1],
                               array);
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


float run_range_constrained(const float obVal0, const float obVal1,
                            const std::array<ConstrainedRange, 3> &ranges) {
  // Coordinates containing values that would change the answer if not excluded
  // via the "range" specified.
  const std::vector<float> varValues0 {2, 4, 2, 4, 6};
  const std::vector<float> varValues1 {3, 2, 4};

  const std::vector<float> data = {1, 6, 11,
                                   2, 7, 12,
                                   3, 8, 13,
                                   4, 9, 14,
                                   5, 10, 15};
  boost::multi_array<float, 3> interpolatedArray(boost::extents[3][5][1]);
  size_t ind = 0;
  for (int j=0; j < interpolatedArray.shape()[1]; j++) {
    for (int i=0; i < interpolatedArray.shape()[0]; i++) {
      interpolatedArray[i][j][0] = data[ind];
      ind++;
    }
  }
  auto array = get2DSlice(interpolatedArray, 1, 0, ranges);
  return bilinearInterpolation("var1", varValues1, obVal1, ranges[0],
                               "var0", varValues0, obVal0, ranges[1],
                               array);
}


CASE("ufo/DataExtractor/bilinearinterp/range_constrain") {
  // Let's constrain the range and change the dimension mapping.
  const float obVal0 = 3.0, obVal1 = 3.0;

  ConstrainedRange con0 = ConstrainedRange(5);
  con0.constrain(2, 5);  // range will ignore 1st two.
  ConstrainedRange con1 = ConstrainedRange(3);
  con1.constrain(1, 3);  // range will ignore 1st.
  ConstrainedRange con2 = ConstrainedRange(1);
  const std::array<ConstrainedRange, 3> ranges {con1, con0, con2};

  const float res = run_range_constrained(obVal0, obVal1, ranges);
  EXPECT(oops::is_close_absolute(res, 11.0f, 1e-5f, 0,
                                 oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE));
}


CASE("ufo/DataExtractor/bilinearinterp/range_constrain_extrapolation") {
  // Constraining the range, such that our extract is now out of bounds.
  const float obVal0 = 5.0, obVal1 = 3.0;

  ConstrainedRange con0 = ConstrainedRange(5);
  con0.constrain(2, 4);  // range will ignore 1st two and the final element.
  ConstrainedRange con1 = ConstrainedRange(3);
  con1.constrain(1, 3);  // range will ignore 1st.
  ConstrainedRange con2 = ConstrainedRange(1);
  const std::array<ConstrainedRange, 3> ranges {con1, con0, con2};

  const std::string msg = "No match found for 'bilinear' interpolation of value '5' of the "
                          "variable 'var0'.  Value is out of bounds.  Consider using "
                          "extrapolation.";
  EXPECT_THROWS_MSG(run_range_constrained(obVal0, obVal1, ranges),
                    msg.c_str());
}


CASE("ufo/DataExtractor/bilinearinterp/range_constrain_3D_array") {
  // Constraining the second dimension. array shape: (5, 6, 3); 2D slice: (:, 2, :)
  const std::vector<float> varValues0 {2, 4, 6, 8, 10};
  const std::vector<float> varValues1 {2, 4, 6};
  const int dimIndex0 = 0;
  const int dimIndex1 = 2;
  const std::string &varName0 = "var0";
  const std::string &varName1 = "var1";
  std::array<ConstrainedRange, 3> ranges {
    ConstrainedRange(varValues0.size()),
    ConstrainedRange(6),
    ConstrainedRange(varValues1.size())};
  ranges[1].constrain(2, 3);

  boost::multi_array<float, 3> interpolatedArray(boost::extents[5][6][3]);
  const std::vector<float> tmp = {1, 2, 3, 4, 5,
                                  6, 7, 8, 9, 10,
                                  11, 12, 13, 14, 15};
  size_t ind = 0;
  for (int j=0; j < interpolatedArray.shape()[2]; j++) {
    for (int i=0; i < interpolatedArray.shape()[0]; i++) {
      interpolatedArray[i][2][j] = tmp[ind];
      ind++;
    }
  }
  auto array = get2DSlice(interpolatedArray, dimIndex0, dimIndex1, ranges);
  const float res = bilinearInterpolation(varName0, varValues0, 4.2f, ranges[dimIndex0],
                                          varName1, varValues1, 4.2f, ranges[dimIndex1],
                                          array);
  EXPECT(oops::is_close_absolute(res, 7.6f, 1e-5f, 0,
                                 oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE));
}


CASE("ufo/DataExtractor/get2DSlice/not_2d_slice") {
  // The range specified means that we don't have a single 2D slice of the array
  // All dimensions here are unconstrained.
  const std::array<ConstrainedRange, 3> ranges {ConstrainedRange(2), ConstrainedRange(2),
                                                ConstrainedRange(2)};
  boost::multi_array<float, 3> interpolatedArray(boost::extents[2][2][2]);

  const std::string msg = "Unable to fetch a 2D array slice with the provided constraints.";
  EXPECT_THROWS_MSG(get2DSlice(interpolatedArray, 0, 1, ranges), msg.c_str());
}


CASE("ufo/DataExtractor/get1DSlice/not_1d_slice") {
  // The range specified means that we don't have a single 1D slice of the array
  // All dimensions here are unconstrained.
  const std::array<ConstrainedRange, 3> ranges {ConstrainedRange(2), ConstrainedRange(2),
                                                ConstrainedRange(2)};
  boost::multi_array<float, 3> interpolatedArray(boost::extents[2][2][2]);

  const std::string msg = "Unable to fetch a 1D array slice with the provided constraints.";
  EXPECT_THROWS_MSG(get1DSlice(interpolatedArray, 0, ranges), msg.c_str());
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
