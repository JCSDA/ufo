/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_INTEGRALIMAGE_H_
#define TEST_UFO_INTEGRALIMAGE_H_

#include <Eigen/Dense>

#include <string>
#include <tuple>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/TestEnvironment.h"

#include "ufo/utils/IntegralImage.h"

namespace ufo {
namespace test {

namespace integralimage {
  // Set up test and reference data sets.
  std::tuple<Eigen::ArrayXXf,
    Eigen::ArrayXXi,
    std::vector<std::vector<float>>,
    std::vector<std::vector<int>>>
    createTestAndReferenceArrays(const int numrows, const int numcols) {
    Eigen::ArrayXXf arrayFloat(numrows, numcols);
    Eigen::ArrayXXi arrayInt(numrows, numcols);
    std::vector<std::vector<float>> referenceIntegralImageFloat;
    std::vector<std::vector<int>> referenceIntegralImageInt;

    // Fill array from which an integral image will be computed.
    for (int i = 0; i < numrows; ++i) {
      for (int j = 0; j < numcols; ++j) {
        arrayFloat(i, j) = static_cast<float>(i * j);
        arrayInt(i, j) = i * j;
      }
    }

    // Create reference integral images.
    for (int i = 0; i < numrows; ++i) {
      referenceIntegralImageFloat.push_back(std::vector<float>(numcols, 0));
      referenceIntegralImageInt.push_back(std::vector<int>(numcols, 0));
    }
    for (int i = 0; i < numrows; ++i) {
      for (int j = 0; j < numcols; ++j) {
        referenceIntegralImageFloat[i][j] = static_cast<float>(i * j);
        referenceIntegralImageInt[i][j] = i * j;
      }
    }

    // Integral image calculation.
    for (std::size_t i = 1; i < numrows; ++i) {
      referenceIntegralImageFloat[i][0] += referenceIntegralImageFloat[i - 1][0];
      referenceIntegralImageInt[i][0] += referenceIntegralImageInt[i - 1][0];
    }
    for (std::size_t j = 1; j < numcols; ++j) {
      referenceIntegralImageFloat[0][j] += referenceIntegralImageFloat[0][j - 1];
      referenceIntegralImageInt[0][j] += referenceIntegralImageInt[0][j - 1];
    }
    for (std::size_t i = 1; i < numrows; ++i) {
      for (std::size_t j = 1; j < numcols; ++j) {
        referenceIntegralImageFloat[i][j] += referenceIntegralImageFloat[i][j - 1] +
          referenceIntegralImageFloat[i - 1][j] -
          referenceIntegralImageFloat[i - 1][j - 1];
        referenceIntegralImageInt[i][j] += referenceIntegralImageInt[i][j - 1] +
          referenceIntegralImageInt[i - 1][j] -
          referenceIntegralImageInt[i - 1][j - 1];
      }
    }

    return {arrayFloat, arrayInt, referenceIntegralImageFloat, referenceIntegralImageInt};
  }
}  // namespace integralimage

CASE("ufo/IntegralImage/createIntegralImage") {
  // Initialise arrays.
  const int numrows = 5;
  const int numcols = 6;
  auto[arrayFloat, arrayInt, referenceIntegralImageFloat, referenceIntegralImageInt] =
    ufo::test::integralimage::createTestAndReferenceArrays(numrows, numcols);

  // Compute integral image for each array.
  ufo::createIntegralImage(arrayFloat);
  ufo::createIntegralImage(arrayInt);

  // Ensure the integralImage has worked as expected.
  for (std::size_t i = 1; i < numrows; ++i) {
    for (std::size_t j = 1; j < numcols; ++j) {
      EXPECT(arrayFloat(i, j) == referenceIntegralImageFloat[i][j]);
      EXPECT(arrayInt(i, j) == referenceIntegralImageInt[i][j]);
    }
  }
}

CASE("ufo/IntegralImage/sumIntegralImagePatch") {
  // Initialise arrays.
  const int numrows = 5;
  const int numcols = 6;
  auto[arrayFloat, arrayInt, referenceIntegralImageFloat, referenceIntegralImageInt] =
    ufo::test::integralimage::createTestAndReferenceArrays(numrows, numcols);

  // Compute integral image for each array.
  ufo::createIntegralImage(arrayFloat);
  ufo::createIntegralImage(arrayInt);

  // Sum integral images between a variety of limits, both with and without bounds checking.
  // Do for a variety of step sizes.
  for (int istep = 0; istep < 3; ++istep) {
    for (int jstep = 0; jstep < 3; ++jstep) {
      for (int i = istep + 1; i < numrows - istep; ++i) {
        for (int j = jstep + 1; j < numcols - jstep; ++j) {
          const float sumFloat = ufo::sumIntegralImagePatch(arrayFloat, i, j, istep, jstep);
          const int sumInt = ufo::sumIntegralImagePatch(arrayInt, i, j, istep, jstep);
          const float sumFloatWithBoundsCheck =
            ufo::sumIntegralImagePatchWithBoundsCheck(arrayFloat, i, j, istep, jstep);
          const int sumIntWithBoundsCheck =
            ufo::sumIntegralImagePatch(arrayInt, i, j, istep, jstep);
          const float expectedSumFloat =
            referenceIntegralImageFloat[i + istep][j + jstep] -
            referenceIntegralImageFloat[i + istep][j - jstep - 1] -
            referenceIntegralImageFloat[i - istep - 1][j + jstep] +
            referenceIntegralImageFloat[i - istep - 1][j - jstep - 1];
          const float expectedSumInt =
            referenceIntegralImageInt[i + istep][j + jstep] -
            referenceIntegralImageInt[i + istep][j - jstep - 1] -
            referenceIntegralImageInt[i - istep - 1][j + jstep] +
            referenceIntegralImageInt[i - istep - 1][j - jstep - 1];
          EXPECT(sumFloat == expectedSumFloat);
          EXPECT(sumFloatWithBoundsCheck == expectedSumFloat);
          EXPECT(sumInt == expectedSumInt);
          EXPECT(sumIntWithBoundsCheck == expectedSumInt);
        }
      }
    }
  }

  // Deliberately throw bounds checking exceptions.
  EXPECT_THROWS_MSG(ufo::sumIntegralImagePatchWithBoundsCheck(arrayFloat, 0, 0, 100, 100),
                    "One or more bounds exceeded");
  EXPECT_THROWS_MSG(ufo::sumIntegralImagePatchWithBoundsCheck(arrayFloat, 100, 100, 0, 0),
                    "One or more bounds exceeded");
}

class IntegralImage : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::IntegralImage";}

  void register_tests() const override {}

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_INTEGRALIMAGE_H_
