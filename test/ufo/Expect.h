/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_EXPECT_H_
#define TEST_UFO_EXPECT_H_

#include <sstream>

// IMPORTANT: The user should include "eckit/testing/Test.h", if necessary after
// defining ECKIT_TESTING_SELF_REGISTER_CASES appropriately.

// Provides more informative output on failure than the raw EXPECT() macro.
#define EXPECT_EQUAL(expr, expected) \
    do { \
        if (!((expr) == (expected))) { \
            std::stringstream str; \
            str << ("EXPECT condition '" #expr " == " #expected "' failed. ") \
                << "(Received: " << expr << "; expected: " << expected << ")"; \
            throw eckit::testing::TestException(str.str(), Here()); \
        } \
    } while (false)

#endif  // TEST_UFO_EXPECT_H_
