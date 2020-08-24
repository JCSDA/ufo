/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_RECURSIVESPLITTER_H_
#define TEST_UFO_RECURSIVESPLITTER_H_

#include <iomanip>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "eckit/testing/Test.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {
namespace test {

typedef std::set<size_t> Group;
typedef std::set<Group> Groups;

Groups getGroups(const RecursiveSplitter &splitter) {
  Groups groups;
  for (const auto &g : splitter.groups()) {
    Group group;
    for (auto index : g)
      group.insert(index);
    groups.insert(group);
  }
  return groups;
}

Groups getMultiElementGroups(const RecursiveSplitter &splitter) {
  Groups groups;
  for (const auto &g : splitter.multiElementGroups()) {
    Group group;
    for (auto index : g)
      group.insert(index);
    groups.insert(group);
  }
  return groups;
}

CASE("ufo/RecursiveSplitter/ZeroIds") {
  RecursiveSplitter splitter(0);
  {
    Groups expectedGroups;
    Groups groups = getGroups(splitter);
    EXPECT_EQUAL(groups, expectedGroups);

    Groups expectedMultiElementGroups;
    Groups multiElementGroups = getMultiElementGroups(splitter);
    EXPECT_EQUAL(multiElementGroups, expectedMultiElementGroups);
  }

  {
    std::vector<size_t> categories;
    splitter.groupBy(categories);

    Groups expectedGroups;
    Groups groups = getGroups(splitter);
    EXPECT_EQUAL(groups, expectedGroups);

    Groups expectedMultiElementGroups;
    Groups multiElementGroups = getMultiElementGroups(splitter);
    EXPECT_EQUAL(multiElementGroups, expectedMultiElementGroups);
  }
}

CASE("ufo/RecursiveSplitter/OneId") {
  RecursiveSplitter splitter(1);
  {
    Groups expectedGroups{{0}};
    Groups groups = getGroups(splitter);
    EXPECT_EQUAL(groups, expectedGroups);

    Groups expectedMultiElementGroups;
    Groups multiElementGroups = getMultiElementGroups(splitter);
    EXPECT_EQUAL(multiElementGroups, expectedMultiElementGroups);
  }

  {
    std::vector<size_t> categories{1};
    splitter.groupBy(categories);

    Groups expectedGroups{{0}};
    Groups groups = getGroups(splitter);
    EXPECT_EQUAL(groups, expectedGroups);

    Groups expectedMultiElementGroups;
    Groups multiElementGroups = getMultiElementGroups(splitter);
    EXPECT_EQUAL(multiElementGroups, expectedMultiElementGroups);
  }
}

CASE("ufo/RecursiveSplitter/TenIds") {
  RecursiveSplitter splitter(10);
  {
    Groups expectedGroups{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}};
    Groups groups = getGroups(splitter);
    EXPECT_EQUAL(groups, expectedGroups);

    Groups expectedMultiElementGroups{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}};
    Groups multiElementGroups = getMultiElementGroups(splitter);
    EXPECT_EQUAL(multiElementGroups, expectedMultiElementGroups);
  }

  {
    std::vector<size_t> categories{1, 2, 1, 2, 1, 2, 1, 2, 1, 2};
    splitter.groupBy(categories);

    Groups expectedGroups{{0, 2, 4, 6, 8}, {1, 3, 5, 7, 9}};
    Groups groups = getGroups(splitter);
    EXPECT_EQUAL(groups, expectedGroups);

    Groups expectedMultiElementGroups{{0, 2, 4, 6, 8}, {1, 3, 5, 7, 9}};
    Groups multiElementGroups = getMultiElementGroups(splitter);
    EXPECT_EQUAL(multiElementGroups, expectedMultiElementGroups);
  }

  {
    // The first and last category in the first group have only a single element
    std::vector<size_t> categories{3, 2, 1, 2, 3, 2, 4, 2, 3, 2};
    splitter.groupBy(categories);

    Groups expectedGroups{{2}, {0, 4, 8}, {6}, {1, 3, 5, 7, 9}};
    Groups groups = getGroups(splitter);
    EXPECT_EQUAL(groups, expectedGroups);

    Groups expectedMultiElementGroups{{0, 4, 8}, {1, 3, 5, 7, 9}};
    Groups multiElementGroups = getMultiElementGroups(splitter);
    EXPECT_EQUAL(multiElementGroups, expectedMultiElementGroups);
  }

  {
    // Two multiple-element categories in the second group.
    // Categories shared by multiple groups.
    std::vector<size_t> categories{1, 3, 1, 3, 1, 2, 1, 1, 1, 1};
    splitter.groupBy(categories);

    Groups expectedGroups{{2}, {0, 4, 8}, {6}, {1, 3}, {5}, {7, 9}};
    Groups groups = getGroups(splitter);
    EXPECT_EQUAL(groups, expectedGroups);

    Groups expectedMultiElementGroups{{0, 4, 8}, {7, 9}, {1, 3}};
    Groups multiElementGroups = getMultiElementGroups(splitter);
    EXPECT_EQUAL(multiElementGroups, expectedMultiElementGroups);
  }
}

CASE("ufo/RecursiveSplitter/IntCategories") {
  RecursiveSplitter splitter(10);

  {
    std::vector<int> categories{1, 2, 1, 2, 1, 2, 1, 2, 1, 2};
    splitter.groupBy(categories);

    Groups expectedGroups{{0, 2, 4, 6, 8}, {1, 3, 5, 7, 9}};
    Groups groups = getGroups(splitter);
    EXPECT_EQUAL(groups, expectedGroups);

    Groups expectedMultiElementGroups{{0, 2, 4, 6, 8}, {1, 3, 5, 7, 9}};
    Groups multiElementGroups = getMultiElementGroups(splitter);
    EXPECT_EQUAL(multiElementGroups, expectedMultiElementGroups);
  }
}

CASE("ufo/RecursiveSplitter/StringCategories") {
  RecursiveSplitter splitter(10);

  {
    std::vector<std::string> categories{"abc", "def", "abc", "def",
                                        "abc", "def", "abc", "def",
                                        "abc", "def"};
    splitter.groupBy(categories);

    Groups expectedGroups{{0, 2, 4, 6, 8}, {1, 3, 5, 7, 9}};
    Groups groups = getGroups(splitter);
    EXPECT_EQUAL(groups, expectedGroups);

    Groups expectedMultiElementGroups{{0, 2, 4, 6, 8}, {1, 3, 5, 7, 9}};
    Groups multiElementGroups = getMultiElementGroups(splitter);
    EXPECT_EQUAL(multiElementGroups, expectedMultiElementGroups);
  }
}

class RecursiveSplitter : public oops::Test {
 public:
  RecursiveSplitter() {}

 private:
  std::string testid() const override {return "ufo::test::RecursiveSplitter";}

  void register_tests() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_RECURSIVESPLITTER_H_
