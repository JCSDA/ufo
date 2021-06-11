/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "../ufo/ParallelPoissonDiskThinning.h"
#include "oops/runs/Run.h"

int main(int argc, char ** argv) {
  oops::Run run(argc, argv);
  ufo::test::ParallelPoissonDiskThinning tests;
  return run.execute(tests);
}
