/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "UfoTrait.h"
#include "oops/runs/Run.h"
#include "test/interface/ObsVector.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::ObsVector<ufo::UfoTrait> tests;
  run.execute(tests);
  return 0;
};

