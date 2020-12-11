---
name: 'OBS UFO Sprint: H(x) for an Instrument'
about: 'OBS UFO Sprint Template: H(x) for Instrument'
title: "[OBS UFO Sprint]: H(x) for Instrument_Name_Goes_Here"
labels: OBS1
assignees: dickdee, huishao-r, rhoneyager

---

# Description

For this instrument, we want to create and test the YAML for H(x) only.

The test: compare against the H(x) produced by GSI.

# Steps

- [x] When using this template, copy all of the text into a new issue. Set an appropriate name for the issue, set assignees, ZenHub release, epics, milestone, and complexity.
- [ ] Checkout the "feature/obs_sprint" branch (```git checkout feature/obs_sprint```) and the make a new branch based on this (```git checkout -b feature/my_obs_sprint_branch_here```).
- [ ] Ensure that we have the input file(s) for your test (in ```/data/users/hofx-dev/CodeSprintDec20/data/JediData/```).
- [ ] Change to the ```test/testinput/instrumentTests``` directory and use the ```template``` subdirectory as a base for your instrument. Copy this directory to a name based on your instrument (```cp -r template atms```).
- [ ] Start with an existing YAML file in ```test/testinput``` that is appropriate to your instrument. Rename it and customize it.
    - An example is given in  ```template/INSTRUMENT_gfs_HofX.yaml.in```.
    - [ ] The naming convention of your YAML should follow: ```instrument_gfs_HofX.yaml.in```.
    - [ ] As we are just testing H(x), strip out anything to do with bias correction and qc.
    - It is okay to hard-code file paths to the test data. 
    - [ ] Write the test.
        - Define tolerance in your yaml file, e.g.: 
        ```vector ref: GsiHofX```
        ```Tolerance: 1.e-6```  (use of the tolerance can be found in the slides for gnssro tests )
- [ ] Run the test:
    - Example script to submit job on S4: /data/users/hshao/JEDI_EMCtest/data/gsi_for_geoval/runctest/driver_test_operator.csh
    - ```/your/path/to/build/bin/test_ObsOperator.x   your_path_to_yaml_here```
- [ ] Examine the console output for success. vectorRef and RMSE.
    - If success, create a PR that points back to the "feature/obs_sprint" branch. Make sure that the PR points back to the ```feature/obs_sprint``` branch and not ```develop```.
        - The PR should set appropriate reviewers and should link against your H(x) issue.
    - If not, create another issue describing the problem and possible actions.
