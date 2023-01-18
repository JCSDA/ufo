---
name: 'Sprint: UFO Naming Conventions'
about: 'Use this template to create a PR issue for the sprint'
title: 'SPRINT: Update UFO (what portion?) to follow Data Conventions'
labels: 'Sprint' 
assignees: ''

---

Within title (above), enter the part of UFO you are changing to the data conventions.


## Description

UFO (what portion) update to follow naming convention

Provide a detailed description if necessary

### Issue(s) addressed

Link the issues to be closed with this PR
- fixes #<issue_number>

If this does not close an existing issue, then [be sure to add Label, Estimate, Assignees, and Epic] as follows:
- Labels: sprint
- Estimate: An estimate of how long this issue will take to resolve. 1 = half a day. 2 = full day. 3 = two days, 5 = three days.
- Epics: "OBSPROC-47: Implement and maintain IODA conventions for UFO via documentation, code repository and tables" 
- Assignees: Anyone **working** on this issue. Anyone who is not working directly on the issue, but
  who should receive notifications, should be listed below.


## Dependencies

If there are PRs that need to be merged before or along with this one, please add "waiting for another PR" label and list the dependencies in the other repositories (example below). Note that the branches in the other repositories should have matching names for the automated tests to pass.
Waiting on the following PRs:
- [ ] waiting on JCSDA/eckit/pull/<pr_number>
- [ ] waiting on JCSDA/atlas/pull/<pr_number>

