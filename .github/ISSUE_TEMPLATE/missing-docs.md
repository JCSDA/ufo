---
name: 'Missing Documentation'
about: 'Use this issue to note missing documentation on ReadTheDocs and/or Doxygen'
title: "[Documentation]: Give a descriptive title here"
labels: documentation
assignees: 'huishao-r', 'rhoneyager', 'mmiesch'

---

## Description

> Edit the text of this pull request to describe what is missing. You should include the file name or the class name.



## ZenHub organization

Please make sure that you have properly set:
- Labels: documentation
- Milestone: OBS - current month
- Estimate: usually 0.5 - 1
- Epics: Select a relevant epic. Here are a few hints.
    - JEDI1.10: For **generic** UFO components and bias correction. Generic means a core UFO component, such as a base class
        or a component of JEDI that more than one agency will use.
    - OBS3.1.2: For **specific** ObsOperators, Filters, and ObsFunctions.
- Assignees: Generally, the person who created the code is responsible for writing the documentation.
    The "git blame" feature may be quite helpful in determining this. @huishao-r, @mmiesch, and @rhoneyager are
    automatically added to help track these issues.


