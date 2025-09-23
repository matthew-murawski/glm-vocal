# glm-vocal

glm-vocal is a matlab toolkit for fitting poisson generalized linear models to single-unit spike counts during natural vocal interactions, packaging timebase construction, event cleanup, design-matrix assembly, and evaluation utilities so session analyses stay reproducible.

## Running Tests

Use the provided script to set up paths, then run the test suite:

```
matlab -batch "addpath('scripts'); startup; results = runtests('tests','IncludeSubfolders',true); disp(results);"
```
