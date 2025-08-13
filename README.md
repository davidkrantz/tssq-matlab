# stabilizing-quad
Implementation of [NAME].

## Setup

To initialize submodules, run
```
git submodule update --init --recursive
```

To compile the submodule `linequad`, `cd` to its `matlab` directory and run
```
make
```

## Usage

* Start Matlab in root directory.
* `init` initializes paths.
* `matlab/examples/` contains the code used to create the manuscript results.
* `matlab/test/` contains bits and pieces of code used when developing functionality.

## Code structure

* `external/` contains third-party packages in the form of Git submodules.
* `matlab/src/` contains all new routines needed for [NAME].

