# Power Dominating Set
This is the implementation of “An Efficient Algorithm for Power Dominating Set” [1].

## Prerequesites
This library requires the following dependencies to be installed before compiling:
* boost
* gurobi and valid `GUROBI_HOME`
* range-v3
* vcpkg (with library tinyxml2)

## Building
This project is built using cmake.
For a building execute the following commands:
```
git clone $GIT/pds-code.git
cd pds-code
cmake -B build -S . --preset=default (or, sudo -E cmake -B build -S . --preset=default)
cmake --build build (or, sudo cmake --build build)
```

## Running
Setup:
```
cp experiments/* build
cd build
ln -s ../inputs
# solver performance
./run_solver_experiment.sh
# upper and lower bounds
./run_bounds_experiment.sh
# reduction rule effectiveness
./run_reduction_experiment.sh
```

### Evaluation
The original data and evaluation scripts for the arXiv version can be found in `evaluation.tar.xz`

[1]: https://arxiv.org/abs/2306.09870
