# ROOT-Sim core 3.0.0

*Brought to you by the [High Performance and Dependable Computing Systems (HPDCS)](https://hpdcs.github.io/)
Research Group*

[![Build Status](https://github.com/ROOT-Sim/core/workflows/ROOT-Sim%20core%20CI/badge.svg)](https://github.com/ROOT-Sim/core/actions)
[![codecov.io](https://codecov.io/gh/ROOT-Sim/branch/master/graphs/badge.svg)](https://codecov.io/gh/ROOT-Sim/core)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/7519f016f3d942b9b12c6ed03ae4ecf8)](https://www.codacy.com/gh/ROOT-Sim/core/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ROOT-Sim/core&amp;utm_campaign=Badge_Grade)
[![doc coverage](https://img.shields.io/endpoint?url=https%3A%2F%2Froot-sim.github.io%2Fcore%2Fdocs%2Fcoverage%2Fmaster.json)](https://root-sim.github.io/core/docs/)
[![GitHub issues](https://img.shields.io/github/issues/ROOT-Sim/core)](https://github.com/ROOT-Sim/core/issues)
[![GitHub](https://img.shields.io/github/license/ROOT-Sim/core)](https://github.com/ROOT-Sim/core/blob/master/LICENSES/GPL-3.0-only.txt)
[![REUSE Compliance Check](https://github.com/ROOT-Sim/core/actions/workflows/reuse_check.yml/badge.svg)](https://github.com/ROOT-Sim/core/actions/workflows/reuse_check.yml)

----------------------------------------------------------------------------------------

## The ROme OpTimistic Simulator

The ROme OpTimistic Simulator is an open source, distributed and parallel simulation framework developed using C/POSIX
technology. It transparently supports all the mechanisms associated with parallelization and distribution of workload
across the nodes (e.g., mapping of simulation objects on different kernel instances) and optimistic synchronization (
e.g., state recoverability). Distributed simulations rely on MPI3. In particular, global synchronization across the
different nodes relies on asynchronous MPI primitives, for increased efficiency.

The programming model supported by ROOT-Sim allows the simulation model developer to use a simple application-callback
function named `ProcessEvent()` as the event handler, whose parameters determine which simulation object is currently
taking control for processing its next event, and where the state of this object is located in memory. An object is a
data structure, whose state can be scattered on dynamically allocated memory chunks, hence the memory address passed to
the callback locates a top level data structure implementing the object state-layout.

ROOT-Sim's development started as a research project late back in 1987, and is currently maintained by the High
Performance and Dependable Computing Systems group, research group of the University of Rome "Tor Vergata".

## ROOT-Sim Core

This repository keeps the sources of the ROOT-Sim core: this is the fundamental library that implements the largest
part of the simulation algorithms used in the simulation framework.

The core can be built and used as a stand-alone low-level library writing C code, or it can be used within other
projects, such as [cROOT-Sim](https://gihub.com/ROOT-Sim/cROOT-Sim), i.e. the C/C++ version of the simulation library.

## Dependencies and platforms

The core successfully compiles on x86 and ARM architectures, using either GCC or Clang compilers, on Linux, Windows,
and macOS.
A compiler supporting the C11 standard is required, such as GCC 8 or later. MSVC on Windows does not properly implement
the full C11 standard (e.g., `stdatomic.h` is not provided), and cannot be therefore used to build the project.

MPI is a mandatory dependency of the project, used to support simulations run on distributed systems.
The core is continuously tested against the following MPI implementations:
*   OpenMPI
*   MPICH
*   Microsoft MPI

Any of the three is required to build the project. A full MPI3 implementation, supporting multithreading, is necessary.

## Building

To build the project, run:

```bash
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=gcc ..
make
```

You can specify a different compiler using the `-DCMAKE_C_COMPILER=` flag.
To run the test suite (which includes a correctness test), run in the `build` folder:

```bash
ctest
```

## Compiling and running a model

The ROOT-Sim core is not expected to be used directly to run models (see, for example,
[cROOT-Sim](https://gihub.com/ROOT-Sim/cROOT-Sim)). Nevertheless, an implementation of a "low-level" model
can be located in `test\integration`.
The test can be compiled using the standard `mpicc` compiler, linking against `librscore` and launching either
locally or using `mpiexec` to run on multiple nodes.
