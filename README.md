# Heat Sim - Parallel 2D blocked

## Introduction
Heat simulates the diffusion of heat in two dimensions. In each time step, a parallel, 2D blocked  Gauss-Seidel iterative solver is used to approximate a solution of the Poisson equation. The contininous problem is discretized with finite differences.

## Available versions and building instructions

The application includes diffrerent versions which can compiled to different binaries using the Makefile. They are:

  * **heat_seq**: Sequential version.
  * **heat_ompss**: Parallel version using OmpSs tasks.
  * **heat_mpi.pure**: Parallel version using MPI.
  * **heat_mpi.omp**: Parallel version using MPI + OmpSs tasking.
  * **heat_mpi.task**: Parallel version using MPI + OmpSs tasking and data-flows.
  * **heat_mpi.interop**: Parallel version using MPI + OmpSs tasking + Interoperability library. *See building instructions, step 1*.


  The simplest way to compile this package is:

  1. Stay in Heat root directory to recursively build all the versions.
     The Heat MPI + OmpSs tasks + Interoperability library version is
     compiled only if the environment variable `INTEROPERABILITY_SRC`
     is set.

  2. Type `make` to compile the selected benchmark's version(s).
     Optionally, you can use a different block size in each dimension
     (BSX and BSY for vertical and horizontal dimensions, respectively)
     when building the benchmark (by default 1024). Type
     `make BSX=MY_BLOCK_SIZE_X BSY=MY_BLOCK_SIZE_Y` in order to change
     this value. If you want the same value in each dimension, type
     `make BSX=MY_BLOCK_SIZE`.

  3. In addition, you can type 'make check' to check the correctness
     of the built versions. By default, the pure MPI version runs with
     4 processes and the hybrid versions run with 2 MPI processes and 2
     hardware threads for each process. You can change these
     parameters in 'scripts/run-tests.sh'.


OmpSs (OmpSs-2) is availible for download at www.pm.bsc.es. 
Please not that the interoperability library is not availible yet.


## OmpSs

OmpSs is a declarative, task-parallel programming model that support data-flow based task execution. It consists of a language specification, a source-to-source compiler for C, C++ and Fortran, and a runtime. The language defines a set of directives that allow a descriptive expression of tasks.  OmpSs allows the programmer to annotate task parameters with in, out and inout clauses that correspond to input, output or input-output access-type semantics of a parameter within that task. This information establishes a producer-consumer relationship between tasks, also called task dependency or data flow. With this information, the runtime is capable of automatic tasks scheduling that maintains correctness of code while relieving the programmer from implementing manual synchronization. In addition, the taskwait construct allows task synchronization and instructs the calling thread to wait on all previously created tasks. While this is similar to tasking in the recent specification of OpenMP, the OmpSs runtime implements a different execution model. 
In OmpSs, an application starts with a predefined set of execution resources and an explicit parallel region does not exist. This view avoids the exposure of threading to the programmer as well as the requirement to handle an additional scope, concretely that of a parallel region. At compile time, the OmpSs compiler processes pragma annotations and generates an intermediate code file. This file includes both user code as well as all required code for task generation, synchronization and error handling. In the final step of compilation, OmpSs invokes the native compiler to create a binary file. 

## Execution instructions

The binaries accept several options. The most relevant options are the size 
of the matrix with `-s` (default: 2048), and the number of timesteps with 
`-t` (default: 100). More options can be seen passing the `-h` option. An example hereof is:

```
$ mpiexec -n 4 -bind-to hwthread:16 heat_mpi.task.1024bs.exe -t 150 -s 8192
```

in which the application will perform 150 timesteps in 4 MPI processes with 16 
hardware threads in each process (used by the OmpSs runtime). The size of the
matrix in each dimension will be 8192 (8192^2 elements in total), this means
that each process will have 2048 * 8192 elements (16 blocks per process).

The provided configuation file (heat.conf) allows to configure the simulation and the topology of MPI processes. 
