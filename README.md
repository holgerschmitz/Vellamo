# Vellamo

Vellamo is a simple fluid code based on the Schnek library and making use of the Huerto algorithm repository.

## Installation

### Prerequisites

Vellamo has currently been tested on Linux only.

You need
* C++ compiler
* MPI, both mpich or Open MPI have been tested
* HDF5, build with MPI support to allow parallel file output
* Schnek build with MPI and HDF5

### Obtaining and Building the Code

There are currently no releases. Obtain the code directly from the GitHub repository.

```bash
git clone --recurse-submodules https://github.com/holgerschmitz/Vellamo.git
```

The `--recurse-submodules` option is important so that you also get the correct version of Huerto.

You may want to look at the `Makefile` before building the code. The variable `HDFBASE` should be set to the location of your HDF5 installation. Once you are happy with everything, simply run

```bash
make
```

This will build three versions of the code and place them into the `bin/` folder.

### Running the Simulation

Run the simulation by executing `vellamo1d`, `vellamo2d`, or `vellamo3d` inside a folder with a `vellamo.setup` file. Look at the `examples/` folder for some examples.

More documentation coming soon.