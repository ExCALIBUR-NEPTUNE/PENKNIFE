# PENKNIFE
Plasma Edge Numerics using Kinetic Neutrals Integrated in Finite Elements

## Building
The easiest way to build the app is to use the [Spack package manager](https://spack.readthedocs.io/en/latest/index.html).
The build process has been tested with **Spack v1.2.0.dev0**; later versions may also work, earlier than **v1.0.0** is not recommended.

### Install Spack via the [official instructions](https://spack.readthedocs.io/en/latest/getting_started.html#installation).
On Debian-based systems (e.g. Ubuntu) the following should work:

```bash
git clone --depth=2 https://github.com/spack/spack.git

# Initialise spack
source $SPACK_ROOT/share/spack/setup-env.sh

# Optionally modify .bashrc such that spack is always initialised in new shells (this only has to be done once)
echo 'source $SPACK_ROOT/share/spack/setup-env.sh' >> $HOME/.bashrc
```

### Install intel compilers
```bash
spack install intel-oneapi-compilers
spack compiler add `spack location -i intel-oneapi-compilers`/compiler/latest/bin
```

Or if they are already installed elsewhere and version < 2024, add them
```bash
spack compiler add `spack location -i intel-oneapi-compilers`/compiler/latest/linux/bin/intel64
spack compiler add `spack location -i intel-oneapi-compilers`/compiler/latest/linux/bin
```
N.B. There is a known problem/conflict with the mkl libraries if the compilers come from the oneapi toolkit - it is recommended not to use this

### Build and clone PENKNIFE:
```bash
cd ..
git clone https://github.com/ExCALIBUR-NEPTUNE/PENKNIFE.git
cd PENKNIFE
git submodule update --init --recursive
spack env activate ./spack -p
spack install
```
N.B.
```bash
spack install -j <nproc>
```
can be used to speed up compilation

Note that some of the dependencies (particularly nektar++) can take some time to install and have high memory usage.

## Running examples
Configuration files for various different examples can be found in the `./examples` directory.
An easy way to run the examples is (from the top level of the repository):

```bash
# Load the nektar spack module, which also loads mpi
spack load nektar
# Optionally source the tab autocomplete script
source ./scripts/run_eg-completion.sh
# Set up and run an example via the helper script
./scripts/run_eg.sh [EquationSystem] [Mesh] [Dimension] [Example] <-n num_MPI>
# e.g. ./scripts/run_eg.sh SingleField MASTU 2D CG -n 16
# If the autocomplete script has been sourced the tab key shows the available examples at each level

```

- This will copy `./examples/[EquationSystem]/[Mesh]/[Dimension]/[Example]` to `./runs/[EquationSystem]/[Mesh]/[Dimension]/[Example]` and run the solver in that directory with num_MPP MPI processes
- The solver executable is assumed to be in the most recently modified `spack-build*` directory. To choose a different build location, use `<-b build_dir_path>`.

To change the number of openMP threads used by each process, use
```bash
# Run an example with OMP and MPI
OMP_NUM_THREADS=[nthreads] ./scripts/run_eg.sh [EquationSystem] [Mesh] [Dimension] [Example]
```

## Postprocessing
- Particle output is generated in an hdf5 file.
- Fluid output is generated as nektar++ checkpoint files. One way to visualise is them is to convert to .vtu using Nektar's `FieldConvert` tool and then use Paraview:
```bash
# Convert to vtu using FieldConvert (requires xml files as args)
cd runs/[EquationSystem]/[Mesh]/[Dimension]/[Example]
spack load nektar # if not already loaded
FieldConvert [config_xml] [mesh_xml] [chk_name] [vtu_name]
# e.g. FieldConvert single_field.xml mastu.xml single_field_100.chk single_field.vtu
paraview single_field.vtu
```

To convert multiple files at once, the following command can be used:
```bash
for i in {0..1000}; do FieldConvert -f single_field.xml mastu.xml single_field_$i.chk single_field_$i.vtu; done;
```
