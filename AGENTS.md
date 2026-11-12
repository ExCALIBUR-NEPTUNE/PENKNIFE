# Working on PENKNIFE

PENKNIFE (Plasma Edge Numerics using Kinetic Neutrals Integrated in Finite
Elements) couples Nektar++ fluid solvers with NESO particle systems and
VANTAGE reactions. The application uses C++17, CMake, MPI, and SYCL.

## Workflow

- Check the worktree and any applicable nested `AGENTS.md` before working.
  Preserve existing user changes, including changes inside submodules.
- For source changes, prepare an unapplied proposed diff, explain its behavior
  and validation plan, and wait for approval before applying it. Explicit
  authorization to edit directly overrides this default. Documentation requested
  by the user may be written directly.
- Keep changes focused. Reuse nearby solver, operator, and storage patterns;
  avoid extra wrappers, descriptors, masks, buffers, or virtual APIs unless the
  task requires them.
- Inspect the local dependency APIs rather than assuming upstream Nektar++ or
  NESO interfaces match this checkout. Change dependencies only when required by
  the task, and report those changes separately.
- After approved edits, run the relevant checks and report what changed, what
  passed, and any validation limits. Do not claim GPU or MPI coverage from a
  serial CPU run.

## Repository map

- `src/PENKNIFE.cpp`: executable entry point.
- `src/EquationSystems/`: plasma systems, including ReducedBraginskii,
  ElectrostaticTurbulence, and anisotropic diffusive fields.
- `src/ParticleSystems/`: particles, reactions, and fluid coupling.
- `src/Closures/`: transport closures and heat flux calculations.
- `src/FluxAdvDiff/`, `src/FluxVolDiff/`, `src/FluxVolAdvDiff/`, `src/FluxTraceAdvDiv/`:
  volume and trace operators, kernels, and explicit factory registrations.
- `src/BndCond/`, `src/MagneticField/`, `src/Pointwise/`: supporting
  boundary and field operators.
- `tests/EquationSystems/`, `tests/TraceFlux/`, `tests/Bohm/`: targeted checks.
- `examples/`: input cases; `runs/`: working copies and simulation output.
- `NESO/`, `NESO/nektar/`, `NESO/neso-particles/`, `VANTAGE-Reactions/`:
  dependency checkouts; inspect their own guidance when working there.
- `repo/packages/penknife/package.py`, `spack.yaml`, `scripts/`: packaging and
  environment/run helpers.
- `docs/sphinx/source/`: user, developer, and theory documentation.

## Build and run

The documented setup uses Spack. Once Spack is available, initialize missing
submodules on a fresh checkout and activate the repository environment:

```bash
git submodule update --init --recursive
source ./scripts/activate.sh
spack install -j <jobs>
```

For an existing checkout, inspect its configured build and dependency paths
before rebuilding. Reuse a compatible compiler and dependency environment;
avoid rebuilding the full dependency stack for a small solver change. The
executable target is `penknife`:

```bash
cmake --build <build-dir> --target penknife -j <jobs>
```

Run examples through the helper with an explicit build directory:

```bash
./scripts/run_eg.sh <EquationSystem> <Mesh> <2D/3D> <Example> -n <ranks> -b <build-dir>
```

The helper copies the case into `runs/` and prompts before replacing an
existing run. Keep generated checkpoints, meshes, plots, and virtual
environments out of source changes.

## Validation

Run commands in the configured dependency environment. Fresh CMake builds need
package discovery configured, for example through `CMAKE_PREFIX_PATH` or package
`*_DIR` settings from a compatible existing build.

Equation-system checks belong to the main project and are off by default:

```bash
cmake -S . -B <build-dir> -DPENKNIFE_BUILD_TESTS=ON
cmake --build <build-dir> --target check_equation_systems -j <jobs>
ctest --test-dir <build-dir> -R '^equations_' --output-on-failure
```

TraceFlux and Bohm are standalone CMake projects; enabling the main project's
tests does not enable them:

```bash
cmake -S tests/TraceFlux -B <trace-build-dir>
cmake --build <trace-build-dir> -j <jobs>
ctest --test-dir <trace-build-dir> --output-on-failure

cmake -S tests/Bohm -B <bohm-build-dir>
cmake --build <bohm-build-dir> -j <jobs>
ctest --test-dir <bohm-build-dir> --output-on-failure
```

- Run TraceFlux for trace/volume diffusion changes and Bohm for sheath or
  boundary-routing changes. TraceFlux includes a two-rank test when an MPI
  launcher is available. AVX2 checks require a capable CPU and
  `-DPENKNIFE_TRACE_TEST_AVX2=ON`.
- Match tests to the behavior changed, including relevant projection, species,
  dimension, and execution backend. Add a regression for new behavior when
  existing checks do not cover it.
- Use `git diff --check` and inspect the final diff. Documentation-only changes
  do not require solver builds.

## Implementation conventions

- Follow `.clang-format`: four-space indentation and Allman braces. Format
  touched code without reformatting unrelated files. Follow nearby member
  naming (`m_` or `this->`).
- Access redesigned Nektar `Field` blocks directly. Preserve padding,
  component strides, interleave widths, and execution-space accessibility;
  do not flatten blocks into legacy layouts. Use direct non-owning coefficient
  pointers where the existing kernels support them.
- Keep coupled anisotropic penalties through `FluxVolDiffPoint` using the
  penalty gradient `tau * n * [U]`. An independent scalar penalty per conserved
  component changes the transport law.
- In fused trace diffusion, `workspace.state` holds the mean state,
  `workspace.gradient` holds the penalty gradient plus the averaged physical
  gradient when required, and `workspace.flux` holds the Cartesian diffusive
  flux before contraction with the normal. Preserve the existing sign and
  boundary-weight conventions.
- Follow the redesigned compressible-flow solver's boundary storage and gather
  patterns. Keep Bohm/Riemann trace classification in setup and preserve the
  adjacent boundary advection and diffusion calls. Diffusion exterior states
  affect `m_gloT1Diff`; advection uses `m_gloT1`.
- Preserve explicit operator factory registration across supported execution
  spaces and data types when extending operators.
- Consult `docs/sphinx/source/developer_guide/contributing.rst` for contribution
  conventions; it directs non-trivial merge requests to `dev`.
