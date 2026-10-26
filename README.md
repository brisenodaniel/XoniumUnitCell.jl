# XoniumUnitCell.jl

Julia simulations of Xonium (coupled fluxonium/transmon) unit-cell parameter sweeps and gates: spectroscopy and flux sweeps, parametric CZ gates, conditional phase and leakage, and flux-noise dephasing. Built on [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl).

## Notebooks at a glance

All notebooks live in `notebooks/` and use the Jupyter kernel **QuantumSims** (set up below).

| Notebook | What it does |
|---|---|
| `parametric_CZ.ipynb` | Reproduces the parametric CZ gate of Ma et al., [arXiv:2308.16040](https://arxiv.org/abs/2308.16040): coupling sweeps over flux bias, spectroscopy (Fig. 2a), and conditional phase from flux pulses (Fig. 2b). Multithreaded sweeps. Conventions and verified parameters are in `parametric_CZ_notes.md`. |
| `xonium_heavy_light.ipynb` | Heavy/light Xonium cell. It sweeps the qubit-subspace couplings against flux detuning, runs the spectroscopy and CZ conditional-phase sweeps (Fig. 2a/2b/2c), and estimates dephasing from flux modulation (Fig. 3). Multithreaded sweeps. |
| `2Q.ipynb` | Two-cell (`ncell`) exploration: a long-fluxonium blueprint cell, and ζ_zz cancellation sweeps for a nearly degenerate 2-cell and for a heavy/light 2-cell. |
| `circuit_quantization.ipynb` | Symbolic and numeric circuit quantization of coupled Xoniums: direct multi-path coupling (with and without diagonal capacitances) and lumped-resonator coupling (LRC). |
| `tai_chi_nodal_analysis.ipynb` | Circuit quantization by the method of nodes for the "tai chi" layout: direct multipath fluxonium-fluxonium coupling and the change to the difference/sum basis. |
| `kagan_notebooks/` | Collaborator notebooks (`tunable_coupler.ipynb`, `04_three-xoniums_flux-to_L.ipynb`), kept for reference. |

`output/` contains static HTML exports of `2Q`, `parametric_CZ` and `xonium_heavy_light`, so you can read the results without running anything.

`src/` holds the reusable code that the notebooks are being migrated onto (see [API.md](API.md) for its reference): circuit Hamiltonians (`circuits.jl`), gate and sweep utilities (`gate_utils.jl`), pulse envelopes (`pulses.jl`), noise and dephasing (`env_psd.jl`, `dephasing.jl`), and plotting (`plotting.jl`).

## Prerequisites

1. **Julia 1.12 or newer** (the kernel is named for the version you install with). The easiest way to install it is [juliaup](https://github.com/JuliaLang/juliaup).
2. **Git**, because two dependencies (`QuantumToolbox` and `QuantumSystems`) are fetched from forks on GitHub via the `[sources]` section of `Project.toml`.
3. **`ProjectRoot` installed in the base (global) Julia environment.** Every notebook runs `using ProjectRoot` *before* activating the project, so it has to be resolvable from the default environment. `make setup GET_PROJECTROOT=1` (or `make projectroot`) installs it for you; plain `make setup` only checks for it.
4. **make** (standard on Linux and macOS; on Windows use WSL, or run the Julia commands in the `Makefile` by hand).
5. **Jupyter** (JupyterLab, or VS Code with the Jupyter and Julia extensions) to open the notebooks. Installing the kernel only needs `IJulia`, which is a project dependency. It does not need Jupyter itself.

## Setup

Clone the repository and run, from its root:

```bash
make setup                      # ProjectRoot must already be in the base Julia environment
make setup GET_PROJECTROOT=1    # also install ProjectRoot there if it is missing
```

`make setup` runs `make kernel` and then `make deps`. By default (`GET_PROJECTROOT=0`) it leaves your base Julia environment alone: it checks that `ProjectRoot` is installed and stops with an error, before installing anything or saving any settings, if it is not. With `GET_PROJECTROOT=1` it runs `make projectroot` first. `THREADS` and `KERNEL_NAME` are passed on to `make kernel`, e.g. `make setup THREADS=8 KERNEL_NAME=MySims`.

Each step can also be run on its own (`make help` lists them):

| Command | What it does | Arguments (default) |
|---|---|---|
| `make setup` | Runs `make kernel` then `make deps`, after checking for `ProjectRoot` (or installing it, see `GET_PROJECTROOT`). | `GET_PROJECTROOT` (`0`; `0` or `1`), plus `THREADS` and `KERNEL_NAME` as for `make kernel` |
| `make projectroot` | Installs `ProjectRoot` in the base environment if it is missing. | none |
| `make deps` | Instantiates and precompiles this project's packages. The first run downloads and precompiles a lot, so it can take several minutes. | none |
| `make kernel` | Installs the Jupyter kernel **QuantumSims** (shown as `QuantumSims 1.12` on Julia 1.12) with `--project=@.`, so it uses the project wherever the repository is cloned. Values given on the command line are saved to `julia_settings.TOML`. | `THREADS` (`threads` in `julia_settings.TOML`, initially `auto`), `KERNEL_NAME` (`kernel_name` in `julia_settings.TOML`, initially `QuantumSims`) |
| `make kernel-remove` | Removes the kernel named by `KERNEL_NAME`. | `KERNEL_NAME` (as for `make kernel`) |
| `make session` | Starts a Julia REPL on this project. | `THREADS` (as for `make kernel`; for this session only, not saved) |
| `make test` | Runs the test suite. | none |

All commands also accept `JULIA` (`julia`), the Julia executable to use, and `SETTINGS` (`julia_settings.TOML`), the settings file.

### Threads and kernel name

The defaults for `make kernel` and `make session` live in `julia_settings.TOML` in the repository root:

```toml
kernel_name = "QuantumSims"
threads = "auto"
```

`threads` is a positive integer or `"auto"`. `auto` means `Sys.CPU_THREADS` (the number of logical CPUs on the machine), resolved when the kernel is installed or the session starts. The file is created with these defaults the first time you run `make kernel` or `make session`. It is git-ignored, since it is a per-machine setting.

Override the defaults on the command line:

```bash
make kernel                                   # settings from julia_settings.TOML
make kernel THREADS=8                         # exactly 8 threads
make kernel THREADS=auto                      # Sys.CPU_THREADS threads
make kernel THREADS=4 KERNEL_NAME=MySims      # a differently named kernel
```

When `make kernel` succeeds with `THREADS` or `KERNEL_NAME` given on the command line, the new values are written to `julia_settings.TOML` and become the defaults for later runs. A value that fails validation (for example `THREADS=0`) stops before anything is installed or saved. You can also edit the file by hand. Re-running `make kernel` with the same name updates that kernel, and a different name installs a separate one next to it.

Use fewer threads than cores if you want to keep the machine usable while a sweep runs, and keep memory in mind: each thread holds its own simulation state.

The notebooks pin BLAS to one thread (`LinearAlgebra.BLAS.set_num_threads(1)`) so Julia's threads don't oversubscribe the cores. The first cell of the sweep notebooks prints `Threads.nthreads()` and warns if only one thread is available.

### Scripts and the Julia REPL

A Jupyter kernel only applies to notebooks. A plain `julia` script or REPL does not use it, so give Julia the same settings the kernel uses (`--project` and `--threads`) on the command line. Run these from the repository root, or point `--project` at it.

Run a script (there is no Make target for this; `--threads=auto` follows Julia's own CPU detection):

```bash
julia --project=. --threads=auto my_script.jl
```

Start an interactive REPL with the project and the thread count from `julia_settings.TOML`:

```bash
make session
```

This is equivalent to `julia --project=. --threads=<N>`, with `<N>` resolved from the settings file (`auto` becomes `Sys.CPU_THREADS`). `make session THREADS=4` uses 4 threads for that session only. It does not change the file. Don't combine it with `make -j`.

Use a number instead of `auto` to fix the thread count. In Julia, `--threads=auto` follows the CPU affinity available to the process, so it can differ from the `Sys.CPU_THREADS` value that `make kernel` uses. Pass the number explicitly if you want them to match. Check the result with `Threads.nthreads()`.

To make the thread count the default for every Julia session in your shell, set it in the environment instead (the `--threads` flag overrides it):

```bash
export JULIA_NUM_THREADS=8
```

Scripts that use the sweep functions should pin BLAS to one thread, as the notebooks do, so the sweeps don't oversubscribe the cores:

```julia
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
```

`using XoniumUnitCell` works once you run with `--project=.`. If your script also uses `ProjectRoot`, as the notebooks do, it must be installed in the base environment (`make projectroot`, or `make setup GET_PROJECTROOT=1`).

To get a REPL that runs on the installed **kernel** itself, which needs Jupyter's console frontend, run:

```bash
jupyter kernelspec list                    # find the name, e.g. quantumsims-1.12
jupyter console --kernel=quantumsims-1.12
```

The kernel's name is the lowercased `KERNEL_NAME` plus the Julia version. This session has the kernel's thread count, but it starts in the directory you launched it from, so run it from the repository root (the kernel's `--project=@.` looks upward from there for `Project.toml`).

### VS Code

Select the kernel **QuantumSims 1.12** in the notebook's kernel picker. The notebook kernel gets its thread count from the kernelspec, not from the VS Code `julia.NumThreads` setting, which only affects the Julia REPL.

## Running the tests

```bash
make test
```

## Troubleshooting

- **`ArgumentError: Package ProjectRoot not found`** in the first cell: `ProjectRoot` is missing from the base environment. Run `make projectroot`.
- **Kernel `QuantumSims` not listed:** run `make kernel`, then restart Jupyter or VS Code. `jupyter kernelspec list` should show `quantumsims-<version>`.
- **Only 1 thread in the notebook:** the notebook is using a different kernel. Pick `QuantumSims` and restart it.
- **`Manifest.toml` is git-ignored**, so `Pkg.instantiate()` resolves versions fresh, and `QuantumToolbox` and `QuantumSystems` track the `main` branch of their forks. If a clean install breaks, ask the maintainer for their `Manifest.toml`.
- **Notebook save hangs in VS Code:** plot outputs are stored as both `image/png` and `text/html`, which can crash the VS Code serializer. Clear the notebook's outputs and save again.

## Conventions

Units are GHz and ns with ħ = 1. Frequencies are often angular (`2π·…`), so check each function's docstring. In the CZ setup the coupling `J` must be negative. See `notebooks/parametric_CZ_notes.md` before touching Hamiltonian, sign or phase code.

## License

See [LICENSE](LICENSE).
