# Setup tasks for XoniumUnitCell.jl. Run `make help` for an overview.

JULIA    ?= julia
# Default kernel settings live in this file (created on first use). Values
# given on the command line, e.g. `make kernel THREADS=8`, override the file
# and are written back to it by `make kernel`.
SETTINGS ?= julia_settings.TOML
SETTINGS_JL = $(JULIA) --startup-file=no scripts/settings.jl $(SETTINGS)

# Number of Julia threads: a positive integer, or "auto", which is replaced by
# Sys.CPU_THREADS, the number of logical CPUs on this machine.
THREADS     ?= $(shell $(SETTINGS_JL) get threads)
# Kernel display name; the installed spec is "<lowercased name>-<julia version>".
KERNEL_NAME ?= $(shell $(SETTINGS_JL) get kernel_name)
NTHREADS     = $(shell $(SETTINGS_JL) resolve $(THREADS))

.PHONY: help setup projectroot deps kernel kernel-remove session test

help:
	@echo "Targets:"
	@echo "  make setup         projectroot + deps + kernel"
	@echo "  make projectroot   install ProjectRoot in the base Julia environment"
	@echo "  make deps          instantiate and precompile this project's packages"
	@echo "  make kernel        install the Jupyter kernel and save THREADS/KERNEL_NAME to $(SETTINGS)"
	@echo "  make kernel-remove remove the Jupyter kernel named KERNEL_NAME"
	@echo "  make session       start a Julia REPL on this project with THREADS threads"
	@echo "  make test          run the test suite"
	@echo ""
	@echo "Variables (defaults from $(SETTINGS)):"
	@echo "  THREADS=$(THREADS)        integer, or 'auto' = Sys.CPU_THREADS"
	@echo "  KERNEL_NAME=$(KERNEL_NAME)"
	@echo "  JULIA=$(JULIA)"
	@echo "Example: make kernel THREADS=8 KERNEL_NAME=MySims"

setup: projectroot deps kernel

$(SETTINGS):
	$(SETTINGS_JL) init

projectroot:
	$(JULIA) -e 'using Pkg; haskey(Pkg.project().dependencies, "ProjectRoot") ? @info("ProjectRoot already installed") : Pkg.add("ProjectRoot")'

deps:
	$(JULIA) --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# `--project=@.` makes the kernel use the Project.toml at or above the
# notebook's directory, so it works wherever the repository is cloned.
# Settings passed on the command line are saved only after the install succeeds.
kernel: $(SETTINGS)
	@$(SETTINGS_JL) resolve "$(THREADS)" > /dev/null
	$(JULIA) --project=. -e 'using IJulia; installkernel(ARGS[1], "--project=@.", "--threads=" * ARGS[2])' -- "$(KERNEL_NAME)" "$(NTHREADS)"
ifeq ($(origin THREADS),command line)
	$(SETTINGS_JL) set threads "$(THREADS)"
endif
ifeq ($(origin KERNEL_NAME),command line)
	$(SETTINGS_JL) set kernel_name "$(KERNEL_NAME)"
endif

kernel-remove:
	$(JULIA) --project=. -e 'using IJulia; d = joinpath(IJulia.kerneldir(), lowercase(replace(ARGS[1], " " => "-")) * "-$$(VERSION.major).$$(VERSION.minor)"); isdir(d) ? (rm(d; recursive=true); @info "Removed" d) : @warn "No such kernel" d' -- "$(KERNEL_NAME)"

# Interactive REPL; not meant to be combined with `make -j`.
session: $(SETTINGS)
	@$(SETTINGS_JL) resolve "$(THREADS)" > /dev/null
	$(JULIA) --project=$(CURDIR) --threads=$(NTHREADS)

test:
	$(JULIA) --project=. -e 'using Pkg; Pkg.test()'
