# CRISPulator

*A pooled genetic screen simulator*

Pooled screens are very useful, but it is difficult to select all of the
necessary parameters *a priori*. This work aims to explore the importance
of various screen parameters and the relevance of design choices on the
downstream analysis and interpretation.

## Setup

To run the simulation, you will need a recent version of
[Julia](http://julialang.org/downloads/) installed and in your PATH. Then
navigate into the root directory of the project and run `julia`. Run the
following command:

```
julia -e 'Pkg.update(); Pkg.add("Crispulator")'
```

this copies `Crispulator` over to the Julia package directory and installs
all of its dependencies.

## Quickstart

From the root directory of the project run

```julia
julia src/run.jl config example_config.yml .
```

## Advanced

See the [Simulation Internals](@ref) page for in-depth documentation needed for more advanced usage
