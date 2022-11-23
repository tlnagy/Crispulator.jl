# CRISPulator

*A pooled genetic screen simulator*

Pooled screens are very useful, but it is difficult to select all of the
necessary parameters *a priori*. This work aims to explore the importance
of various screen parameters and the relevance of design choices on the
downstream analysis and interpretation.

!!! compat "Paper"
    Nagy, T. & Kampmann, M. CRISPulator: a discrete simulation tool for pooled
    genetic screens. BMC Bioinformatics 18, 347 (2017).
    <https://doi.org/10.1186/s12859-017-1759-9> (Preprint also
    [available](https://doi.org/10.1101/119131)) 

## Setup

To run the simulation, you will need a recent version of
[Julia](http://julialang.org/downloads/) installed and in your PATH. Then
navigate into the root directory of the project and run `julia`. Run the
following command:

```sh
$ julia -e 'using Pkg; Pkg.develop(url="https://github.com/tlnagy/Crispulator.jl.git")'
```

this downloads `Crispulator` and installs it into the Julia package directory
with all of its dependencies. Note the `$` just represents the prompt, you don't
need to type it. Make sure to run the tests and verify that everything is
passing 

```sh
$ julia -e 'Pkg.test("Crispulator")'
```

## Quickstart

For most simple cases, no writing of Julia code is necessary. See the
[Tutorial](@ref) for using Crispulator in this manner.


## Advanced

See the [Custom Simulations](@ref) for a step-by-step guide to writing a custom
simulation on top of Crispulator. Additionally, the [Crispulator Internals](@ref)
page has in-depth documentation needed for more advanced usage
