# CRISPulator

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://tlnagy.github.io/Crispulator.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://tlnagy.github.io/Crispulator.jl/latest)
[![Build
Status](https://travis-ci.com/tlnagy/Crispulator.jl.svg?token=MCUYuFeh1dFnAvCDpb4q&branch=master)](https://travis-ci.com/tlnagy/Crispulator.jl)
[![Build
status](https://ci.appveyor.com/api/projects/status/nr33nmqu6hvowni5/branch/master?svg=true)](https://ci.appveyor.com/project/tlnagy/crispulator-jl/branch/master)
[![License](http://img.shields.io/:license-apache-blue.svg?style=flat-square)](http://www.apache.org/licenses/LICENSE-2.0.html)

A pooled genetic screen simulator

## Motivation

Pooled screens are very useful, but they are also extraordinarily complex.
This work aims to explore the importance of various screen parameters and
the relevance of design choices on the downstream analysis and
interpretation.

## Setup

To run the simulation, you will need a recent version of
[Julia](http://julialang.org/downloads/) (v0.5+) installed. Start
Julia and enter the following command in the Julia command line interface (REPL):

```julia
Pkg.update(); Pkg.clone("https://github.com/tlnagy/Crispulator.jl.git"); Pkg.build("Crispulator")
```

this downloads `Crispulator` and installs all of its dependencies. We also
recommend that you run:

```julia
Pkg.test("Crispulator")
```

which loads and runs all of the project's internal tests and makes sure
that everything is ready.

## Usage via Julia REPL/Jupyter Notebooks

Once installed, run the following commands into the Julia command line:

```julia
include(joinpath(Pkg.dir("Crispulator"), "src", "run.jl"))
```

this loads the runscript into the REPL. Then run

```julia
cp(joinpath(Pkg.dir("Crispulator"), "example_config.yml"), joinpath(pwd(), "custom_config.yml"))
```

this copies the YAML configuration file over into your current working
directory. Edit this file to change the parameters of the simulation. Any
standard text editor should be able to open it (e.g. Vim, Atom,
SublimeText, TextEdit, Notepad, etc). If you are using a Jupyter Notebook
then it should appear under the Files tab, clicking on it will let you
edit it.

Once you've made your changes, run:

```julia
bootstrap_config(joinpath(pwd(), "custom_config.yml"), "test_output/", false)
```

which runs the simulation according to the parameters in
`custom_config.yml` and saves the output in the `test_output` directory.
You can go back and edit the `custom_config.yml` and rerun the code, but
we recommend that you run

```julia
workspace()
```

to clear the workspace so no overwrite warnings get displayed.

## Usage via terminal

While REPL usage is convenient, it can be beneficial to interact with
`Crispulator` via the terminal. First, make sure that julia is in your
PATH (i.e. check that `julia -e "versioninfo()"` runs properly and shows
0.5+ as the version number). Then from the `Crispulator` package root
directory, run:

```
julia src/run.jl config example_config.yml test_output
```

where `config` tells `CRISPulator` to use the provided config
`example_config.yml` and `test_output` is the directory where the results
will be saved. This directory will be created if it doesn't exist.
