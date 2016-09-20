# CRISPulator

[![Build Status](https://travis-ci.com/tlnagy/Crispulator.jl.svg?token=MCUYuFeh1dFnAvCDpb4q&branch=master)](https://travis-ci.com/tlnagy/Crispulator.jl)
[![License](http://img.shields.io/:license-apache-blue.svg?style=flat-square)](http://www.apache.org/licenses/LICENSE-2.0.html)

A pooled genetic screen simulator

## Motivation

Pooled screens are very useful, but they are also extraordinarily complex.
This work aims to explore the importance of various screen parameters and
the relevance of design choices on the downstream analysis and
interpretation.

## Setup

To run the simulation, you will need a recent version of
[Julia](http://julialang.org/downloads/) installed and in your PATH. Start
Julia and enter the following command in the REPL:

```
Pkg.clone("git://github.com/tlnagy/Crispulator.jl.git"); Pkg.build("Crispulator")
```

this downloads `Crispulator` and installs all of its dependencies. We also
recommend that you run:

```
Pkg.test("Crispulator")
```

which loads and runs all of the project's internal tests and makes sure
that everything is ready.

## Usage

Once installed, run the following from the REPL:

```
cd(Pkg.dir("Crispulator"))
run(`julia src/run.jl config example_config.yml test_output`)
```

where `config` tells CRISPulator to use the provided config `example_config.yml`
and `test_output` is the directory where the results will be saved. This
directory will be created if it doesn't exist.
