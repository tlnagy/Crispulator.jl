# Pooled Screen Optimization

[![Build Status](https://travis-ci.com/tlnagy/pooled-screen-optimization.svg?token=MCUYuFeh1dFnAvCDpb4q&branch=master)](https://travis-ci.com/tlnagy/pooled-screen-optimization)
[![License](http://img.shields.io/:license-apache-blue.svg?style=flat-square)](http://www.apache.org/licenses/LICENSE-2.0.html)

## Motivation

Pooled screens are very useful, but they are also extraordinarily complex.
This work aims to explore the importance of various screen parameters and
the relevance of design choices on the downstream analysis and
interpretation. 

## Usage

To run the simulation, you will need a recent version of
[Julia](http://julialang.org) installed and in your PATH. Then navigate
into the root directory of the project and run

```
julia -p N src/run.jl {{experiment_file.jl}} {{output_filepath.csv}}
```

where `N+1` is the number of total workers to use for simulation,
`{{experiment_file.jl` is the filename of the simulation to load from
`src/exps/` and `{{output_filepath.csv}}` is the directory and filename
that will be used for storing the results of the simulation.

## Testing

To run the unit tests do

```
julia test/runtests.jl
```
