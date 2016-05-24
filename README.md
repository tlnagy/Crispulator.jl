# Pooled Screen Optimization

## Motivation

Pooled screens are very useful, but they are also extraordinarily complex.
This work aims to explore the importance of various screen parameters and
the relevance of design choices on the downstream analysis and
interpretation. 

## Usage

To run the simulation, you will need a recent version of
[Julia](http://julialang.org) installed and in your PATH. Then navigate into
the `src/` directory and run

```
julia -p N run.jl scan_best_methods ../data/output.csv
```

## Testing

To run the unit tests do

```
julia test/runtests.jl
```
