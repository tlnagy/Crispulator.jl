# Multiprocessing

If you want to scan a wider range of conditions then `Crispulator.jl`'s
multiprocessing capabilities can come in handy. Say you wanted to re-run one of
the experiments from the paper, e.g. 
[Fig
11](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1759-9/figures/11), 
we need to scan a decent number of conditions to generate the data. Accelerating
this process by allowing the simulation to run on many cores at once can help
make it more feasible. We can do this by supplying Julia with the `-p` flag
followed by an integer of the number of additional cores we want to use.

```julia
julia -p 1 run.jl exp growth_sensitivity_library.jl growth_sensitivity_lib.csv
```

```@example
io = IOBuffer(); #hide
run(pipeline(`julia -p 1 ../../run.jl exp growth_sensitivity_library.jl growth_sensitivity_lib.csv --debug`; stdout = io)) #hide
s = String(take!(io))
println(join(filter(x->occursin("[ Info:",x), split(s, "\n")), "\n")) #hide
```

So in this case we'll be using 1+1 cores to run the simulation. Naturally, on a
big workstation/HPC you can scale up much more.

!!! warning
    `Crispulator.jl` can start to gobble up a lot of RAM with many cores. I
    would try increasing by one or two to see how the memory usage scales prior
    to trying to use all cores with the `-p auto` flag.

The output will be saved to the provided filename (with a suffix describing the
state of the git repository).