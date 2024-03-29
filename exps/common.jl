using IterTools
using Distributed

is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")

function build_parameter_space(::T, parameters::Dict{Symbol, Vector}, num_runs::Int) where {T <: ScreenSetup}
    fields = collect(keys(parameters))
    n_fields = length(fields)
    runs = []
    for vals in Iterators.product([parameters[field] for field in fields]...)
        for run in 1:num_runs
            setup = T()
            for idx in 1:n_fields
                setfield!(setup, fields[idx], vals[idx])
            end
            push!(runs, (setup, run))
        end
    end
    runs
end

function build_parameter_space(::T, parameters::Dict{Symbol, Vector},
                                    libs::Vector{Library}, num_runs::Int) where {T <: ScreenSetup}
    fields = collect(keys(parameters))
    n_fields = length(fields)
    runs = []
    for vals in IterTools.product([parameters[field] for field in fields]...)
        for lib in libs
            for run in 1:num_runs
                setup = T()
                for idx in 1:n_fields
                    setfield!(setup, fields[idx], vals[idx])
                end
                push!(runs, (setup, lib, run))
            end
        end
    end
    runs
end

function grouped_param_space(::T, parameters::Dict{Symbol, Vector},
                                  libs::Vector{Library},
                                  dists::Vector{Symbol}, num_runs::Int) where {T <: ScreenSetup}
    fields = collect(keys(parameters))
    n_fields = length(fields)
    deleteat!(fields, findall(in(dists), fields))
    runs = []
    grouped_params = zip([parameters[field] for field in fields]...)
    (length(dists) > 0) && push!(fields, dists...)
    for vals in Iterators.product(grouped_params, [parameters[dist] for dist in dists]...)
        vals = Any[vals[1]...; [vals[i+1] for i in 1:length(dists)]...]
        for lib in libs
            for run in 1:num_runs
                setup = T()
                for idx in 1:n_fields
                    setfield!(setup, fields[idx], vals[idx])
                end
                push!(runs, (setup, lib, run))
            end
        end
    end
    runs
end

function grouped_param_space(::T, parameters::Dict{Symbol, Vector}, 
                                  dists::Vector{Symbol}, num_runs::Int) where {T <: ScreenSetup}
    fields = collect(keys(parameters))
    n_fields = length(fields)
    deleteat!(fields, findall(in(dists), fields))
    runs = []
    grouped_params = zip([parameters[field] for field in fields]...)
    (length(dists) > 0) && push!(fields, dists...)
    for vals in Iterators.product(grouped_params, [parameters[dist] for dist in dists]...)
        vals = Any[vals[1]...; [vals[i+1] for i in 1:length(dists)]...]
        for run in 1:num_runs
            setup = T()
            for idx in 1:n_fields
                setfield!(setup, fields[idx], vals[idx])
            end
            push!(runs, (setup, run))
        end
    end
    runs
end

function test_methods(genes, methods, measures, genetypes)
    local results = []
    for (method, measure, genetype) in Iterators.product(methods, measures, genetypes)
        if genetype != :all
            subgene = genes[genes[!, :behavior] .== genetype, :]
        else
            subgene = genes
        end
        local result = 0.0
        if measure == :incdec
            result = method(abs.(subgene[!, :pvalmeanprod_bin2_div_bin1]), subgene[!, :class], Set([:increasing, :decreasing]))
        elseif measure == :dec
            result = method(subgene[!, :pvalmeanprod_bin2_div_bin1], subgene[!, :class], Set([:decreasing]), rev=false)
        else
            result = method(subgene[!, :pvalmeanprod_bin2_div_bin1], subgene[!, :class], Set([:increasing]))
        end
        (typeof(result) <: Tuple) && (result = result[1])
        push!(results, result)
    end
    (results..., )
end

function test_methods_snr(bc_counts, genes, methods, measures, genetypes)
    (test_methods(genes, methods, measures, genetypes)..., Crispulator.signal(bc_counts), Crispulator.noise(bc_counts))
end

function compute_snr(bc_counts::DataFrame, genes::DataFrame)
    sig, noi = Crispulator.signal(bc_counts), Crispulator.noise(bc_counts)
    (sig/noi, sig, noi)
end

function compute_name(filename::AbstractString)
    front, back = splitext(filename)
    # commit of the code that generated this data
    commit = strip(read(`git rev-parse --short HEAD`, String))
    # whether there are any uncommitted changes
    n = length(split(read(`git status --porcelain`, String), "\n"))
    status = n > 1 ? "dirty" : "clean"
    string(join([front; commit; status], "_"), back)
end

"""
    construct_hierarchical_label(hierarchy::Array{Symbol, 2}, df::DataFrame, renames::Vector{Symbol})

Creates a hierarchical index in the given dataframe
"""
function construct_hierarchical_label(hierarchy::Array{Symbol, 2}, df::DataFrame, renames::Vector{Symbol})
    procs = []
    for i in 1:size(hierarchy, 1)
        row_labels = repeat(hierarchy[i:i, :], inner=[size(df, 1), 1])
        data = DataFrame(hcat(row_labels, df[!, i], Array(df[!, size(hierarchy, 1)+1:end])), :auto)
        rename!(data, renames)
        push!(procs, data)
    end
    vcat(procs..., )
end
