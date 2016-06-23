function build_parameter_space{T <: ScreenSetup}(::T, parameters::Dict{Symbol, Vector}, num_runs::Int)
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

function build_parameter_space{T <: ScreenSetup}(::T, libs::Vector{Library},
                                                 parameters::Dict{Symbol, Vector}, num_runs::Int)
    fields = collect(keys(parameters))
    n_fields = length(fields)
    runs = []
    for vals in Iterators.product([parameters[field] for field in fields]...)
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

function grouped_param_space{T <: ScreenSetup}(::T, parameters::Dict{Symbol, Vector}, dists::Vector{Symbol}, num_runs::Int)
    fields = collect(keys(parameters))
    n_fields = length(fields)
    deleteat!(fields, findin(fields, dists))
    runs = []
    grouped_params = zip([parameters[field] for field in fields]...)
    push!(fields, dists...)
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

@everywhere function test_methods(genes, methods, measures, genetypes)
    local results = []
    for (method, measure, genetype) in Iterators.product(methods, measures, genetypes)
        if genetype != :all
            subgene = genes[genes[:behavior] .== genetype, :]
        else
            subgene = genes
        end
        local result = 0.0
        if measure == :incdec
            result = method(abs(subgene[:pvalmeanprod]), subgene[:class], Set([:increasing, :decreasing]))
        elseif measure == :dec
            result = method(subgene[:pvalmeanprod], subgene[:class], Set([:decreasing]), rev=false)
        else
            result = method(subgene[:pvalmeanprod], subgene[:class], Set([:increasing]))
        end
        (typeof(result) <: Tuple) && (result = result[1])
        push!(results, result)
    end
    (results...)
end

function compute_name(filename::AbstractString)
    front, back = splitext(filename)
    # commit of the code that generated this data
    commit = strip(readall(`git rev-parse --short HEAD`))
    # whether there are any uncommitted changes
    n = length(split(readall(`git status --porcelain`), "\n"))
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
        row_labels = repeat(hierarchy[i, :], inner=[size(df, 1), 1])
        data = DataFrame(hcat(row_labels, df[i], Array(df[size(hierarchy, 1)+1:end])))
        names!(data, renames)
        push!(procs, data)
    end
    vcat(procs...)
end
