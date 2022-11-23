if Base.current_project() != joinpath(@__DIR__, "Project.toml") || (get(ENV, "CI", nothing) == "true")
    @info "Activating simulation environment"
    using Pkg
    Pkg.activate(@__DIR__)
    @info "Instantiating environment"
    Pkg.instantiate()
end
using Distributed
using ProgressMeter

@everywhere begin
    using Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
end

@info "Loading simulation framework"
@everywhere begin
    using Crispulator

    using ArgParse
    using DataFrames
    using Distributions
    using DataStructures
    using CSV
    using Gadfly
    using YAML
end

include("parsing.jl")
include("commands.jl")

function main()
    parsed_args = parse_args(build_arg_table())

    command = parsed_args["%COMMAND%"]

    if command == "ls"
        foreach(x -> println(x), ls())
    elseif command == "exp"
        bootstrap_exp(
            parsed_args[command]["analysis_file"],
            parsed_args[command]["output_file"],
            parsed_args[command]["addprocs"],
            parsed_args[command]["debug"]
        )
    else
        bootstrap_config(
            parsed_args[command]["config_file"],
            parsed_args[command]["output_dir"],
            parsed_args[command]["no-graph"]
        )
    end
end

# fire up simulation if run using command line
if !isinteractive()
    main()
end
