using Pkg
@info "Activating simulation environment"
Pkg.activate(@__DIR__)
@info "Instantiating environment"
Pkg.instantiate()

@info "Loading simulation framework"
using Crispulator

using ArgParse
using Distributed
using DataFrames
using Distributions
using Gadfly
using YAML

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
