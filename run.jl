# Sane behavior when run from the REPL
using Pkg
source_dir = typeof(Base.source_dir()) === nothing ? joinpath(Pkg.dir("Crispulator")) : Base.source_dir()
Pkg.activate(source_dir)

using ArgParse
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
