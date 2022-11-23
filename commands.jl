"""
List the custom simulation scripts that are available
"""
function ls()
    files = readdir(joinpath(@__DIR__, "exps"))
    filter(file -> file != "common.jl", files)
end

"""
    bootstrap_exp(analysis_file::String, output_dir::String, proc_count::Int, debug::Bool)

Run custom simulation using one of the scripts in the `exps/` folder
"""
function bootstrap_exp(analysis_file::String,
                       output_dir::String,
                       proc_count::Int,
                       debug::Bool)

    # add additional threads
    addprocs(proc_count)
    @info "Using $(nprocs()) threads"

    # load simulation code on all cores
    @info "Loading analysis code"
    analysis_full_path = joinpath(@__DIR__, "exps", analysis_file)

    (!isfile(analysis_full_path)) && error("Please provide a valid analysis file to run")
    (!isdir(dirname(abspath(output_dir)))) && error("Please provide a valid directory for output")

    # load common analysis code
    include(joinpath(@__DIR__, "exps", "common.jl"))

    # load specific analysis
    include(analysis_full_path)

    # compute output filename
    output_file = Base.invokelatest(compute_name, output_dir)

    # Run
    @info "Running $analysis_file and saving output in $output_file"
    func_name = Symbol(splitext(analysis_file)[1])
    # getfield(Main, func_name)(output_file; debug=debug, quiet = false)
    Base.invokelatest(getfield(Main, func_name), output_file; debug = debug)
    @info "Done"
end

"""
bootstrap_config(config_file::String, output_dir::String, suppress_graph::Bool)

Run simulation using parameters supplied in a YAML configuration file
"""
function bootstrap_config(config_file::String,
                          output_dir::String,
                          suppress_graph::Bool)
    # load config
    config = YAML.load_file(config_file)

    if !isdir(output_dir)
        @info "Directory $output_dir does not exist, attempting to create"
        try
            mkdir(output_dir)
        catch ex
            if isa(ex, SystemError)
                error("Unable to create $output_dir")
            end
        end
    end

    @info "Using $(nprocs()) thread(s)"

    # load analysis files and simulation
    include(joinpath(@__DIR__, "utils", "parse_yaml.jl"))
    include(joinpath(@__DIR__, "utils", "runconfig.jl"))

    # run
    Base.invokelatest(runconfig, Base.invokelatest(parse, config)..., output_dir, suppress_graph)
end
