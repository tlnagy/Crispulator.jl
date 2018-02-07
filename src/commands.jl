# Sane behavior when run from the REPL
source_dir = typeof(Base.source_dir()) == Void ? joinpath(Pkg.dir("Crispulator"), "src") : Base.source_dir()
"""
List the custom simulation scripts that are available
"""
function ls()
    files = readdir(joinpath(source_dir, "exps"))
    filter(file -> file != "common.jl", files)
end

"""
bootstrap_exp(analysis_file::ASCIIString, output_dir::ASCIIString, proc_count::Int, debug::Bool)

Run custom simulation using one of the scripts in the `exps/` folder
"""
function bootstrap_exp(analysis_file::Compat.String,
                       output_dir::Compat.String,
                       proc_count::Int,
                       debug::Bool)

    # add additional threads
    addprocs(proc_count)
    println("Using $(nprocs()) threads")

    # load simulation code on all cores
    include(joinpath(source_dir, "simulation", "load.jl"))
    print("Loading analysis code...")
    analysis_full_path = joinpath(source_dir, "exps", analysis_file)

    (!isfile(analysis_full_path)) && error("Please provide a valid analysis file to run")
    (!isdir(dirname(abspath(output_dir)))) && error("Please provide a valid directory for output")

    # load common analysis code
    include(joinpath(source_dir, "exps", "common.jl"))

    # load specific analysis
    include(analysis_full_path)
    println("Done.")

    # compute output filename
    output_file = compute_name(output_dir)

    # Run
    println("Running $analysis_file and saving output in $output_file")
    func_name = Symbol(splitext(analysis_file)[1])
    getfield(Main, func_name)(output_file, debug=debug)
end

"""
bootstrap_config(config_file::String, output_dir::String, suppress_graph::Bool)

Run simulation using parameters supplied in a YAML configuration file
"""
function bootstrap_config(config_file::Compat.String,
                          output_dir::Compat.String,
                          suppress_graph::Bool)
    # load config
    config = YAML.load_file(config_file)

    if !isdir(output_dir)
        info("Directory $output_dir does not exist, attempting to create")
        try
            mkdir(output_dir)
        catch ex
            if isa(ex, SystemError)
                error("Unable to create $output_dir")
            end
        end
    end

    info("Using $(nprocs()) thread(s)")

    # load analysis files and simulation
    include(joinpath(source_dir, "simulation", "load.jl"))
    include(joinpath(source_dir, "utils", "parse_yaml.jl"))
    include(joinpath(source_dir, "utils", "runconfig.jl"))

    # run
    Base.invokelatest(runconfig, Base.invokelatest(parse, config)..., output_dir, suppress_graph)
end
