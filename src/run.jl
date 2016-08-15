using ArgParse
using YAML
using Compat
import Compat: UTF8String, ASCIIString, view, readstring

function parse_cmdline()
    settings = ArgParseSettings(
        description="\033[32mPooled Screen Optimizer\033[0m\n\n\n\n" *
        "\033[31mHigh performance simulation and exploratory analysis code "*
        "for studying the behavior of pooled screens under common " *
        "experimental conditions.\033[0m"
    )

    @add_arg_table settings begin
        "config"
            action = :command
            help = "run simulation using a YAML config file (easy)"
        "exp"
            action = :command
            help = "run modified simulation from src/exps/ (hard)"
        "ls"
            action = :command
            help = "list all modified simulations available"
    end

    settings["ls"].description = "Prints out the available modified simulation " *
    "code. Simply prints out files in src/exps/."

    @add_arg_table settings["exp"] begin
        "--debug"
            help = "Use smaller parameter space for diagnosing issues"
            action = :store_true
        "analysis_file"
            help = "Run a modified simulation from the src/exps directory"
            required = true
        "output_file"
            help = "File to output results to [.CSV, .TSV, etc]"
            required = true
    end

    settings["exp"].description = "Run a modified simulation. Much more" *
    " flexible than using the configuration files, but also more complex"

    @add_arg_table settings["config"] begin
        "config_file"
            help = "path to YAML configuration file to load"
            required = true
        "output_dir"
            help = "path to output directory"
            required = true
    end

    settings["config"].description

    settings.epilog = readstring(normpath(joinpath(Base.source_dir(),"..","LICENSE")))

    return parse_args(settings)
end


function main()
    parsed_args = parse_cmdline()

    command = parsed_args["%COMMAND%"]

    if command == "ls"
        files = readdir(joinpath(Base.source_dir(), "exps"))
        for file in files
            (file == "common.jl") && continue
            println(file)
        end
    elseif command == "exp"
        println("Using $(nprocs()) threads")
        include(joinpath(Base.source_dir(), "simulation", "load.jl"))
        print("Loading analysis code...")
        target_file = parsed_args[command]["analysis_file"]
        target_full = joinpath(Base.source_dir(), "exps", target_file)
        output_file = parsed_args[command]["output_file"]
        (!isfile(target_full)) && error("Please provide a valid analysis file to run")
        (!isdir(dirname(abspath(output_file)))) && error("Please provide a valid directory for output")
        include(joinpath(Base.source_dir(), "exps", "common.jl"))
        include(target_full)
        println("Done.")
        output_file = compute_name(output_file)
        println("Running $target_file and saving output in $output_file")
        main(output_file, debug=parsed_args[command]["debug"])
    else
        config = YAML.load_file(parsed_args[command]["config_file"])
        (!isdir(parsed_args[command]["output_dir"])) && error("Please provide a valid directory for output")
        println("Using $(nprocs()) threads")
        include(joinpath(Base.source_dir(), "simulation", "load.jl"))
        include(joinpath(Base.source_dir(), "utils", "parse_yaml.jl"))
        include(joinpath(Base.source_dir(), "utils", "runconfig.jl"))
        runconfig(parse(config)...)
    end
end

# fire up simulation if run using command line
if !isinteractive()
    main()
end
