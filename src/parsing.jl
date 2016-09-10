function build_arg_table()
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
        "--addprocs", "-p"
            help = "Add additional processors"
            arg_type = Int
            default = 0
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
        "--no-graph", "-g"
            help = "Suppress graphical output"
            action = :store_true
        "config_file"
            help = "path to YAML configuration file to load"
            required = true
        "output_dir"
            help = "path to output directory"
            required = true
    end

    settings["config"].description

    if typeof(Base.source_dir()) != Void
        settings.epilog = readstring(normpath(joinpath(Base.source_dir(),"..","LICENSE")))
    end

    return settings
end
