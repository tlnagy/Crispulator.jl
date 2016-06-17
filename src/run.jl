(!isdir(Pkg.dir("ArgParse"))) && Pkg.add("ArgParse")
using ArgParse

function parse_cmdline()
    settings = ArgParseSettings(
        description="\033[32mPooled Screen Optimizer\033[0m\n\n\n\n" *
        "\033[31mHigh performance simulation and exploratory analysis code "*
        "for studying the behavior of pooled screens under common " *
        "experimental conditions.\033[0m"
    )

    @add_arg_table settings begin
        "--debug"
            help = "Use smaller parameter space for diagnosing issues"
            action = :store_true
        "analysis_file"
            help = "File from src/analyses/ to run"
            required = true
        "output_file"
            help = "File to output results to [.CSV, .TSV, etc]"
            required = true
    end

    settings.epilog = readall(normpath(joinpath(Base.source_dir(),"..","LICENSE")))

    return parse_args(settings)
end


function main()
    parsed_args = parse_cmdline()
    # target_file = parsed_args[:analysis_file]
    # filename = ARGS[2]
    println("Using $(nprocs()) threads")
    include(joinpath(Base.source_dir(), "simulation", "load.jl"))
    print("Loading analysis code...")
    target_file = parsed_args["analysis_file"]
    target_full = joinpath(Base.source_dir(), "exps", target_file)
    output_file = parsed_args["output_file"]
    (!isfile(target_full)) && error("Please provide a valid analysis file to run")
    (!isdir(dirname(abspath(output_file)))) && error("Please provide a valid directory for output")
    include(joinpath(Base.source_dir(), "exps", "common.jl"))
    include(target_full)
    println("Done.")
    output_file = compute_name(output_file)
    println("Running $target_file and saving output in $output_file")
    main(output_file, debug=parsed_args["debug"])
end

# fire up simulation if run using command line
if !isinteractive()
    main()
end
