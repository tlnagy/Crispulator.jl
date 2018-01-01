using ArgParse
include(joinpath(Pkg.dir("Crispulator"), "src", "simulation", "load.jl"))

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--screen_type"
            help = "Type of screening experiment"
			arg_type = String
			default = "FACS"
		"--num_replicates"
			help = "The number of replicates for each group"
			arg_type = Int
			default = 4
		"--num_genes"
			help = "The number of genes on genome"
			arg_type = Int
			default = 1000
		"--num_guides"
			help = "The number of guides for each target gene"
			arg_type = Int
			default = 10
		"--bin_ratio"
			help = "FACS ratio"
			arg_type = Float64
			default = 0.25
		"--positive_ratio"
			help = "Ratio of positive phenotype genes"
			arg_type = Float64
			default = 0.2
		"--noise_ratio"
			help = "Noise ratio"
			arg_type = Float64
			default = 0.1
		"--seq_depth"
			help = "Sequencing Depth"
			arg_type = Int
			default = 1000
		"--representation"
			help = "Bottleneck representation"
			arg_type = Int
			default = 100
        "--num_bottlenecks"
            help = "Bottleneck level for growth-screening"
            arg_type = Int
            default = 10
        "--random_seed"
            help = "Random seed"
            arg_type = Int
            default = 0
        "--output_path"
            help = "Path of output directory"
            arg_type = String
            default = "output"
    end
    return parse_args(s)
end

function generate(args)
    srand(args["random_seed"])

    if args["screen_type"] == "FACS"
        facs_param = FacsScreen()
        facs_param.num_genes = args["num_genes"]
        facs_param.coverage = args["num_guides"]
        facs_param.representation = args["representation"]
        bin_prob = args["bin_ratio"]
        facs_param.bin_info[:bin1] = (0, bin_prob)
        facs_param.bin_info[:bin2] = (1-bin_prob, 1)
        facs_param.seq_depth = args["seq_depth"]
        facs_param.σ = args["noise_ratio"]
        pheno_prob = args["positive_ratio"]

        lib = Library(CRISPRn())
        lib.phenotype_probs = Categorical([pheno_prob/2,0.05,0.95-pheno_prob,pheno_prob/2])
        guides, guide_freqs_dist = construct_library(facs_param, lib)
    else
        growth_param = GrowthScreen()
        growth_param.num_genes = args["num_genes"]
        growth_param.coverage = args["num_guides"]
        growth_param.representation = args["representation"]
        growth_param.seq_depth = args["seq_depth"]
        growth_param.noise = args["noise_ratio"]
        growth_param.num_bottlenecks = args["num_bottlenecks"]
        pheno_prob = args["positive_ratio"]
        max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
            :inactive => (1-0.05-pheno_prob, Delta(0.0)),
            :negcontrol => (0.05, Delta(0.0)),
            :increasing => (pheno_prob/2, TruncatedNormal(0.1, 0.1, 0.025, 1)),
            :decreasing => (pheno_prob/2, TruncatedNormal(-0.55, 0.2, -1, -0.1))
        )

        lib = Library(max_phenotype_dists, CRISPRn())
        guides, guide_freqs_dist = construct_library(growth_param, lib)
    end

    dir_path = string(args["output_path"], "/")
    println(dir_path)
    mkpath(dir_path)
    num_sample = args["num_replicates"]
    for rep_no in 1:num_sample
        println("Running $(rep_no)...")

        if args["screen_type"] == "FACS"
            cells, cell_phenotypes = transfect(facs_param, lib, guides, guide_freqs_dist)
            bin_cells = select(facs_param, cells, cell_phenotypes, guides)
            freqs = counts_to_freqs(bin_cells, length(guides))

            seq_depths = Dict{Symbol, Int}()
            for binname in keys(bin_cells)
                seq_depths[binname] = rand(Poisson(facs_param.seq_depth))
            end
        else
            cells, cell_phenotypes = transfect(growth_param, lib, guides, guide_freqs_dist)
            bin_cells = select(growth_param, cells, cell_phenotypes, guides)
            freqs = counts_to_freqs(bin_cells, length(guides))
            seq_depths = Dict{Symbol, Int}()
            for binname in keys(bin_cells)
                seq_depths[binname] = rand(Poisson(growth_param.seq_depth))
            end
        end

        raw_data = sequencing(seq_depths, guides, freqs)
        bc_counts, genes = differences_between_bins(raw_data)

        file_name = string(dir_path,"/bc_count_$(rep_no).csv")
        writetable(file_name, bc_counts)
    end
end

function main()
    parsed_args = parse_commandline()
    generate(parsed_args)
end

main()
