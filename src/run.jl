using DataFrames
using Distributions
using HypothesisTests
using StatsBase
@everywhere using DataFrames
@everywhere using Distributions
@everywhere using HypothesisTests
@everywhere using StatsBase

@everywhere function include_all()
    filenames = ["common.jl", "utils.jl", "library.jl", "transfection.jl",
                 "selection.jl", "sequencing.jl", "analysis.jl"]
    for filename in filenames
        include(filename)
    end
end
@everywhere include_all()

@everywhere function run_exp(; num_genes = 500, # number of target genes
                               coverage = 5, # number of guides per gene
                               representation = 100, # Number of cells with each guide
                               moi = 0.25, # multiplicity of infection
                               σ = 1.0, # std dev expected for cells during facs sorting
                               bin_info = Dict(:bin1 => (0.0, 1/3), :bin2 => (2/3, 1.0)),
                               seq_depth = 10^7
                            )

    lib = Library()
    guides, guide_freqs = construct_library(lib, num_genes, coverage)

    guide_count = num_genes * coverage
    cell_count = guide_count*representation

    guide_freqs_dist = Categorical(guide_freqs)

    cells = transfect(guides, guide_freqs_dist, cell_count, moi)

    bin_cells = facs_sort(cells, guides, bin_info, σ)

    freqs = counts_to_freqs(bin_cells, guide_count)
    raw_data = sequencing(Dict(:bin1=>seq_depth,:bin2=>seq_depth), guides, freqs)

    auroc = analyze(raw_data, gen_plots=false)

    auroc, representation
end

function run_wrapper()
    representations = [1000]
    # representations = [1, 5, 10, 50, 100, 500, 1000]
    num_runs = 10

    runs = vec([(rep,run) for rep in representations, run in 1:num_runs])

    results = @time pmap(args -> run_exp(; representation = args[1]), runs)

    for val in results; println(val) end
end

# fire up simulation if run using command line
if !isinteractive()
    run_wrapper()
end
