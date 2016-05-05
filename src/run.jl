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
                               num_cells_per_bin = 2e6,
                               seq_depth = 10^7
                            )

    lib = Library()
    guides, guide_freqs = construct_library(lib, num_genes, coverage)

    guide_count = num_genes * coverage
    cell_count = guide_count*representation

    guide_freqs_dist = Categorical(guide_freqs)

    min_perc = minimum([range[2] - range[1] for (binname, range) in bin_info])
    expand_to = round(Int64, num_cells_per_bin/min_perc)

    cells = transfect(guides, guide_freqs_dist, cell_count, moi, expand_to)

    bin_cells = facs_sort(cells, guides, bin_info, σ)

    freqs = counts_to_freqs(bin_cells, guide_count)
    raw_data = sequencing(Dict(:bin1=>seq_depth,:bin2=>seq_depth), guides, freqs)

    auroc = analyze(raw_data, gen_plots=false)

    auroc, representation
end

function run_wrapper()
    representations = [1, 5, 10, 50, 100, 500, 1000]
    num_runs = 10

    runs = vec([(rep,run) for run in 1:num_runs, rep in representations])

    results = @time pmap(args -> run_exp(; representation = args[1]), runs)
    results = collect(zip(results...))
    writetable("../data/output.csv", DataFrame(auroc=[results[1]...], representation=[results[2]...]))
end

# fire up simulation if run using command line
if !isinteractive()
    run_wrapper()
end
