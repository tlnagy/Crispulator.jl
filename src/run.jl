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

@everywhere function run_exp(setup::ScreenSetup; run_idx=-1)

    lib = Library()
    guides, guide_freqs = construct_library(lib, setup.num_genes, setup.coverage)

    guide_count = setup.num_genes * setup.coverage
    cell_count = guide_count*setup.representation

    guide_freqs_dist = Categorical(guide_freqs)

    min_perc = minimum([range[2] - range[1] for (binname, range) in setup.bin_info])
    expand_to = round(Int64, setup.num_cells_per_bin/min_perc)

    cells = transfect(guides, guide_freqs_dist, cell_count, setup.moi, expand_to)

    bin_cells = facs_sort(cells, guides, setup.bin_info, setup.σ)

    freqs = counts_to_freqs(bin_cells, guide_count)
    raw_data = sequencing(Dict(:bin1=>setup.seq_depth,:bin2=>setup.seq_depth), guides, freqs)

    auroc = analyze(raw_data, gen_plots=false)

    [auroc; as_array(setup)...; run_idx]
end

function run_wrapper(filepath)
    representations = logspace(0, 3, 10)
    bin_sizes = 2.5*logspace(4,6,10)
    noises = [0.5, 1, 2]
    num_runs = 10

    runs = []
    for rep in representations, min_bin in bin_sizes, noise in noises
        for run in 1:num_runs
            setup = ScreenSetup()
            setup.representation = round(Int64, rep)
            setup.num_cells_per_bin = round(Int64, min_bin)
            setup.σ = noise
            push!(runs, (setup, run))
        end
    end

    results = @time pmap(args -> run_exp(args[1]; run_idx=args[2]), runs)
    results = collect(zip(results...))

    col_names = [:auroc; fieldnames(ScreenSetup)...; :run]

    results = DataFrame(Any[map(collect, results)...], col_names)
    delete!(results, :bin_info)

    writetable(filepath, results)
end

# fire up simulation if run using command line
if !isinteractive()
    run_wrapper("../data/output.csv")
end
