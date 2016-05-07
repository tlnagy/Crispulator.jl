using DataFrames
using Distributions
using HypothesisTests
using StatsBase
using Iterators
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

@everywhere function run_exp(setup::ScreenSetup; run_idx=-1, testalg=false)

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

    auroc = analyze(raw_data, gen_plots=false, testmethod=testalg)

    [auroc...; as_array(setup)...; run_idx]
end

function build_parameter_space(parameters::Dict{Symbol, Vector}, num_runs::Int)
    fields = collect(keys(parameters))
    n_fields = length(fields)
    runs = []
    for vals in Iterators.product([parameters[field] for field in fields]...)
        for run in 1:num_runs
            setup = ScreenSetup()
            for idx in 1:n_fields
                setfield!(setup, fields[idx], vals[idx])
            end
            push!(runs, (setup, run))
        end
    end
    runs
end

function run_wrapper(filepath)
    parameters = Dict{Symbol, Vector}(
        :representation => map(x->round(Int64, x), logspace(0, 3, 2)),
        :num_cells_per_bin => map(x->round(Int64, x),  2.5*logspace(3,6,2)),
        :seq_depth => map(x->round(Int64, x),  logspace(5,7,2)),
        :σ => [0.5, 1, 2]
    )
    runs = build_parameter_space(parameters, 10)

    results = @time pmap(args -> run_exp(args[1]; run_idx=args[2], testalg=true), runs)
    results = DataFrame(hcat(results...)')

    names!(results, [:pval_auroc; :mean_auroc; :pvalmeanprod_auroc ; fieldnames(ScreenSetup)...; :run])
    # remove bin_info for now because there isn't a good way to encode
    # that information
    delete!(results, :bin_info)

    writetable(filepath, results)
end

# fire up simulation if run using command line
if !isinteractive()
    run_wrapper("../data/test_alg.csv")
end
