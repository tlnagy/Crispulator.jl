include(normpath(joinpath(@__FILE__, "..", "..", "exps", "common.jl")))
using Iterators

function runconfig(setup::ScreenSetup, lib::Library, num_runs::Int64)

    # setup up runs
    runs = [(setup, lib, i) for i in 1:num_runs]

    methods = [venn, auprc]
    measures = [:inc, :dec, :incdec]
    genetypes = [:sigmoidal, :linear, :all]
    test_method_wrapper = (bc_counts, genes) -> test_methods(genes, methods, measures, genetypes)

    results = pmap(args -> run_exp(args[1], args[2],
                   test_method_wrapper;
                   run_idx=args[3]), runs)
    results = DataFrame(permutedims(hcat(results...), [2, 1]))
    hierarchy = vcat([hcat(item...) for item in Iterators.product(map(Symbol, methods), measures, genetypes)]...)
    new_names = [[:method, :measure, :genetype, :score]...; fieldnames(setup)...; :run_idx]
    results = construct_hierarchical_label(hierarchy, results, new_names)

    grouped_info = by(results, [:method, :measure, :genetype]) do grouped_df
        n = size(grouped_df, 1)
        mean_score = mean(grouped_df[:score])
        std_score = std(grouped_df[:score])
        conf_int = 2.58 * std_score./sqrt(n)
        DataFrame(
            std_score = std_score,
            mean_score = mean_score,
            conf_max = mean_score + conf_int,
            conf_min = mean_score - conf_int,
            n = n
        )
    end
    println(grouped_info)
end
