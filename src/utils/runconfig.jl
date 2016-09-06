include(normpath(joinpath(@__FILE__, "..", "..", "exps", "common.jl")))
using Iterators
using Gadfly

function runconfig(setup::ScreenSetup,
                   lib::Library,
                   num_runs::Int64,
                   output_dir::ASCIIString)


    info("Running config")
    methods = [venn, auprc]
    measures = [:inc, :dec, :incdec]
    genetypes = [:sigmoidal, :linear, :all]

    crisprtype = lib.cas9_behavior

    test_method_wrapper = (bc_counts, genes) -> test_methods(genes, methods, measures, genetypes)

    # a little abuse of lexical scoping
    gen_plots = (bc_counts, genes) -> begin
        info("Generating plots")
        new_filename = joinpath(output_dir, "counts.svg")
        nopseudo = bc_counts[(bc_counts[:counts_bin1] .> 0.5) & (bc_counts[:counts_bin2] .> 0.5), :]
        draw(SVG(new_filename, 10cm, 10cm), plot(nopseudo,
        x=:counts_bin1, y=:counts_bin2, color=:class, Scale.x_log10, Scale.y_log10,
        Theme(highlight_width=0pt), Coord.cartesian(fixed=true),
        Guide.title("$(typeof(crisprtype)) $(typeof(setup))"),
        Guide.xlabel("log counts t0"), Guide.ylabel("log counts t1")))
        new_filename = joinpath(output_dir, "volcano.svg")
        draw(SVG(new_filename, 10cm, 10cm), plot(genes,
        x=:mean, y=:pvalue, color=:class, Theme(highlight_width=0pt),
        Guide.title("$(typeof(crisprtype)) $(typeof(setup))"),
        Guide.xlabel("mean log2 fold change"), Guide.ylabel("-log10 pvalue")))
        test_methods(genes, methods, measures, genetypes)
    end

    funcs = [gen_plots, [test_method_wrapper for i in 1:(num_runs - 1)]...]

    # setup up runs
    runs = [(setup, lib, funcs[i], i) for i in 1:num_runs]

    results = pmap(args -> run_exp(args[1], args[2], args[3]; run_idx=args[4]), runs)

    info("Analyzing results")

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
    info("Saving results in $output_dir")
    writetable(joinpath(output_dir, "results_table.csv"), grouped_info)
end
