include(normpath(joinpath(@__FILE__, "..", "..", "exps", "common.jl")))
using IterTools
using CSV

function runconfig(setup::ScreenSetup,
                   lib::Library,
                   num_runs::Int64,
                   output_dir::Compat.String,
                   suppress_graph::Bool)


    info("Running config")
    methods = [venn, auprc]
    measures = [:inc, :dec, :incdec]
    genetypes = [:sigmoidal, :linear, :all]

    crisprtype = lib.cas9_behavior

    test_method_wrapper = (bc_counts, genes) -> test_methods_snr(bc_counts, genes, methods, measures, genetypes)

    # a little abuse of lexical scoping
    gen_plots = (bc_counts, genes) -> begin
        info("Generating plots")
        new_filename = joinpath(output_dir, "counts.svg")
        nopseudo = bc_counts[(bc_counts[:counts_bin1] .> 0.5) .& (bc_counts[:counts_bin2] .> 0.5), :]
        draw(SVG(new_filename, 10cm, 10cm), plot(nopseudo,
        x=:counts_bin1, y=:counts_bin2, color=:class, Scale.x_log10, Scale.y_log10,
        Theme(highlight_width=0pt), Coord.cartesian(fixed=true),
        Guide.title("$(typeof(crisprtype)) $(typeof(setup))"),
        Guide.xlabel("log counts t0"), Guide.ylabel("log counts t1")))
        new_filename = joinpath(output_dir, "volcano.svg")
        draw(SVG(new_filename, 10cm, 10cm), plot(genes,
        x=:mean_bin2_div_bin1, y=:pvalue_bin2_div_bin1, color=:class, Theme(highlight_width=0pt),
        Guide.title("$(typeof(crisprtype)) $(typeof(setup))"),
        Guide.xlabel("mean log2 fold change"), Guide.ylabel("-log10 pvalue")))
        test_methods_snr(bc_counts, genes, methods, measures, genetypes)
    end

    if suppress_graph
        funcs = [test_method_wrapper for i in 1:num_runs]
    else
        eval(:(using Gadfly))
        funcs = [gen_plots, [test_method_wrapper for i in 1:(num_runs - 1)]...]
    end

    # setup up runs
    runs = [(setup, lib, funcs[i], i) for i in 1:num_runs]

    results = pmap(args -> run_exp(args[1], args[2], args[3]; run_idx=args[4]), runs)

    info("Analyzing results")

    results = DataFrame(permutedims(hcat(results...), [2, 1]))
    mean_snr = mean(results[:x19]./results[:x20])
    std_snr = std(results[:x19]./results[:x20])
    snr_result = "SNR score = $(round(mean_snr, 3)) +/- $(round(std_snr, 3))"
    delete!(results, [:x19, :x20])
    hierarchy = vcat([hcat(item...) for item in IterTools.product(map(Symbol, methods), measures, genetypes)]...)
    new_names = [[:method, :measure, :genetype, :score]...; fieldnames(setup)...; :run_idx]
    results = construct_hierarchical_label(hierarchy, results, new_names)

    venn_result, auprc_result = "", ""

    grouped_info = by(results, [:method, :measure, :genetype]) do grouped_df
        n = size(grouped_df, 1)
        mean_score = mean(grouped_df[:score])
        std_score = std(float.(grouped_df[:score]))
        conf_int = 2.58 * std_score./sqrt(n)

        if (grouped_df[1, :measure] == :incdec) &&
           (grouped_df[1, :genetype] == :all)

            disp_score = round(mean_score, 3)
            conf_max = round(mean_score + conf_int, 3)
            conf_min = round(mean_score - conf_int, 3)

            if (grouped_df[1, :method] == :venn)
                venn_result = "Venn score = $disp_score, 95% conf int ($conf_min, $conf_max)"
            else
                auprc_result = "AUPRC score = $disp_score, 95% conf int ($conf_min, $conf_max)"
            end
        end

        DataFrame(
            std_score = std_score,
            mean_score = mean_score,
            conf_max = mean_score + conf_int,
            conf_min = mean_score - conf_int,
            n = n
        )
    end
    info("Saving results in $output_dir")
    CSV.write(joinpath(output_dir, "results_table.csv"), grouped_info)
    println("\nQuick results:\n##############\n$venn_result\n$auprc_result\n$snr_result")
end
