function analyze(raw_data::Dict{Symbol, DataFrame}; gen_plots=false)
    for (bin, seq_data) in raw_data
        sort!(seq_data, cols=[:barcodeid])
        # add a pseudocount of 0.5 to every value to prevent -Inf's when
        # taking the log
        seq_data[:counts] += 0.5
        seq_data[:freqs] = seq_data[:counts]./sum(seq_data[:counts])
    end
    @assert length(keys(raw_data)) == 2 "exactly two bins needed"

    combined = copy(raw_data[:bin1])
    rename!(combined, Dict(:freqs => :freqs1))
    combined[:freqs2] = raw_data[:bin2][:freqs]
    combined[:log2fc] = log2(combined[:freqs1]./combined[:freqs2])

    nonnegs = combined[combined[:class] .!= :negcontrol, :]
    negcontrols = combined[combined[:class] .== :negcontrol, :log2fc]

    genes = by(nonnegs, [:gene, :behavior, :class]) do barcodes
        result = MannWhitneyUTest(barcodes[:log2fc], negcontrols)
        DataFrame(pvalue = -log10(pvalue(result)), mean= mean(barcodes[:log2fc]))
    end
    genes[:absmean] = abs(genes[:mean])
    genes[:pvalmeanprod] = genes[:absmean] .* genes[:pvalue]

    auroc = compute_auroc(Vector(genes[:pvalmeanprod]), Vector(genes[:class]), 50)

    if gen_plots
        draw(PNG("plots/volcano_plot_by_behavior.png", 12cm, 10cm, dpi=300),
        plot(genes, x=:mean, y=:pvalue, color=:behavior, Guide.xlabel("mean log2 FC"),
        Guide.ylabel("-log10 pvalue"), Theme(highlight_width=0pt)))
        draw(PNG("plots/volcano_plot_by_class.png", 12cm, 10cm, dpi=300),
        plot(genes, x=:mean, y=:pvalue, color=:class, Theme(highlight_width=0pt)))

        # remove guides that weren't observed during sorting
        filtered = combined[isfinite(combined[:obs_phenotype]), :]

        draw(PNG("plots/distributions_of_observed_phenotypes.png", 12cm, 6cm, dpi=300),
        plot(filtered, x=:obs_phenotype, color=:class, Geom.density))
        draw(PNG("plots/freq_scatter.png", 12cm, 10cm, dpi=300),
        plot(filtered, x=:freqs1, y=:freqs2, color=:class, Scale.x_log10, Scale.y_log10,
        Theme(highlight_width=0pt)))

        # ROC plots

        _, s_tprs, s_fprs = compute_roc(genes[genes[:class] .== :sigmoid, :], :pvalmeanprod, 50)
        _, l_tprs, l_fprs = compute_roc(genes[genes[:class] .== :linear, :], :pvalmeanprod, 50)

        data = DataFrame(tprs=vcat(s_tprs, l_tprs), fprs = vcat(s_fprs, l_fprs),
               class=vcat(rep(:sigmoid,50), rep(:linear, 50)))

        draw(PNG("plots/roc.png", 12cm, 10cm, dpi=300),
        plot(data, x=fprs, y=tprs, color=class, Geom.line, Coord.cartesian(fixed=true),
        Guide.xlabel("fpr"), Guide.ylabel("tpr")))
    end
    auroc
end
