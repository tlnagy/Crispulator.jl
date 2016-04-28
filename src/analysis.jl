function analyze(raw_data::Vector{DataFrame}; gen_plots=false)
    map(x -> sort!(x, cols=[:barcodeid]), raw_data);
    # add a pseudocount of 0.5 to every value to prevent -Inf's when
    # taking the log
    map(x -> x[:counts] += 0.5, raw_data);
    map(x -> x[:freqs] = x[:counts]./sum(x[:counts]), raw_data);

    combined = copy(raw_data[1])
    rename!(combined, Dict(:freqs => :freqs1))
    combined[:freqs2] = raw_data[2][:freqs]
    combined[:log2fc] = log2(combined[:freqs1]./combined[:freqs2])

    nonnegs = combined[combined[:class] .!= :negcontrol, :]
    negcontrols = combined[combined[:class] .== :negcontrol, :log2fc]

    genes = by(nonnegs, [:gene, :behavior, :class]) do barcodes
        result = MannWhitneyUTest(barcodes[:log2fc], negcontrols)
        DataFrame(pvalue = -log10(pvalue(result)), mean= mean(barcodes[:log2fc]))
    end
    genes[:absmean] = abs(genes[:mean])
    genes[:pvalmeanprod] = genes[:absmean] .* genes[:pvalue]

    auroc, tprs, fprs = compute_roc(genes, :pvalmeanprod, 50)

    if gen_plots
        draw(PNG("plots/volcano_plot_by_behavior.png", 12cm, 10cm, dpi=300),
        plot(genes, x=:mean, y=:pvalue, color=:behavior, Guide.xlabel("mean log2 FC"),
        Guide.ylabel("-log10 pvalue"), Theme(highlight_width=0pt)))
        draw(PNG("plots/volcano_plot_by_class.png", 12cm, 10cm, dpi=300),
        plot(genes, x=:mean, y=:pvalue, color=:class, Theme(highlight_width=0pt)))
        draw(PNG("plots/distributions_of_observed_phenotypes.png", 12cm, 6cm, dpi=300),
        plot(combined, x=:obs_phenotype, color=:class, Geom.density))
        draw(PNG("plots/roc.png", 12cm, 10cm, dpi=300),
        plot(x=fprs, y=tprs, Geom.line, Coord.cartesian(fixed=true),
        Guide.xlabel("fpr"), Guide.ylabel("tpr")))
    end
    auroc
end
