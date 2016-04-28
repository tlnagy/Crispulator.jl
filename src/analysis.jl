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
        draw(PDF("plots/volcano_plot_by_behavior.pdf", 12cm, 10cm),
        plot(genes, x=:mean, y=:pvalue, color=:behavior, Guide.xlabel("mean log2 FC"),
        Guide.ylabel("-log10 pvalue"), Theme(highlight_width=0pt)))
        draw(PDF("plots/volcano_plot_by_class.pdf", 12cm, 10cm),
        plot(genes, x=:mean, y=:pvalue, color=:class, Theme(highlight_width=0pt)))
        draw(PDF("plots/distributions_of_observed_phenotypes.pdf", 12cm, 6cm),
        plot(combined, x=:obs_phenotype, color=:class, Geom.density))
        draw(PDF("plots/roc.pdf", 12cm, 10cm),
        plot(x=fprs, y=tprs, Geom.line, Coord.cartesian(fixed=true),
        Guide.xlabel("fpr"), Guide.ylabel("tpr")))
    end
    auroc
end
