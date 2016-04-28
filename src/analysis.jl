function analyze(raw_data::Vector{DataFrame}; gen_plots=false)
    map(x -> sort!(x, cols=[:barcodeid]), raw_data);
    map(x -> x[:freqs] = x[:counts]./sum(x[:counts]), raw_data);

    combined = copy(raw_data[1])
    rename!(combined, Dict(:freqs => :freqs1))
    combined[:freqs2] = raw_data[2][:freqs]
    combined[:log2fc] = log2(combined[:freqs1]./combined[:freqs2])

    nonnegs = combined[combined[:class] .!= :negcontrol, :]
    negcontrols = combined[combined[:class] .== :negcontrol, :log2fc]

    if gen_plots
        df = by(nonnegs, [:gene, :behavior, :class]) do barcodes
            result = MannWhitneyUTest(barcodes[:log2fc], negcontrols)
            DataFrame(pvalue = -log10(pvalue(result)), mean= mean(barcodes[:log2fc]))
        end
        draw(PDF("plots/volcano_plot_by_behavior.pdf", 12cm, 10cm),
        plot(df, x=:mean, y=:pvalue, color=:behavior, Guide.xlabel("mean log2 FC"),
        Guide.ylabel("-log10 pvalue"), Theme(highlight_width=0pt)))
        draw(PDF("plots/volcano_plot_by_class.pdf", 12cm, 10cm),
        plot(df, x=:mean, y=:pvalue, color=:class, Theme(highlight_width=0pt)))
        draw(PDF("plots/distributions_of_observed_phenotypes.pdf", 12cm, 6cm),
        plot(combined, x=:obs_phenotype, color=:class, Geom.density))
    end
end
