
# This analysis generates several different common plots used to
# diagnose and explore screen data including: (1) volcano plot of gene
# pvalues and effect sizes (2) precision-recall curve for the hits
# (3) raw count plots
using Gadfly
using ColorBrewer

function gen_plots(filepath; debug=false, quiet=false)
    colors = palette("Accent", 4)

    for screentype in [FacsScreen(), GrowthScreen()]
        for crisprtype in [CRISPRi(), CRISPRn()]
            if typeof(screentype) == GrowthScreen
                max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
                    :inactive => (0.83, Delta(0.0)),
                    :negcontrol => (0.05, Delta(0.0)),
                    :increasing => (0.02, TruncatedNormal(0.1, 0.1, 0.025, 1)),
                    :decreasing => (0.1, TruncatedNormal(-0.55, 0.2, -1, -0.1))
                )
                lib = Library(max_phenotype_dists, crisprtype)
            else
                lib = Library(crisprtype)
            end

            # redefining function each time and abusing lexical scoping
            gen_plots = (bc_counts, genes) -> begin
                curr_dir = dirname(filepath)
                new_filename = joinpath(curr_dir, "$(basename(filepath))_$(typeof(crisprtype))_$(typeof(screentype))_counts.svg")
                nopseudo = bc_counts[(bc_counts[:counts_bin1] .> 0.5) .& (bc_counts[:counts_bin2] .> 0.5), :]
                draw(SVG(new_filename, 10cm, 10cm), plot(nopseudo,
                x=:counts_bin1, y=:counts_bin2, color=:class, Scale.x_log10,
                Scale.y_log10, Scale.color_discrete_manual(colors..., levels=sort(unique(nopseudo[:class]))),
                Theme(highlight_width=0pt),
                Coord.cartesian(fixed=true),
                Guide.title("$(typeof(crisprtype)) $(typeof(screentype))"),
                Guide.xlabel("log counts t0"), Guide.ylabel("log counts t1")))
                new_filename = joinpath(curr_dir, "$(basename(filepath))_$(typeof(crisprtype))_$(typeof(screentype))_volcano.svg")
                draw(SVG(new_filename, 10cm, 10cm), plot(genes,
                x=:mean_bin2_div_bin1, y=:pvalue_bin2_div_bin1, color=:class,
                Scale.color_discrete_manual(colors..., levels=sort(unique(genes[:class]))),
                Theme(highlight_width=0pt),
                Guide.title("$(typeof(crisprtype)) $(typeof(screentype))"),
                Guide.xlabel("mean log2 fold change"), Guide.ylabel("-log10 pvalue")))
                (1)
            end

            run_exp(screentype, lib, gen_plots; run_idx=1)
        end
    end
end
