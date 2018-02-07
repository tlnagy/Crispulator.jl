function differences_between_bins(raw_data::Associative{Symbol, DataFrame};
                                  first_bin=:bin1,
                                  last_bin=maximum(keys(raw_data)))

    for (bin, seq_data) in raw_data
        sort!(seq_data, cols=[:barcodeid])
        # add a pseudocount of 0.5 to every value to prevent -Inf's when
        # taking the log
        seq_data[:counts] += 0.5
        seq_data[:freqs] = seq_data[:counts]./sum(seq_data[:counts])
        # normalize to median of negative controls, fixes #19
        # TODO: consider normalizing the std dev
        negcontrol_freqs = seq_data[seq_data[:class] .== :negcontrol, :freqs]
        (length(negcontrol_freqs) == 0) && error("No negative control guides found. Try increasing "*
            "the frequency of negative controls or increase the number of genes.")
        med = median(negcontrol_freqs)
        seq_data[:rel_freqs] = seq_data[:freqs] ./ med
    end

    combined = copy(raw_data[first_bin])
    rename!(combined, Dict(:freqs => Symbol("freqs_", first_bin),
                           :counts => Symbol("counts_", first_bin),
                           :rel_freqs => Symbol("rel_freqs_", first_bin)))

    for (bin, seq_data) in raw_data
        (bin == first_bin) && continue

        combined[Symbol("freqs_", bin)] = seq_data[:freqs]
        combined[Symbol("counts_", bin)] = seq_data[:counts]
        combined[Symbol("rel_freqs_", bin)] = seq_data[:rel_freqs]
        combined[Symbol("log2fc_", bin)] =
        log2(combined[Symbol("rel_freqs_", bin)]./combined[Symbol("rel_freqs_", first_bin)])
    end

    nonnegs = combined[combined[:class] .!= :negcontrol, :]
    negcontrols = combined[combined[:class] .== :negcontrol, Symbol("log2fc_", last_bin)]

    genes = by(nonnegs, [:gene, :behavior, :class]) do barcodes
        log2fcs = barcodes[Symbol("log2fc_", last_bin)]
        result = MannWhitneyUTest(log2fcs, negcontrols)
        DataFrame(pvalue = -log10(pvalue(result)), mean= mean(log2fcs))
    end
    genes[:absmean] = abs(genes[:mean])
    genes[:pvalmeanprod] = genes[:mean] .* genes[:pvalue]

    combined, genes
end