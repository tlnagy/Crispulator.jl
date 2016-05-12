function difference_between_two_bins(raw_data::Dict{Symbol, DataFrame})
    for (bin, seq_data) in raw_data
        sort!(seq_data, cols=[:barcodeid])
        # add a pseudocount of 0.5 to every value to prevent -Inf's when
        # taking the log
        seq_data[:counts] += 0.5
        seq_data[:freqs] = seq_data[:counts]./sum(seq_data[:counts])
    end
    @assert length(keys(raw_data)) == 2 "exactly two bins needed"

    combined = copy(raw_data[:bin1])
    rename!(combined, Dict(:freqs => :freqs1, :counts => :counts1))
    combined[:freqs2] = raw_data[:bin2][:freqs]
    combined[:counts2] = raw_data[:bin2][:counts]
    combined[:log2fc] = log2(combined[:freqs2]./combined[:freqs1])

    nonnegs = combined[combined[:class] .!= :negcontrol, :]
    negcontrols = combined[combined[:class] .== :negcontrol, :log2fc]

    genes = by(nonnegs, [:gene, :behavior, :class]) do barcodes
        result = MannWhitneyUTest(barcodes[:log2fc], negcontrols)
        DataFrame(pvalue = -log10(pvalue(result)), mean= mean(barcodes[:log2fc]))
    end
    genes[:absmean] = abs(genes[:mean])
    genes[:pvalmeanprod] = genes[:absmean] .* genes[:pvalue]

    genes
end
