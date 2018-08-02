"""
$(SIGNATURES)

Given the raw data from [`Simulation.sequencing`](@ref) returns two DataFrames

1. `guide_data`: This DataFrame contains the per-guide level data including the
    log2 fold change in the normalized frequencies of each guide between each
    pairwise combination of bins. Thus, if there are `n` bins, then it computes
    the log2 fold changes for the ``\\frac{n!}{2(n-2)!}`` combinations

2. `gene_data`: This DataFrame contains the same information but grouped by
    gene. The log2 fold change data from the first DataFrame is used to calculate
    the average log2 fold change per gene and a pvalue computed using a
    [Mann-Whitney U-test](https://en.wikipedia.org/wiki/Mann-Whitney_U_test) as
    measure of how consistently shifted the guides are of this gene versus the
    population of negative control guides. (see below for more info)

A typical 2 bin `guide_data` DataFrame contains the following columns:

| Column Name  | Meaning  |
|:--|:--|
| `gene` | the gene ID of that this guide targets|
| `knockdown` | activity of the guide on 0 to 1 scale, where 1 is complete knockout|
| `barcodeid` | the ID of this specific guide|
| `theo_phenotype` | expected phenotype of this guide, generally a -1 to 1 scale|
| `behavior` | whether the target gene displays a linear or sigmoidal response to incomplete knockdown (see [`Simulation.Library`](@ref) for more details)|
| `class` | which phenotype distribution the target gene was drawn from (see [`Simulation.Library`](@ref) for more details). Serves as the "ground truth" label against which screen performance is evaluated, e.g. with [`Simulation.auprc`](@ref) |
| `initial_freq` | frequency of guide post-transfection (see [`Simulation.transfect`](@ref))|
| `counts_bin1` | the raw number of reads for each guide in the first bin|
| `freqs_bin1` | the number of reads for each guide divided by the total number of reads in this bin|
| `rel_freqs_bin1` | the frequency of each guide divided by the median frequency of negative control guides|
| `counts_bin2` | the raw number of reads for each guide in the second bin|
| `freqs_bin2` | the number of reads for each guide divided by the total number of reads in this bin|
| `rel_freqs_bin2` | the frequency of each guide divided by the median frequency of negative control guides for this bin|
| `log2fc_bin2_div_bin1` | the log2 fold change in relative guide frequencies between `bin2` and `bin1`, equivalent to `log2(rel_freqs_bin2/rel_freqs_bin1)` |

A typical `gene_data` DataFrame contains the following data:

| Column Name  | Meaning  |
|:--|:--|
| `gene` | this gene's ID|
| `class` | see above |
| `behavior` | see above |
| `mean_bin2_div_bin1` | the mean log 2 fold change in relative frequencies between from `bin1` to `bin2` for all the guides targeting this gene. Calculated as ``\\frac{1}{k}\\sum_k \\text{log2fc_bin2_div_bin1}_k`` for the ``k`` guides targeting each gene |
| `pvalue_bin2_div_bin1` | the -log10 pvalue of the log2 fold changes of all guides targeting this gene as computed by the non-parametric Mann-Whitney U-test. A measure of the consistency of the log 2 fold changes[^1]|
| `absmean_bin2_div_bin1` | absolute value of `mean_bin2_div_bin1` per-gene |
| `pvalmeanprod_bin2_div_bin1` | `mean_bin2_div_bin1` multiplied with the `pvalue_bin2_div_bin1` per-gene|

# Further reading

[^1]: Kampmann M, Bassik MC, Weissman JS. Integrated platform for genome-wide screening
    and construction of high-density genetic interaction maps in mammalian cells.
    *Proc Natl Acad Sci U S A*. 2013;110:E2317–26.

"""
function differences_between_bins(raw_data::Associative{Symbol, DataFrame})
    for (bin, seq_data) in raw_data
        sort!(seq_data, [:barcodeid])
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

    bin1 = first(raw_data)[2]
    cols_to_copy = [col for col in names(bin1) if !(col in (:counts, :freqs, :rel_freqs))]
    guide_data = bin1[cols_to_copy]
    gene_data = DataFrame()

    # processed bins; so we don't re-add bins when doing all combos
    proc_bins = Set{Symbol}()

    # compute pairwise log 2 fold changes between every combination of bins
    for bins in combinations(collect(keys(raw_data)), 2)
        for bin in bins
            (bin in proc_bins) && continue
            push!(proc_bins, bin)
            guide_data[[Symbol("counts_", bin), Symbol("freqs_", bin), Symbol("rel_freqs_", bin)]] =
                raw_data[bin][[:counts, :freqs, :rel_freqs]]
        end
        suffix = "_$(bins[2])_div_$(bins[1])"
        log2fc_col = Symbol(:log2fc, suffix)
        guide_data[log2fc_col] = log2.(guide_data[Symbol(:rel_freqs_, bins[2])]
                                        ./guide_data[Symbol(:rel_freqs_, bins[1])])

        nonnegs = guide_data[guide_data[:class] .!= :negcontrol, :]
        negcontrols = guide_data[guide_data[:class] .== :negcontrol, log2fc_col]

        tmp = by(nonnegs, [:gene, :behavior, :class]) do barcodes
            log2fcs = barcodes[log2fc_col]
            result = MannWhitneyUTest(log2fcs, negcontrols)
            DataFrame(pvalue = -log10(pvalue(result)), mean= mean(log2fcs))
        end
        gene_data[[:gene, :behavior, :class, Symbol(:pvalue, suffix), Symbol(:mean, suffix)]] = tmp[[:gene, :behavior, :class, :pvalue, :mean]]
        gene_data[Symbol(:absmean, suffix)] = abs.(gene_data[Symbol(:mean, suffix)])
        gene_data[Symbol(:pvalmeanprod, suffix)] = gene_data[Symbol(:mean, suffix)] .* gene_data[Symbol(:pvalue, suffix)]

    end

    guide_data, gene_data
end
