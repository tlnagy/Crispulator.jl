# Custom Simulations

This is an example of building a custom simulation of a FACS-based screen.

!!! tip
    This section is complementary to the Implementation section of the
    [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1759-9#Sec2)

## Setup

```@setup 1
using Simulation
using Gadfly
using DataFrames
using Distributions
```

Lets first start the Julia REPL or a Julia session inside of a Jupyter Notebook
and load the packages we'll need:

```julia
using Gadfly
include(joinpath(Pkg.dir("Crispulator"), "src", "simulation", "load.jl"))
```

## Basic screen parameters

First lets design a simple [`Simulation.FacsScreen`](@ref) with 250 genes with
5 guides per gene. Lets say that we make sure we have 1000x as many cells as
guides during transfection (`representation`) and sorting
(`bottleneck_representation`) and 1000x as many reads as guides during sequencing
(`seq_depth`). We'll leave the rest of the values as the defaults and print
the object

```@example 1
s = FacsScreen()
s.num_genes = 250
s.coverage = 5
s.representation = 1000
s.bottleneck_representation = 1000
s.seq_depth = 1000
println(s)
```

## Construction of true phenotype distribution

Next, lets make our distribution of true phenotypes. The basic layout is a `Dict`
mapping a class name to a tuple of the probability of selecting this class and
then the [`Distributions.Sampleable`](https://juliastats.github.io/Distributions.jl/latest/types.html#Sampleable-1)
from which to draw a random phenotype from this class. The probabilities across
all the classes should add up to 1.

For example, here we are making three different classes of "genes": the
first group are `:inactive`, i.e. they have no phenotype, so we'll set their
phenotypes to 0.0 using a [`Simulation.Delta`](@ref). We'll also make them 60%
of all the genes. The second group are the negative controls `:negcontrol`
(the only required group) which make up 10% of the population of genes and also
have no effect. The final group is `:increasing` which makes up 30% of all
genes and which are represented by a Normal(μ=0.1, σ=0.1) distribution clamped
between 0.025 and 1.

```@example 1
max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
    :inactive => (0.60, Delta(0.0)),
    :negcontrol => (0.1, Delta(0.0)),
    :increasing => (0.3, TruncatedNormal(0.1, 0.1, 0.025, 1)),
);
```

!!! note

    The `:negcontrol` class needs to be present because Crispulator normalizes
    the frequencies of all other guides against the median frequency of the
    negative control guides. Also the distribution of `:negcontrol` guides serve
    as the null distribution against which the log2 fold changes of guides
    targeting a specific gene are assayed to calculate a statistical
    significance of the shift for each gene. See
    [`Simulation.differences_between_bins`](@ref) for more details.


## Library construction

Now, we actually build the library. Here we're making a
[`Simulation.CRISPRi`](@ref) library and then getting the guides that were built
from the true phenotype distribution that we constructed above and we also get
the frequency of each guide in the library.

```@example 1
lib = Library(max_phenotype_dists, CRISPRi())
guides, guide_freqs_dist = construct_library(s, lib);
```

Lets first look at what the true phenotype distribution of our different classes
of guides looks like

```@example 1
df = DataFrame(Dict(
    :phenotype=>map(x->x.theo_phenotype, guides),
    :class=>map(x->x.class, guides),
    :freq=>pdf.(guide_freqs_dist, 1:(length(guides)))
))
plot(df, x=:phenotype, color=:class, Geom.histogram, Guide.ylabel("Number of guides"),
Guide.title("Guide phenotype distribution"))
```

As you can see, most guides should have a phenotype of 0. In FACS Screens this is
equivalent to having no preference to being in either the left (`bin1`) or right
(`bin2`) bins. The `:increasing` genes have a small preference to be in the right
bin.

We can also look at the frequency of each guide in the library, which follows a
Log-Normal distribution.

```@example 1
plot(df, x=:freq, color=:class, Geom.histogram(position=:stack),
    Guide.xlabel("Frequency"), Guide.ylabel("Number of guides"),
    Guide.title("Frequencies of guides in simulated library"))
```

## Performing the screen

Now, we'll actually perform the screen. We'll first perform the transection via
[`Simulation.transfect`](@ref), followed by the selection process via
[`Simulation.select`](@ref):

```@example 1
cells, cell_phenotypes = transfect(s, lib, guides, guide_freqs_dist)
bin_cells = select(s, cells, cell_phenotypes, guides)
freqs = counts_to_freqs(bin_cells, length(guides));
```

Lets look at what the observed phenotype distribution looks like when the
selection was performed:

```@example 1
df = DataFrame(Dict(
    :phenotype=>map(x->x.theo_phenotype, guides),
    :class=>map(x->x.class, guides),
    :obs_freq=>map(x->x.obs_phenotype, guides)
))
plot(df, x=:obs_freq, Geom.density, Guide.xlabel("Observed phenotype on FACS machine"),
Guide.title("Kernel density estimate of guide observed phenotypes"), Guide.ylabel("ρ"))
```

As you can see, this looks like many FACS plots, e.g. when looking at density
along the fluorescence channel. A quick sanity check is that we should see a
slight enrichment of the frequency of `:increasing` genes on the right side

```@example 1
plot(df, x=:obs_freq, color=:class, Geom.density, Guide.xlabel("Observed phenotype on FACS machine"),
Guide.title("Kernel density estimate of guide observed phenotypes"), Guide.ylabel("ρ"))
```

And that is what we see. The change is really small (this is pretty usual), but
the later analysis will be able to pull out the `increasing` genes.

## Sequencing and Analysis

Now we'll use [`Simulation.sequencing`](@ref) to simulate sequencing by
transforming the guide frequencies into a Categorical distribution and drawing
a random sample of reads from this distribution. Finally, we'll use
the [`Simulation.differences_between_bins`](@ref) function to compute the
differences between bins on a per-guide level (`guide_data`) and per-gene level
(`gene_data`).

```@example 1
raw_data = sequencing(Dict(:bin1=>s.seq_depth, :bin2=>s.seq_depth), guides, freqs)
guide_data, gene_data = differences_between_bins(raw_data);
```

Here's what the per-guide data looks like:

```@example 1
head(guide_data) # hide
```

!!! tip
    See [`Simulation.differences_between_bins`](@ref) for details on what each
    column means.

And the gene level data

```@example 1
head(gene_data) # hide
```

We can generate standard pooled screen plots from this dataset. Like a count
scatterplot:

```@example 1
nopseudo = guide_data[(guide_data[:counts_bin1] .> 0.5) .& (guide_data[:counts_bin2] .> 0.5), :]
plot(nopseudo, x=:counts_bin1, y=:counts_bin2, color=:class, Scale.x_log10,
Scale.y_log10, Theme(highlight_width=0pt), Coord.cartesian(fixed=true),
Guide.xlabel("log counts bin1"), Guide.ylabel("log counts bin2"))
```

And a volcano plot:

```@example 1
plot(gene_data, x=:mean_bin2_div_bin1, y=:pvalue_bin2_div_bin1, color=:class, Theme(highlight_width=0pt),
Guide.xlabel("mean log2 fold change"), Guide.ylabel("-log10 pvalue"))
```

And finally we can see how well we can differentiate between the different
classes using Area Under the Precision-Recall Curve ([`Simulation.auprc`](@ref))

```@example 1
auprc(gene_data[:pvalmeanprod_bin2_div_bin1], gene_data[:class], Set([:increasing]))[1]
```

[`Simulation.auroc`](@ref) and [`Simulation.venn`](@ref) are also good summary
statistics.
