var documenterSearchIndex = {"docs":
[{"location":"custom/#Custom-Simulations","page":"Custom Simulations","title":"Custom Simulations","text":"","category":"section"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"This is an example of building a custom simulation of a FACS-based screen.","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"tip: Tip\nThis section is complementary to the Implementation section of the paper","category":"page"},{"location":"custom/#Setup","page":"Custom Simulations","title":"Setup","text":"","category":"section"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"using Crispulator\nusing Gadfly\nusing DataFrames\nusing Distributions","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"Lets first start the Julia REPL or a Julia session inside of a Jupyter Notebook and load the packages we'll need:","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"using Crispulator\nusing Gadfly","category":"page"},{"location":"custom/#Basic-screen-parameters","page":"Custom Simulations","title":"Basic screen parameters","text":"","category":"section"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"First lets design a simple Crispulator.FacsScreen with 250 genes with 5 guides per gene. Lets say that we make sure we have 1000x as many cells as guides during transfection (representation) and sorting (bottleneck_representation) and 1000x as many reads as guides during sequencing (seq_depth). We'll leave the rest of the values as the defaults and print the object","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"s = FacsScreen()\ns.num_genes = 250\ns.coverage = 5\ns.representation = 1000\ns.bottleneck_representation = 1000\ns.seq_depth = 1000\nprintln(s)","category":"page"},{"location":"custom/#Construction-of-true-phenotype-distribution","page":"Custom Simulations","title":"Construction of true phenotype distribution","text":"","category":"section"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"Next, lets make our distribution of true phenotypes. The basic layout is a Dict mapping a class name to a tuple of the probability of selecting this class and then the Distributions.Sampleable from which to draw a random phenotype from this class. The probabilities across all the classes should add up to 1.","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"For example, here we are making three different classes of \"genes\": the first group are :inactive, i.e. they have no phenotype, so we'll set their phenotypes to 0.0 using a Crispulator.Delta. We'll also make them 60% of all the genes. The second group are the negative controls :negcontrol (the only required group) which make up 10% of the population of genes and also have no effect. The final group is :increasing which makes up 30% of all genes and which are represented by a Normal(μ=0.1, σ=0.1) distribution clamped between 0.025 and 1.","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(\n    :inactive => (0.60, Delta(0.0)),\n    :negcontrol => (0.1, Delta(0.0)),\n    :increasing => (0.3, truncated(Normal(0.1, 0.1), 0.025, 1)),\n);","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"note: Note\nThe :negcontrol class needs to be present because Crispulator normalizes the frequencies of all other guides against the median frequency of the negative control guides. Also the distribution of :negcontrol guides serve as the null distribution against which the log2 fold changes of guides targeting a specific gene are assayed to calculate a statistical significance of the shift for each gene. See Crispulator.differences_between_bins for more details.","category":"page"},{"location":"custom/#Library-construction","page":"Custom Simulations","title":"Library construction","text":"","category":"section"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"Now, we actually build the library. Here we're making a Crispulator.CRISPRi library and then getting the guides that were built from the true phenotype distribution that we constructed above and we also get the frequency of each guide in the library.","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"lib = Library(max_phenotype_dists, CRISPRi())\nguides, guide_freqs_dist = construct_library(s, lib);","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"Lets first look at what the true phenotype distribution of our different classes of guides looks like","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"df = DataFrame(Dict(\n    :phenotype=>map(x->x.theo_phenotype, guides),\n    :class=>map(x->x.class, guides),\n    :freq=>pdf.(guide_freqs_dist, 1:(length(guides)))\n))\nplot(df, x=:phenotype, color=:class, Geom.histogram, Guide.ylabel(\"Number of guides\"),\nGuide.title(\"Guide phenotype distribution\"))","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"As you can see, most guides should have a phenotype of 0. In FACS Screens this is equivalent to having no preference to being in either the left (bin1) or right (bin2) bins. The :increasing genes have a small preference to be in the right bin.","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"We can also look at the frequency of each guide in the library, which follows a Log-Normal distribution.","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"plot(df, x=:freq, color=:class, Geom.histogram(position=:stack),\n    Guide.xlabel(\"Frequency\"), Guide.ylabel(\"Number of guides\"),\n    Guide.title(\"Frequencies of guides in simulated library\"))","category":"page"},{"location":"custom/#Performing-the-screen","page":"Custom Simulations","title":"Performing the screen","text":"","category":"section"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"Now, we'll actually perform the screen. We'll first perform the transection via Crispulator.transfect, followed by the selection process via Crispulator.select:","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"cells, cell_phenotypes = transfect(s, lib, guides, guide_freqs_dist)\nbin_cells = Crispulator.select(s, cells, cell_phenotypes, guides)\nfreqs = counts_to_freqs(bin_cells, length(guides));","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"Lets look at what the observed phenotype distribution looks like when the selection was performed:","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"df = DataFrame(Dict(\n    :phenotype=>map(x->x.theo_phenotype, guides),\n    :class=>map(x->x.class, guides),\n    :obs_freq=>map(x->x.obs_phenotype, guides)\n))\nplot(df, x=:obs_freq, Geom.density, Guide.xlabel(\"Observed phenotype on FACS machine\"),\nGuide.title(\"Kernel density estimate of guide observed phenotypes\"), Guide.ylabel(\"ρ\"))","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"As you can see, this looks like many FACS plots, e.g. when looking at density along the fluorescence channel. A quick sanity check is that we should see a slight enrichment of the frequency of :increasing genes on the right side","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"plot(df, x=:obs_freq, color=:class, Geom.density, Guide.xlabel(\"Observed phenotype on FACS machine\"),\nGuide.title(\"Kernel density estimate of guide observed phenotypes\"), Guide.ylabel(\"ρ\"))","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"And that is what we see. The change is really small (this is pretty usual), but the later analysis will be able to pull out the increasing genes.","category":"page"},{"location":"custom/#Sequencing-and-Analysis","page":"Custom Simulations","title":"Sequencing and Analysis","text":"","category":"section"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"Now we'll use Crispulator.sequencing to simulate sequencing by transforming the guide frequencies into a Categorical distribution and drawing a random sample of reads from this distribution. Finally, we'll use the Crispulator.differences_between_bins function to compute the differences between bins on a per-guide level (guide_data) and per-gene level (gene_data).","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"raw_data = sequencing(Dict(:bin1=>s.seq_depth, :bin2=>s.seq_depth), guides, freqs)\nguide_data, gene_data = differences_between_bins(raw_data);","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"Here's what the per-guide data looks like:","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"first(guide_data, 10) # hide","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"tip: Tip\nSee Crispulator.differences_between_bins for details on what each column means.","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"And the gene level data","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"first(gene_data, 10) # hide","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"We can generate standard pooled screen plots from this dataset. Like a count scatterplot:","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"nopseudo = guide_data[(guide_data[!, :counts_bin1] .> 0.5) .& (guide_data[!, :counts_bin2] .> 0.5), :]\nplot(nopseudo, x=:counts_bin1, y=:counts_bin2, color=:class, Scale.x_log10,\nScale.y_log10, Theme(highlight_width=0pt), Coord.cartesian(fixed=true),\nGuide.xlabel(\"log counts bin1\"), Guide.ylabel(\"log counts bin2\"))","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"And a volcano plot:","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"plot(gene_data, x=:mean_bin2_div_bin1, y=:pvalue_bin2_div_bin1, color=:class, Theme(highlight_width=0pt),\nGuide.xlabel(\"mean log2 fold change\"), Guide.ylabel(\"-log10 pvalue\"))","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"And finally we can see how well we can differentiate between the different classes using Area Under the Precision-Recall Curve (Crispulator.auprc)","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"auprc(gene_data[!, :pvalmeanprod_bin2_div_bin1], gene_data[!, :class], Set([:increasing]))[1]","category":"page"},{"location":"custom/","page":"Custom Simulations","title":"Custom Simulations","text":"Crispulator.auroc and Crispulator.venn are also good summary statistics.","category":"page"},{"location":"#CRISPulator","page":"Home","title":"CRISPulator","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A pooled genetic screen simulator","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pooled screens are very useful, but it is difficult to select all of the necessary parameters a priori. This work aims to explore the importance of various screen parameters and the relevance of design choices on the downstream analysis and interpretation.","category":"page"},{"location":"#Setup","page":"Home","title":"Setup","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To run the simulation, you will need a recent version of Julia installed and in your PATH. Then navigate into the root directory of the project and run julia. Run the following command:","category":"page"},{"location":"","page":"Home","title":"Home","text":"$ julia -e 'Pkg.update(); Pkg.add(\"Crispulator\")'","category":"page"},{"location":"","page":"Home","title":"Home","text":"this copies Crispulator over to the Julia package directory and installs all of its dependencies. Note the $ just represents the prompt, you don't need to type it. Make sure to run the tests and verify that everything is passing","category":"page"},{"location":"","page":"Home","title":"Home","text":"$ julia -e 'Pkg.test(\"Crispulator\")'","category":"page"},{"location":"#Quickstart","page":"Home","title":"Quickstart","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For most simple cases, no writing of Julia code is necessary. See the Tutorial for using Crispulator in this manner.","category":"page"},{"location":"#Advanced","page":"Home","title":"Advanced","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"See the Custom Simulations for a step-by-step guide to writing a custom simulation on top of Crispulator. Additionally, the Crispulator Internals page has in-depth documentation needed for more advanced usage","category":"page"},{"location":"internals/#Crispulator-Internals","page":"Crispulator Internals","title":"Crispulator Internals","text":"","category":"section"},{"location":"internals/","page":"Crispulator Internals","title":"Crispulator Internals","text":"More in-depth documentation on specific types and functions.","category":"page"},{"location":"internals/#Contents","page":"Crispulator Internals","title":"Contents","text":"","category":"section"},{"location":"internals/","page":"Crispulator Internals","title":"Crispulator Internals","text":"Pages = [\"internals.md\"]","category":"page"},{"location":"internals/#Index","page":"Crispulator Internals","title":"Index","text":"","category":"section"},{"location":"internals/","page":"Crispulator Internals","title":"Crispulator Internals","text":"Pages = [\"internals.md\"]","category":"page"},{"location":"internals/#Simulation-Types","page":"Crispulator Internals","title":"Simulation Types","text":"","category":"section"},{"location":"internals/","page":"Crispulator Internals","title":"Crispulator Internals","text":"Crispulator.FacsScreen\nCrispulator.GrowthScreen","category":"page"},{"location":"internals/#Crispulator.FacsScreen","page":"Crispulator Internals","title":"Crispulator.FacsScreen","text":"A type representing the parameters used in a typical FACS screen.\n\nnum_genes\nNumber of genes targeted by the screen\ncoverage\nNumber of guides per gene\nrepresentation\nNumber of cells with each guide. representation=10 means that there are 10 times as many cells as guides during transfection. The actual number of cells per guide post-transfection will be less depending on the MOI\nmoi\nThe multiplicity of infection, lambda, of the screen. We model this as a Poisson process during transfection (see Crispulator.transfect).\nnote: Note\nWe do not model multiple infections. We assume that the MOI is properly selected and less than half the cells are transfected by any virus, i.e. lambda lt 05 and then select only the cells that have a single transfection occurrence:P(x = 1 Poisson(λ))For lambda = 025 this works out to being approx 195 of the number of cells (num_genes * coverage * representation).\n\nσ\nThe standard deviation expected for cells during FACS sorting. This should be set according to the biological variance experimentally observed, e.g. in the fluorescence intensity of isogenic cells\nbin_info\nRange of guide phenotypes to collect in each bin\nExample\nIn the following example\np = 0.05\nbin_info = Dict(:bin1 => (0.0, p), :bin2 => (1.0-p, 1.0))\nThe 5th percentile of cells sorted according to their phenotype (fluorescence, size, etc) will be compared to the 95th percentile.\n\nbottleneck_representation\nNumber of cells sorted expressed as an integer multiple of the number of guides\nseq_depth\nSequencing depth as a integer multiple of the number of guides, i.e. seq_depth=10 is equivalent to 10 * num_genes * coverage reads.\n\n\n\n\n\n","category":"type"},{"location":"internals/#Crispulator.GrowthScreen","page":"Crispulator Internals","title":"Crispulator.GrowthScreen","text":"A type representing the parameters used in a typical growth-based screen.\n\nnum_genes\nNumber of genes targeted by the screen\ncoverage\nNumber of guides per gene\nrepresentation\nNumber of cells with each guide. representation=10 means that there are 10 times as many cells as guides during transfection. The actual number of cells per guide post-transfection will be less depending on the MOI\nmoi\nThe multiplicity of infection, lambda, of the screen. We model this as a Poisson process during transfection (see Crispulator.transfect).\nnote: Note\nWe do not model multiple infections. We assume that the MOI is properly selected and less than half the cells are transfected by any virus, i.e. lambda lt 05 and then select only the cells that have a single transfection occurrence:P(x = 1 Poisson(λ))For lambda = 025 this works out to being approx 195 of the number of cells (num_genes * coverage * representation).\n\nseq_depth\nSequencing depth as a integer multiple of the number of guides, i.e. seq_depth=10 is equivalent to 10 * num_genes * coverage reads.\n\nbottleneck_representation\nFor growth screens, how much of a bottleneck is applied. This is minimum number of cells that is kept when passaging the pool of cells. This is expressed as an integer multiple of the number of guides.\nnum_bottlenecks\nFor growth screens, how many bottlenecks are applied. This is the integer number of passages during the growth screen.\nnoise\nBefore each passage, the theoretical phenotype of each cell is convolved with a normal noise distribution with a standard deviation, σ, dictated by this parameter. This should be set based on an expectation of noisiness of the subsampling\n\n\n\n\n\n","category":"type"},{"location":"internals/#Key-functions","page":"Crispulator Internals","title":"Key functions","text":"","category":"section"},{"location":"internals/","page":"Crispulator Internals","title":"Crispulator Internals","text":"Crispulator.construct_library\nCrispulator.transfect\nCrispulator.select\nCrispulator.sequencing\nCrispulator.counts_to_freqs\nCrispulator.differences_between_bins","category":"page"},{"location":"internals/#Crispulator.construct_library","page":"Crispulator Internals","title":"Crispulator.construct_library","text":"construct_library(setup, lib)\n\n\nConstructs the guide library for N genes with coverage number of guides per gene. Returns a tuple of guides and their relative frequencies (assigned randomly).\n\n\n\n\n\n","category":"function"},{"location":"internals/#Crispulator.transfect","page":"Crispulator Internals","title":"Crispulator.transfect","text":"transfect(setup, lib, guides, guide_freqs_dist)\n\n\nSimulates the transfection the of the library of guides (guides) according to the parameters in Crispulator.Library lib and the frequencies represented by the Categorical distribution guide_freqs_dist.\n\nThis function creates a population of cells post-transfection where each cell has single guide RNA present according to the frequencies of the guides in the present in the library. The cells are then grown up so that there are setup.bottleneck_representation * num_guides of them. For FACS screens we assume a simple linear expansion of the cells, while for Growth screens the cells grow according as a function of their phenotype and degree of knockdown.\n\n\n\n\n\n","category":"function"},{"location":"internals/#Crispulator.select","page":"Crispulator Internals","title":"Crispulator.select","text":"select(setup, cells, cell_phenotypes, guides; debug)\n\n\nGiven cells, a vector of integers, and guides, a vector of barcodes, performs simulated FACS sorting of the cells into bins with the given cutoffs. σ refers to the spread of the phenotype during the FACS screen.\n\n\n\n\n\nselect(setup, initial_cells, initial_cell_phenotypes, guides; debug)\n\n\nGrowth Screen selection\n\n\n\n\n\n","category":"function"},{"location":"internals/#Crispulator.sequencing","page":"Crispulator Internals","title":"Crispulator.sequencing","text":"sequencing(depths, guides, samples)\n\n\nSimulates next-gen sequencing by generating a Categorical distribution based on frequencies of each guide in each bin and then randomly samples from this distributions up to the depth provided in depths for each bin. Returns a dict of DataFrames with the bin names as the keys. This object can then be passed through Crispulator.counts_to_freqs followed by Crispulator.differences_between_bins.\n\n\n\n\n\n","category":"function"},{"location":"internals/#Crispulator.counts_to_freqs","page":"Crispulator Internals","title":"Crispulator.counts_to_freqs","text":"counts_to_freqs(bin_cells, guide_count)\n\n\nConverts raw counts to frequencies by dividing by the total number of reads for each sample.\n\n\n\n\n\n","category":"function"},{"location":"internals/#Crispulator.differences_between_bins","page":"Crispulator Internals","title":"Crispulator.differences_between_bins","text":"differences_between_bins(raw_data)\n\n\nGiven the raw data from Crispulator.sequencing returns two DataFrames\n\nguide_data: This DataFrame contains the per-guide level data including the  log2 fold change in the normalized frequencies of each guide between each  pairwise combination of bins. Thus, if there are n bins, then it computes  the log2 fold changes for the fracn2(n-2) combinations\ngene_data: This DataFrame contains the same information but grouped by  gene. The log2 fold change data from the first DataFrame is used to calculate  the average log2 fold change per gene and a pvalue computed using a  Mann-Whitney U-test as  measure of how consistently shifted the guides are of this gene versus the  population of negative control guides. (see below for more info)\n\nA typical 2 bin guide_data DataFrame contains the following columns:\n\nColumn Name Meaning\ngene the gene ID of that this guide targets\nknockdown activity of the guide on 0 to 1 scale, where 1 is complete knockout\nbarcodeid the ID of this specific guide\ntheo_phenotype expected phenotype of this guide, generally a -1 to 1 scale\nbehavior whether the target gene displays a linear or sigmoidal response to incomplete knockdown (see Crispulator.Library for more details)\nclass which phenotype distribution the target gene was drawn from (see Crispulator.Library for more details). Serves as the \"ground truth\" label against which screen performance is evaluated, e.g. with Crispulator.auprc\ninitial_freq frequency of guide post-transfection (see Crispulator.transfect)\ncounts_bin1 the raw number of reads for each guide in the first bin\nfreqs_bin1 the number of reads for each guide divided by the total number of reads in this bin\nrel_freqs_bin1 the frequency of each guide divided by the median frequency of negative control guides\ncounts_bin2 the raw number of reads for each guide in the second bin\nfreqs_bin2 the number of reads for each guide divided by the total number of reads in this bin\nrel_freqs_bin2 the frequency of each guide divided by the median frequency of negative control guides for this bin\nlog2fc_bin2_div_bin1 the log2 fold change in relative guide frequencies between bin2 and bin1, equivalent to log2(rel_freqs_bin2/rel_freqs_bin1)\n\nA typical gene_data DataFrame contains the following data:\n\nColumn Name Meaning\ngene this gene's ID\nclass see above\nbehavior see above\nmean_bin2_div_bin1 the mean log 2 fold change in relative frequencies between from bin1 to bin2 for all the guides targeting this gene. Calculated as frac1ksum_k textlog2fc_bin2_div_bin1_k for the k guides targeting each gene\npvalue_bin2_div_bin1 the -log10 pvalue of the log2 fold changes of all guides targeting this gene as computed by the non-parametric Mann-Whitney U-test. A measure of the consistency of the log 2 fold changes[1]\nabsmean_bin2_div_bin1 absolute value of mean_bin2_div_bin1 per-gene\npvalmeanprod_bin2_div_bin1 mean_bin2_div_bin1 multiplied with the pvalue_bin2_div_bin1 per-gene\n\nFurther reading\n\n[1]: Kampmann M, Bassik MC, Weissman JS. Integrated platform for genome-wide screening and construction of high-density genetic interaction maps in mammalian cells. Proc Natl Acad Sci U S A. 2013;110:E2317–26.\n\n\n\n\n\n","category":"function"},{"location":"internals/#Miscellaneous-Types","page":"Crispulator Internals","title":"Miscellaneous Types","text":"","category":"section"},{"location":"internals/","page":"Crispulator Internals","title":"Crispulator Internals","text":"Crispulator.Library\nCrispulator.Barcode\nCrispulator.KDPhenotypeRelationship\nCrispulator.Cas9Behavior\nCrispulator.CRISPRn\nCrispulator.CRISPRi\nCrispulator.Delta","category":"page"},{"location":"internals/#Crispulator.Library","page":"Crispulator Internals","title":"Crispulator.Library","text":"Wrapper containing all library construction parameters\n\nknockdown_dist\nDistribution of guide knockdown efficiencies\nknockdown_probs\nmax_phenotype_dists\nMaximum phenotype categories mapped to their probability of being selected and the distribution to draw from if they are selected.\nExample\nThe basic layout is a Dict mapping a class name to a tuple of the probability of selecting this class and then the Distributions.Sampleable from which to draw a random phenotype from this class. The probabilities across all the classes should add up to 1.\nmax_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(\n    :inactive => (0.60, Delta(0.0)),\n    :negcontrol => (0.1, Delta(0.0)),\n    :increasing => (0.3, truncated(Normal(0.1, 0.1), 0.025, 1)),\n);\nLibrary(max_phenotype_dists, CRISPRi());\nFor example, here we are making three different classes of \"genes\": the first group are :inactive, i.e. they have no phenotype, so we'll set their phenotypes to 0.0 using a Crispulator.Delta. We'll also make them 60% of all the genes. The second group are the negative controls :negcontrol (the only required group) which make up 10% of the population of genes and also have no effect. The final group is :increasing which makes up 30% of all genes and which are represented by a Normal(μ=0.1, σ=0.1) distribution clamped between 0.025 and 1.\n\nphenotype_probs\nkd_phenotype_relationships\nKnockdown-phenotype relationships mapped to their probability of being selected and their respective Crispulator.KDPhenotypeRelationship\n\nrelationship_probs\ncas9_behavior\nWhether this library is Crispulator.CRISPRi or Crispulator.CRISPRn\n\n\n\n\n\n","category":"type"},{"location":"internals/#Crispulator.Barcode","page":"Crispulator Internals","title":"Crispulator.Barcode","text":"Any entity that is tracked through the pooled experiment. For CRISPR screens, this is equivalent to the sgRNA. This object stores properties relating to the performance of this entity in the screen.\n\ngene\nThe target gene id\nknockdown\nThe knockdown efficiency\ntheo_phenotype\nThe theoretical phenotype exerted\nobs_phenotype\nThe observed phenotype (only relevant for FACS screens)\nbehavior\nThe gene behavior, either sigmoidal or linear\nclass\nThe gene activity class, either inactive, up or down\ninitial_freq\nInitial frequency after transfection\n\n\n\n\n\n","category":"type"},{"location":"internals/#Crispulator.KDPhenotypeRelationship","page":"Crispulator Internals","title":"Crispulator.KDPhenotypeRelationship","text":"A type representing a relationship between degree of knockdown and effect on phenotype\n\n\n\n\n\n","category":"type"},{"location":"internals/#Crispulator.Cas9Behavior","page":"Crispulator Internals","title":"Crispulator.Cas9Behavior","text":"A type representing the behavior of different Cas9s\n\n\n\n\n\n","category":"type"},{"location":"internals/#Crispulator.CRISPRn","page":"Crispulator Internals","title":"Crispulator.CRISPRn","text":"struct CRISPRn <: Crispulator.Cas9Behavior\n\nCRISPR KO behavior is more complex since sgRNA-directed DNA damage repair is stochastic. We assume that 2/3 of repair events at a given locus lead to a frameshift, and that the screen is carried out in diploid cells. The assumption that only bi-allelic frame-shift mutations lead to a phenotype in CRISPRn screens for most sgRNAs is supported by the empirical finding that in-frame deletions mostly do not show strong phenotypes, unless they occur in regions encoding conserved residues or domains[2]\n\n[2]: Horlbeck MA, Gilbert LA, Villalta JE, Adamson B, Pak RA, Chen Y, Fields AP, Park CY, Corn JE, Kampmann M, Weissman JS: Compact and highly active next- generation libraries for CRISPR-mediated gene repression and activation. Elife 2016, 5.\n\n\n\n\n\n","category":"type"},{"location":"internals/#Crispulator.CRISPRi","page":"Crispulator Internals","title":"Crispulator.CRISPRi","text":"struct CRISPRi <: Crispulator.Cas9Behavior\n\nCRISPRi behavior is simply determined by the activity of the guide\n\n\n\n\n\n","category":"type"},{"location":"internals/#Crispulator.Delta","page":"Crispulator Internals","title":"Crispulator.Delta","text":"struct Delta <: Distributions.Sampleable{Distributions.Univariate, Distributions.Discrete}\n\nConstructs a delta function at a given δ value. This distribution always emits the same value.\n\n\n\n\n\n","category":"type"},{"location":"internals/#Miscellaneous-functions","page":"Crispulator Internals","title":"Miscellaneous functions","text":"","category":"section"},{"location":"internals/","page":"Crispulator Internals","title":"Crispulator Internals","text":"Crispulator.signal\nCrispulator.noise\nCrispulator.linear\nCrispulator.sigmoid\nCrispulator.auprc\nCrispulator.auroc\nCrispulator.venn","category":"page"},{"location":"internals/#Crispulator.signal","page":"Crispulator Internals","title":"Crispulator.signal","text":"signal(guide_data; log2fc_col)\n\n\nComputes the signal in an experiment, where the experimental signal is defined to be the average of signal of true hit genes. That value, the true hit gene signal, is the average ratio of the log2 fold change for a guide targeting a specific gene to the guide's theoretical phenotype.\n\nfrac1N_true times k sum_i=1^N_true sum_j=1^k\nfraclog_2 fc_ijtexttheo phenotype_ij\n\nwhere N_true is the number of true hit genes and k is the number of genes.\n\n\n\n\n\n","category":"function"},{"location":"internals/#Crispulator.noise","page":"Crispulator Internals","title":"Crispulator.noise","text":"noise(guide_data; log2fc_col)\n\n\nComputes the noise in an experiment, where noise is defined to be the standard deviation of log2 fold change in the negative controls. If more than 2 bins are used then the name of log2 fold change information can be provided to differentiate between multiple log 2 fold changes.\n\n\n\n\n\n","category":"function"},{"location":"internals/#Crispulator.linear","page":"Crispulator Internals","title":"Crispulator.linear","text":"linear(x, l)\n\n\nGiven value(s) x apply a simple linear function with maximum value l\n\n\n\n\n\n","category":"function"},{"location":"internals/#Crispulator.sigmoid","page":"Crispulator Internals","title":"Crispulator.sigmoid","text":"sigmoid(x, l, k, p)\n\n\nGiven value(s) x apply the sigmoidal function with maximum value l, a steepness k, and an inflection point p.\n\n\n\n\n\n","category":"function"},{"location":"internals/#Crispulator.auprc","page":"Crispulator Internals","title":"Crispulator.auprc","text":"auprc(scores, classes, pos_labels; rev)\n\n\nComputes the area under the Precision-Recall curve using a lower trapezoidal estimator, which is more accurate for skewed datasets.\n\nK. Boyd, K. H. Eng, and C. D. Page, “Area under the Precision-Recall Curve: Point Estimates and Confidence Intervals,” in Machine Learning and Knowledge Discovery in Databases, H. Blockeel, K. Kersting, S. Nijssen, and F. Železný, Eds. Springer Berlin Heidelberg, 2013, pp. 451–466.\n\n\n\n\n\n","category":"function"},{"location":"internals/#Crispulator.auroc","page":"Crispulator Internals","title":"Crispulator.auroc","text":"auroc(scores, classes, pos_labels; rev)\n\n\nOptimized function for computing the area under the receiver operator characteristic curve.\n\n\n\n\n\n","category":"function"},{"location":"internals/#Crispulator.venn","page":"Crispulator Internals","title":"Crispulator.venn","text":"venn(scores, classes, pos_labels; rev)\n\n\nGiven N positive examples, computes the percentage of the top N/2 hits that are correct\n\n\n\n\n\n","category":"function"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The easiest way to use Crispulator is via the command line config file. If this format is too constraining, the Custom Simulations has a detailed walk- through of writing a custom simulation where each step can be modified according to need.","category":"page"},{"location":"tutorial/#Graphical-overview","page":"Tutorial","title":"Graphical overview","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The simulation is laid out in the following manner:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: )","category":"page"},{"location":"tutorial/#Getting-started","page":"Tutorial","title":"Getting started","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"First, navigate to the Crispulator directory.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"tip: Tip\nYou can find the directory by running$ julia -e 'println(Pkg.dir(\"Crispulator\"))'","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"There should be a YAML file called example_config.yml. Open this is in a text editor and it should look like this","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"# This is an example configuration file. Whitespace is important.\n\n# Settings pertaining to the library design\nlibrary:\n    genome:\n        num-genes: 500\n        num-guides-per-gene: 5\n        frac-increasing-genes: 0.02 # fraction of genes with a positive phenotype\n        frac-decreasing-genes: 0.1 # fraction of genes with a negative phenotype\n\n    guides:\n        crispr-type: CRISPRn # either CRISPRi or CRISPRn\n        frac-high-quality: 0.9 # fraction of high quality guides\n        mean-high-quality-kd: 0.85 # mean knockdown by a high quality guide (CRISPRi only)\n\nscreen:\n    type: facs # either facs or growth\n    num-runs: 10 # how many independent runs\n\n    representation: # integer value, how much larger are samples than the library\n        - transfection: 100\n        - selection: 100\n        - sequencing: 100\n\n# screen-type specific parameters\n\n    bin-size: 0.25 # size of tail to sample from, must be between 0 and 0.5 (FACS only)\n    std-noise: 1 # (FACS only)\n    num-bottlenecks: 10 # (Growth only)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"This gives access to most dials in the simulation, if something is missing than see Custom Simulations.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Now, lets remove all genes that have a positive phenotype by changing line 8 to 0.0:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"        frac-increasing-genes: 0.0 # fraction of genes with a positive phenotype","category":"page"},{"location":"tutorial/#Running-simulation","page":"Tutorial","title":"Running simulation","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Now, we can actually run the code by executing the following command","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"julia run.jl config example_config.yml test_output","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"tip: Tip\nHere config tells CRISPulator to use the provided config example_config.yml and test_output is the directory where the results will be saved. This directory will be created if it doesn't exist.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The output should look like","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"readdir(\".\")","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"readdir(\"..\")","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"readdir(joinpath(\"..\"))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"include(joinpath(\"..\", \"..\", \"commands.jl\")) #hide\nbootstrap_config(\"example_config.yml\", \"test_output\", false) #hide","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The test_output/ directory should now be populated with all the files","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"foreach(println, readdir(\"test_output\")) # hide","category":"page"},{"location":"tutorial/#Output","page":"Tutorial","title":"Output","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The folder contains one of the raw count scatterplots (left) and a volcano plot of mean log2 fold change versus significance of each gene (right)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: ) (Image: )","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"It also has a useful table that contains all the summary statistic information.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using DataFrames # hide\nusing CSV # hide\nfirst(CSV.read(joinpath(\"test_output\", \"results_table.csv\"), DataFrame), 10) # hide","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The table below describes each column","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Column Name Meaning\nmethod Which summary statistic was used (e.g. Crispulator.auprc)\nmeasure Whether the score is only for increasing genes (inc), decreasing (dec) or both (incdec). Allows independent evaluation on which type of genes the screen can accurately evaluate.\ngenetype Whether the score is for linear, sigmoidal, or all genes (see Crispulator.KDPhenotypeRelationship). Helps determine if CRISPRn or CRISPRi is better for this design.\nmean_score Average score\nstd_score Standard deviation in scores\nconf_max Upper limit of 95% confidence interval\nconf_min Lower limit of 95% confidence interval\nn Number of independent replicates","category":"page"}]
}
