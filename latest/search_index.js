var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#CRISPulator-1",
    "page": "Home",
    "title": "CRISPulator",
    "category": "section",
    "text": "A pooled genetic screen simulatorPooled screens are very useful, but it is difficult to select all of the necessary parameters a priori. This work aims to explore the importance of various screen parameters and the relevance of design choices on the downstream analysis and interpretation."
},

{
    "location": "index.html#Setup-1",
    "page": "Home",
    "title": "Setup",
    "category": "section",
    "text": "To run the simulation, you will need a recent version of Julia installed and in your PATH. Then navigate into the root directory of the project and run julia. Run the following command:$ julia -e \'Pkg.update(); Pkg.add(\"Crispulator\")\'this copies Crispulator over to the Julia package directory and installs all of its dependencies. Note the $ just represents the prompt, you don\'t need to type it. Make sure to run the tests and verify that everything is passing$ julia -e \'Pkg.test(\"Crispulator\")\'"
},

{
    "location": "index.html#Quickstart-1",
    "page": "Home",
    "title": "Quickstart",
    "category": "section",
    "text": "For most simple cases, no writing of Julia code is necessary. See the Tutorial for using Crispulator in this manner."
},

{
    "location": "index.html#Advanced-1",
    "page": "Home",
    "title": "Advanced",
    "category": "section",
    "text": "See the Custom Simulations for a step-by-step guide to writing a custom simulation on top of Crispulator. Additionally, the Simulation Internals page has in-depth documentation needed for more advanced usage"
},

{
    "location": "tutorial.html#",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial.html#Tutorial-1",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "section",
    "text": "The easiest way to use Crispulator is via the command line config file. If this format is too constraining, the Custom Simulations has a detailed walk- through of writing a custom simulation where each step can be modified according to need."
},

{
    "location": "tutorial.html#Graphical-overview-1",
    "page": "Tutorial",
    "title": "Graphical overview",
    "category": "section",
    "text": "The simulation is laid out in the following manner:(Image: )"
},

{
    "location": "tutorial.html#Getting-started-1",
    "page": "Tutorial",
    "title": "Getting started",
    "category": "section",
    "text": "First, navigate to the Crispulator directory.tip: Tip\nYou can find the directory by running$ julia -e \'println(Pkg.dir(\"Crispulator\"))\'There should be a YAML file called example_config.yml. Open this is in a text editor and it should look like this# This is an example configuration file. Whitespace is important.\n\n# Settings pertaining to the library design\nlibrary:\n    genome:\n        num-genes: 500\n        num-guides-per-gene: 5\n        frac-increasing-genes: 0.02 # fraction of genes with a positive phenotype\n        frac-decreasing-genes: 0.1 # fraction of genes with a negative phenotype\n\n    guides:\n        crispr-type: CRISPRn # either CRISPRi or CRISPRn\n        frac-high-quality: 0.9 # fraction of high quality guides\n        mean-high-quality-kd: 0.85 # mean knockdown by a high quality guide (CRISPRi only)\n\nscreen:\n    type: facs # either facs or growth\n    num-runs: 10 # how many independent runs\n\n    representation: # integer value, how much larger are samples than the library\n        - transfection: 100\n        - selection: 100\n        - sequencing: 100\n\n# screen-type specific parameters\n\n    bin-size: 0.25 # size of tail to sample from, must be between 0 and 0.5 (FACS only)\n    std-noise: 1 # (FACS only)\n    num-bottlenecks: 10 # (Growth only)This gives access to most dials in the simulation, if something is missing than see Custom Simulations.Now, lets remove all genes that have a positive phenotype by changing line 8 to 0.0:        frac-increasing-genes: 0.0 # fraction of genes with a positive phenotype"
},

{
    "location": "tutorial.html#Running-simulation-1",
    "page": "Tutorial",
    "title": "Running simulation",
    "category": "section",
    "text": "Now, we can actually run the code by executing the following commandjulia src/run.jl config example_config.yml test_outputtip: Tip\nHere config tells CRISPulator to use the provided config example_config.yml and test_output is the directory where the results will be saved. This directory will be created if it doesn\'t exist.The output should look likeoutput = readstring(`julia ../../src/run.jl config example_config.yml test_output`) # hide\nprintln(output) # hideThe test_output/ directory should now be populated with all the filesforeach(println, readdir(\"test_output\")) # hide"
},

{
    "location": "tutorial.html#Output-1",
    "page": "Tutorial",
    "title": "Output",
    "category": "section",
    "text": "The folder contains one of the raw count scatterplots (left) and a volcano plot of mean log2 fold change versus significance of each gene (right)(Image: ) (Image: )It also has a useful table that contains all the summary statistic information.using DataFrames # hide\nhead(readtable(joinpath(\"test_output\", \"results_table.csv\"))) # hideThe table below describes each columnColumn Name Meaning\nmethod Which summary statistic was used (e.g. Simulation.auprc)\nmeasure Whether the score is only for increasing genes (inc), decreasing (dec) or both (incdec). Allows independent evaluation on which type of genes the screen can accurately evaluate.\ngenetype Whether the score is for linear, sigmoidal, or all genes (see Simulation.KDPhenotypeRelationship). Helps determine if CRISPRn or CRISPRi is better for this design.\nmean_score Average score\nstd_score Standard deviation in scores\nconf_max Upper limit of 95% confidence interval\nconf_min Lower limit of 95% confidence interval\nn Number of independent replicates"
},

{
    "location": "custom.html#",
    "page": "Custom Simulations",
    "title": "Custom Simulations",
    "category": "page",
    "text": ""
},

{
    "location": "custom.html#Custom-Simulations-1",
    "page": "Custom Simulations",
    "title": "Custom Simulations",
    "category": "section",
    "text": "This is an example of building a custom simulation of a FACS-based screen.tip: Tip\nThis section is complementary to the Implementation section of the paper"
},

{
    "location": "custom.html#Setup-1",
    "page": "Custom Simulations",
    "title": "Setup",
    "category": "section",
    "text": "using Simulation\nusing Gadfly\nusing DataFrames\nusing DistributionsLets first start the Julia REPL or a Julia session inside of a Jupyter Notebook and load the packages we\'ll need:using Gadfly\ninclude(joinpath(Pkg.dir(\"Crispulator\"), \"src\", \"simulation\", \"load.jl\"))"
},

{
    "location": "custom.html#Basic-screen-parameters-1",
    "page": "Custom Simulations",
    "title": "Basic screen parameters",
    "category": "section",
    "text": "First lets design a simple Simulation.FacsScreen with 250 genes with 5 guides per gene. Lets say that we make sure we have 1000x as many cells as guides during transfection (representation) and sorting (bottleneck_representation) and 1000x as many reads as guides during sequencing (seq_depth). We\'ll leave the rest of the values as the defaults and print the objects = FacsScreen()\ns.num_genes = 250\ns.coverage = 5\ns.representation = 1000\ns.bottleneck_representation = 1000\ns.seq_depth = 1000\nprintln(s)"
},

{
    "location": "custom.html#Construction-of-true-phenotype-distribution-1",
    "page": "Custom Simulations",
    "title": "Construction of true phenotype distribution",
    "category": "section",
    "text": "Next, lets make our distribution of true phenotypes. The basic layout is a Dict mapping a class name to a tuple of the probability of selecting this class and then the Distributions.Sampleable from which to draw a random phenotype from this class. The probabilities across all the classes should add up to 1.For example, here we are making three different classes of \"genes\": the first group are :inactive, i.e. they have no phenotype, so we\'ll set their phenotypes to 0.0 using a Simulation.Delta. We\'ll also make them 60% of all the genes. The second group are the negative controls :negcontrol (the only required group) which make up 10% of the population of genes and also have no effect. The final group is :increasing which makes up 30% of all genes and which are represented by a Normal(μ=0.1, σ=0.1) distribution clamped between 0.025 and 1.max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(\n    :inactive => (0.60, Delta(0.0)),\n    :negcontrol => (0.1, Delta(0.0)),\n    :increasing => (0.3, TruncatedNormal(0.1, 0.1, 0.025, 1)),\n);note: Note\nThe :negcontrol class needs to be present because Crispulator normalizes the frequencies of all other guides against the median frequency of the negative control guides. Also the distribution of :negcontrol guides serve as the null distribution against which the log2 fold changes of guides targeting a specific gene are assayed to calculate a statistical significance of the shift for each gene. See Simulation.differences_between_bins for more details."
},

{
    "location": "custom.html#Library-construction-1",
    "page": "Custom Simulations",
    "title": "Library construction",
    "category": "section",
    "text": "Now, we actually build the library. Here we\'re making a Simulation.CRISPRi library and then getting the guides that were built from the true phenotype distribution that we constructed above and we also get the frequency of each guide in the library.lib = Library(max_phenotype_dists, CRISPRi())\nguides, guide_freqs_dist = construct_library(s, lib);Lets first look at what the true phenotype distribution of our different classes of guides looks likedf = DataFrame(Dict(\n    :phenotype=>map(x->x.theo_phenotype, guides),\n    :class=>map(x->x.class, guides),\n    :freq=>pdf.(guide_freqs_dist, 1:(length(guides)))\n))\nplot(df, x=:phenotype, color=:class, Geom.histogram, Guide.ylabel(\"Number of guides\"),\nGuide.title(\"Guide phenotype distribution\"))As you can see, most guides should have a phenotype of 0. In FACS Screens this is equivalent to having no preference to being in either the left (bin1) or right (bin2) bins. The :increasing genes have a small preference to be in the right bin.We can also look at the frequency of each guide in the library, which follows a Log-Normal distribution.plot(df, x=:freq, color=:class, Geom.histogram(position=:stack),\n    Guide.xlabel(\"Frequency\"), Guide.ylabel(\"Number of guides\"),\n    Guide.title(\"Frequencies of guides in simulated library\"))"
},

{
    "location": "custom.html#Performing-the-screen-1",
    "page": "Custom Simulations",
    "title": "Performing the screen",
    "category": "section",
    "text": "Now, we\'ll actually perform the screen. We\'ll first perform the transection via Simulation.transfect, followed by the selection process via Simulation.select:cells, cell_phenotypes = transfect(s, lib, guides, guide_freqs_dist)\nbin_cells = select(s, cells, cell_phenotypes, guides)\nfreqs = counts_to_freqs(bin_cells, length(guides));Lets look at what the observed phenotype distribution looks like when the selection was performed:df = DataFrame(Dict(\n    :phenotype=>map(x->x.theo_phenotype, guides),\n    :class=>map(x->x.class, guides),\n    :obs_freq=>map(x->x.obs_phenotype, guides)\n))\nplot(df, x=:obs_freq, Geom.density, Guide.xlabel(\"Observed phenotype on FACS machine\"),\nGuide.title(\"Kernel density estimate of guide observed phenotypes\"), Guide.ylabel(\"ρ\"))As you can see, this looks like many FACS plots, e.g. when looking at density along the fluorescence channel. A quick sanity check is that we should see a slight enrichment of the frequency of :increasing genes on the right sideplot(df, x=:obs_freq, color=:class, Geom.density, Guide.xlabel(\"Observed phenotype on FACS machine\"),\nGuide.title(\"Kernel density estimate of guide observed phenotypes\"), Guide.ylabel(\"ρ\"))And that is what we see. The change is really small (this is pretty usual), but the later analysis will be able to pull out the increasing genes."
},

{
    "location": "custom.html#Sequencing-and-Analysis-1",
    "page": "Custom Simulations",
    "title": "Sequencing and Analysis",
    "category": "section",
    "text": "Now we\'ll use Simulation.sequencing to simulate sequencing by transforming the guide frequencies into a Categorical distribution and drawing a random sample of reads from this distribution. Finally, we\'ll use the Simulation.differences_between_bins function to compute the differences between bins on a per-guide level (guide_data) and per-gene level (gene_data).raw_data = sequencing(Dict(:bin1=>s.seq_depth, :bin2=>s.seq_depth), guides, freqs)\nguide_data, gene_data = differences_between_bins(raw_data);Here\'s what the per-guide data looks like:head(guide_data) # hidetip: Tip\nSee Simulation.differences_between_bins for details on what each column means.And the gene level datahead(gene_data) # hideWe can generate standard pooled screen plots from this dataset. Like a count scatterplot:nopseudo = guide_data[(guide_data[:counts_bin1] .> 0.5) .& (guide_data[:counts_bin2] .> 0.5), :]\nplot(nopseudo, x=:counts_bin1, y=:counts_bin2, color=:class, Scale.x_log10,\nScale.y_log10, Theme(highlight_width=0pt), Coord.cartesian(fixed=true),\nGuide.xlabel(\"log counts bin1\"), Guide.ylabel(\"log counts bin2\"))And a volcano plot:plot(gene_data, x=:mean_bin2_div_bin1, y=:pvalue_bin2_div_bin1, color=:class, Theme(highlight_width=0pt),\nGuide.xlabel(\"mean log2 fold change\"), Guide.ylabel(\"-log10 pvalue\"))And finally we can see how well we can differentiate between the different classes using Area Under the Precision-Recall Curve (Simulation.auprc)auprc(gene_data[:pvalmeanprod_bin2_div_bin1], gene_data[:class], Set([:increasing]))[1]Simulation.auroc and Simulation.venn are also good summary statistics."
},

{
    "location": "internals.html#",
    "page": "Simulation Internals",
    "title": "Simulation Internals",
    "category": "page",
    "text": ""
},

{
    "location": "internals.html#Simulation-Internals-1",
    "page": "Simulation Internals",
    "title": "Simulation Internals",
    "category": "section",
    "text": "More in-depth documentation on specific types and functions."
},

{
    "location": "internals.html#Contents-1",
    "page": "Simulation Internals",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"internals.md\"]"
},

{
    "location": "internals.html#Index-1",
    "page": "Simulation Internals",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"internals.md\"]"
},

{
    "location": "internals.html#Simulation.FacsScreen",
    "page": "Simulation Internals",
    "title": "Simulation.FacsScreen",
    "category": "type",
    "text": "A type representing the parameters used in a typical FACS screen.\n\nnum_genes\nNumber of genes targeted by the screen\ncoverage\nNumber of guides per gene\nrepresentation\nNumber of cells with each guide. representation=10 means that there are 10 times as many cells as guides during transfection. The actual number of cells per guide post-transfection will be less depending on the MOI\nmoi\nThe multiplicity of infection, lambda, of the screen. We model this as a Poisson process during transfection (see Simulation.transfect).\nnote: Note\nWe do not model multiple infections. We assume that the MOI is properly selected and less than half the cells are transfected by any virus, i.e. lambda lt 05 and then select only the cells that have a single transfection occurrence:P(x=1Poisson())For lambda = 025 this works out to being approx 195 of the number of cells (num_genes * coverage * representation).\n\nσ\nThe standard deviation expected for cells during FACS sorting. This should be set according to the biological variance experimentally observed, e.g. in the fluorescence intensity of isogenic cells\nbin_info\nRange of guide phenotypes to collect in each bin\nExample\nIn the following example\np = 0.05\nbin_info = Dict(:bin1 => (0.0, p), :bin2 => (1.0-p, 1.0))\nThe 5th percentile of cells sorted according to their phenotype (fluorescence, size, etc) will be compared to the 95th percentile.\n\nbottleneck_representation\nNumber of cells sorted expressed as an integer multiple of the number of guides\nseq_depth\nSequencing depth as a integer multiple of the number of guides, i.e. seq_depth=10 is equivalent to 10 * num_genes * coverage reads.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.GrowthScreen",
    "page": "Simulation Internals",
    "title": "Simulation.GrowthScreen",
    "category": "type",
    "text": "A type representing the parameters used in a typical growth-based screen.\n\nnum_genes\nNumber of genes targeted by the screen\ncoverage\nNumber of guides per gene\nrepresentation\nNumber of cells with each guide. representation=10 means that there are 10 times as many cells as guides during transfection. The actual number of cells per guide post-transfection will be less depending on the MOI\nmoi\nThe multiplicity of infection, lambda, of the screen. We model this as a Poisson process during transfection (see Simulation.transfect).\nnote: Note\nWe do not model multiple infections. We assume that the MOI is properly selected and less than half the cells are transfected by any virus, i.e. lambda lt 05 and then select only the cells that have a single transfection occurrence:P(x=1Poisson())For lambda = 025 this works out to being approx 195 of the number of cells (num_genes * coverage * representation).\n\nseq_depth\nSequencing depth as a integer multiple of the number of guides, i.e. seq_depth=10 is equivalent to 10 * num_genes * coverage reads.\n\nbottleneck_representation\nFor growth screens, how much of a bottleneck is applied. This is minimum number of cells that is kept when passaging the pool of cells. This is expressed as an integer multiple of the number of guides.\nnum_bottlenecks\nFor growth screens, how many bottlenecks are applied. This is the integer number of passages during the growth screen.\nnoise\nBefore each passage, the theoretical phenotype of each cell is convolved with a normal noise distribution with a standard deviation, σ, dictated by this parameter. This should be set based on an expectation of noisiness of the subsampling\n\n\n\n"
},

{
    "location": "internals.html#Simulation-Types-1",
    "page": "Simulation Internals",
    "title": "Simulation Types",
    "category": "section",
    "text": "Simulation.FacsScreen\nSimulation.GrowthScreen"
},

{
    "location": "internals.html#Simulation.construct_library",
    "page": "Simulation Internals",
    "title": "Simulation.construct_library",
    "category": "function",
    "text": "construct_library(setup, lib)\n\n\nConstructs the guide library for N genes with coverage number of guides per gene. Returns a tuple of guides and their relative frequencies (assigned randomly).\n\n\n\n"
},

{
    "location": "internals.html#Simulation.transfect",
    "page": "Simulation Internals",
    "title": "Simulation.transfect",
    "category": "function",
    "text": "transfect(setup, lib, guides, guide_freqs_dist)\n\n\nSimulates the transfection the of the library of guides (guides) according to the parameters in Simulation.Library lib and the frequencies represented by the Categorical distribution guide_freqs_dist.\n\nThis function creates a population of cells post-transfection where each cell has single guide RNA present according to the frequencies of the guides in the present in the library. The cells are then grown up so that there are setup.bottleneck_representation * num_guides of them. For FACS screens we assume a simple linear expansion of the cells, while for Growth screens the cells grow according as a function of their phenotype and degree of knockdown.\n\n\n\n"
},

{
    "location": "internals.html#Base.Sort.select",
    "page": "Simulation Internals",
    "title": "Base.Sort.select",
    "category": "function",
    "text": "select(setup, cells, cell_phenotypes, guides; debug)\n\n\nGiven cells, a vector of integers, and guides, a vector of barcodes, performs simulated FACS sorting of the cells into bins with the given cutoffs. σ refers to the spread of the phenotype during the FACS screen.\n\n\n\nselect(setup, initial_cells, initial_cell_phenotypes, guides; debug)\n\n\nGrowth Screen selection\n\n\n\n"
},

{
    "location": "internals.html#Simulation.sequencing",
    "page": "Simulation Internals",
    "title": "Simulation.sequencing",
    "category": "function",
    "text": "sequencing(depths, guides, samples)\n\n\nSimulates next-gen sequencing by generating a Categorical distribution based on frequencies of each guide in each bin and then randomly samples from this distributions up to the depth provided in depths for each bin. Returns a dict of DataFrames with the bin names as the keys. This object can then be passed through Simulation.counts_to_freqs followed by Simulation.differences_between_bins.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.counts_to_freqs",
    "page": "Simulation Internals",
    "title": "Simulation.counts_to_freqs",
    "category": "function",
    "text": "counts_to_freqs(bin_cells, guide_count)\n\n\nConverts raw counts to frequencies by dividing by the total number of reads for each sample.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.differences_between_bins",
    "page": "Simulation Internals",
    "title": "Simulation.differences_between_bins",
    "category": "function",
    "text": "differences_between_bins(raw_data)\n\n\nGiven the raw data from Simulation.sequencing returns two DataFrames\n\nguide_data: This DataFrame contains the per-guide level data including the  log2 fold change in the normalized frequencies of each guide between each  pairwise combination of bins. Thus, if there are n bins, then it computes  the log2 fold changes for the fracn2(n-2) combinations\ngene_data: This DataFrame contains the same information but grouped by  gene. The log2 fold change data from the first DataFrame is used to calculate  the average log2 fold change per gene and a pvalue computed using a  Mann-Whitney U-test as  measure of how consistently shifted the guides are of this gene versus the  population of negative control guides. (see below for more info)\n\nA typical 2 bin guide_data DataFrame contains the following columns:\n\nColumn Name Meaning\ngene the gene ID of that this guide targets\nknockdown activity of the guide on 0 to 1 scale, where 1 is complete knockout\nbarcodeid the ID of this specific guide\ntheo_phenotype expected phenotype of this guide, generally a -1 to 1 scale\nbehavior whether the target gene displays a linear or sigmoidal response to incomplete knockdown (see Simulation.Library for more details)\nclass which phenotype distribution the target gene was drawn from (see Simulation.Library for more details). Serves as the \"ground truth\" label against which screen performance is evaluated, e.g. with Simulation.auprc\ninitial_freq frequency of guide post-transfection (see Simulation.transfect)\ncounts_bin1 the raw number of reads for each guide in the first bin\nfreqs_bin1 the number of reads for each guide divided by the total number of reads in this bin\nrel_freqs_bin1 the frequency of each guide divided by the median frequency of negative control guides\ncounts_bin2 the raw number of reads for each guide in the second bin\nfreqs_bin2 the number of reads for each guide divided by the total number of reads in this bin\nrel_freqs_bin2 the frequency of each guide divided by the median frequency of negative control guides for this bin\nlog2fc_bin2_div_bin1 the log2 fold change in relative guide frequencies between bin2 and bin1, equivalent to log2(rel_freqs_bin2/rel_freqs_bin1)\n\nA typical gene_data DataFrame contains the following data:\n\nColumn Name Meaning\ngene this gene\'s ID\nclass see above\nbehavior see above\nmean_bin2_div_bin1 the mean log 2 fold change in relative frequencies between from bin1 to bin2 for all the guides targeting this gene. Calculated as frac1ksum_k textlog2fc_bin2_div_bin1_k for the k guides targeting each gene\npvalue_bin2_div_bin1 the -log10 pvalue of the log2 fold changes of all guides targeting this gene as computed by the non-parametric Mann-Whitney U-test. A measure of the consistency of the log 2 fold changes[1]\nabsmean_bin2_div_bin1 absolute value of mean_bin2_div_bin1 per-gene\npvalmeanprod_bin2_div_bin1 mean_bin2_div_bin1 multiplied with the pvalue_bin2_div_bin1 per-gene\n\nFurther reading\n\n[1]: Kampmann M, Bassik MC, Weissman JS. Integrated platform for genome-wide screening and construction of high-density genetic interaction maps in mammalian cells. Proc Natl Acad Sci U S A. 2013;110:E2317–26.\n\n\n\n"
},

{
    "location": "internals.html#Key-functions-1",
    "page": "Simulation Internals",
    "title": "Key functions",
    "category": "section",
    "text": "Simulation.construct_library\nSimulation.transfect\nSimulation.select\nSimulation.sequencing\nSimulation.counts_to_freqs\nSimulation.differences_between_bins"
},

{
    "location": "internals.html#Simulation.Library",
    "page": "Simulation Internals",
    "title": "Simulation.Library",
    "category": "type",
    "text": "Wrapper containing all library construction parameters\n\nknockdown_dist\nDistribution of guide knockdown efficiencies\nknockdown_probs\nmax_phenotype_dists\nMaximum phenotype categories mapped to their probability of being selected and the distribution to draw from if they are selected.\nExample\nThe basic layout is a Dict mapping a class name to a tuple of the probability of selecting this class and then the Distributions.Sampleable from which to draw a random phenotype from this class. The probabilities across all the classes should add up to 1.\nmax_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(\n    :inactive => (0.60, Delta(0.0)),\n    :negcontrol => (0.1, Delta(0.0)),\n    :increasing => (0.3, TruncatedNormal(0.1, 0.1, 0.025, 1)),\n);\nLibrary(max_phenotype_dists, CRISPRi());\nFor example, here we are making three different classes of \"genes\": the first group are :inactive, i.e. they have no phenotype, so we\'ll set their phenotypes to 0.0 using a Simulation.Delta. We\'ll also make them 60% of all the genes. The second group are the negative controls :negcontrol (the only required group) which make up 10% of the population of genes and also have no effect. The final group is :increasing which makes up 30% of all genes and which are represented by a Normal(μ=0.1, σ=0.1) distribution clamped between 0.025 and 1.\n\nphenotype_probs\nkd_phenotype_relationships\nKnockdown-phenotype relationships mapped to their probability of being selected and their respective Simulation.KDPhenotypeRelationship\n\nrelationship_probs\ncas9_behavior\nWhether this library is Simulation.CRISPRi or Simulation.CRISPRn\n\n\n\n"
},

{
    "location": "internals.html#Simulation.Barcode",
    "page": "Simulation Internals",
    "title": "Simulation.Barcode",
    "category": "type",
    "text": "Any entity that is tracked through the pooled experiment. For CRISPR screens, this is equivalent to the sgRNA. This object stores properties relating to the performance of this entity in the screen.\n\ngene\nThe target gene id\nknockdown\nThe knockdown efficiency\ntheo_phenotype\nThe theoretical phenotype exerted\nobs_phenotype\nThe observed phenotype (only relevant for FACS screens)\nbehavior\nThe gene behavior, either sigmoidal or linear\nclass\nThe gene activity class, either inactive, up or down\ninitial_freq\nInitial frequency after transfection\n\n\n\n"
},

{
    "location": "internals.html#Simulation.KDPhenotypeRelationship",
    "page": "Simulation Internals",
    "title": "Simulation.KDPhenotypeRelationship",
    "category": "type",
    "text": "A type representing a relationship between degree of knockdown and effect on phenotype\n\n\n\n"
},

{
    "location": "internals.html#Simulation.Cas9Behavior",
    "page": "Simulation Internals",
    "title": "Simulation.Cas9Behavior",
    "category": "type",
    "text": "A type representing the behavior of different Cas9s\n\n\n\n"
},

{
    "location": "internals.html#Simulation.CRISPRn",
    "page": "Simulation Internals",
    "title": "Simulation.CRISPRn",
    "category": "type",
    "text": "type CRISPRn <: Simulation.Cas9Behavior\n\nCRISPR KO behavior is more complex since sgRNA-directed DNA damage repair is stochastic. We assume that 2/3 of repair events at a given locus lead to a frameshift, and that the screen is carried out in diploid cells. The assumption that only bi-allelic frame-shift mutations lead to a phenotype in CRISPRn screens for most sgRNAs is supported by the empirical finding that in-frame deletions mostly do not show strong phenotypes, unless they occur in regions encoding conserved residues or domains[2]\n\n[2]: Horlbeck MA, Gilbert LA, Villalta JE, Adamson B, Pak RA, Chen Y, Fields AP, Park CY, Corn JE, Kampmann M, Weissman JS: Compact and highly active next- generation libraries for CRISPR-mediated gene repression and activation. Elife 2016, 5.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.CRISPRi",
    "page": "Simulation Internals",
    "title": "Simulation.CRISPRi",
    "category": "type",
    "text": "type CRISPRi <: Simulation.Cas9Behavior\n\nCRISPRi behavior is simply determined by the activity of the guide\n\n\n\n"
},

{
    "location": "internals.html#Simulation.Delta",
    "page": "Simulation Internals",
    "title": "Simulation.Delta",
    "category": "type",
    "text": "type Delta <: Distributions.Sampleable{Distributions.Univariate,Distributions.Discrete}\n\nConstructs a delta function at a given δ value. This distribution always emits the same value.\n\n\n\n"
},

{
    "location": "internals.html#Miscellaneous-Types-1",
    "page": "Simulation Internals",
    "title": "Miscellaneous Types",
    "category": "section",
    "text": "Simulation.Library\nSimulation.Barcode\nSimulation.KDPhenotypeRelationship\nSimulation.Cas9Behavior\nSimulation.CRISPRn\nSimulation.CRISPRi\nSimulation.Delta"
},

{
    "location": "internals.html#Simulation.signal",
    "page": "Simulation Internals",
    "title": "Simulation.signal",
    "category": "function",
    "text": "signal(guide_data; log2fc_col)\n\n\nComputes the signal in an experiment, where the experimental signal is defined to be the average of signal of true hit genes. That value, the true hit gene signal, is the average ratio of the log2 fold change for a guide targeting a specific gene to the guide\'s theoretical phenotype.\n\nfrac1N_true times k sum_i=1^N_true sum_j=1^k\nfraclog_2 fc_ijtexttheo phenotype_ij\n\nwhere N_true is the number of true hit genes and k is the number of genes.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.noise",
    "page": "Simulation Internals",
    "title": "Simulation.noise",
    "category": "function",
    "text": "noise(guide_data; log2fc_col)\n\n\nComputes the noise in an experiment, where noise is defined to be the standard deviation of log2 fold change in the negative controls. If more than 2 bins are used then the name of log2 fold change information can be provided to differentiate between multiple log 2 fold changes.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.linear",
    "page": "Simulation Internals",
    "title": "Simulation.linear",
    "category": "function",
    "text": "linear(x, l)\n\n\nGiven value(s) x apply a simple linear function with maximum value l\n\n\n\n"
},

{
    "location": "internals.html#Simulation.sigmoid",
    "page": "Simulation Internals",
    "title": "Simulation.sigmoid",
    "category": "function",
    "text": "sigmoid(x, l, k, p)\n\n\nGiven value(s) x apply the sigmoidal function with maximum value l, a steepness k, and an inflection point p.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.auprc",
    "page": "Simulation Internals",
    "title": "Simulation.auprc",
    "category": "function",
    "text": "auprc(scores, classes, pos_labels; rev)\n\n\nComputes the area under the Precision-Recall curve using a lower trapezoidal estimator, which is more accurate for skewed datasets.\n\nK. Boyd, K. H. Eng, and C. D. Page, “Area under the Precision-Recall Curve: Point Estimates and Confidence Intervals,” in Machine Learning and Knowledge Discovery in Databases, H. Blockeel, K. Kersting, S. Nijssen, and F. Železný, Eds. Springer Berlin Heidelberg, 2013, pp. 451–466.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.auroc",
    "page": "Simulation Internals",
    "title": "Simulation.auroc",
    "category": "function",
    "text": "auroc(scores, classes, pos_labels; rev)\n\n\nOptimized function for computing the area under the receiver operator characteristic curve.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.venn",
    "page": "Simulation Internals",
    "title": "Simulation.venn",
    "category": "function",
    "text": "venn(scores, classes, pos_labels; rev)\n\n\nGiven N positive examples, computes the percentage of the top N/2 hits that are correct\n\n\n\n"
},

{
    "location": "internals.html#Miscellaneous-functions-1",
    "page": "Simulation Internals",
    "title": "Miscellaneous functions",
    "category": "section",
    "text": "Simulation.signal\nSimulation.noise\nSimulation.linear\nSimulation.sigmoid\nSimulation.auprc\nSimulation.auroc\nSimulation.venn"
},

]}
