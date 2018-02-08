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
    "text": "To run the simulation, you will need a recent version of Julia installed and in your PATH. Then navigate into the root directory of the project and run julia. Run the following command:julia -e 'Pkg.update(); Pkg.add(\"Crispulator\")'this copies Crispulator over to the Julia package directory and installs all of its dependencies."
},

{
    "location": "index.html#Quickstart-1",
    "page": "Home",
    "title": "Quickstart",
    "category": "section",
    "text": "From the root directory of the project runjulia src/run.jl config example_config.yml ."
},

{
    "location": "index.html#Advanced-1",
    "page": "Home",
    "title": "Advanced",
    "category": "section",
    "text": "See the Simulation Internals page for in-depth documentation needed for more advanced usage"
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
    "category": "Type",
    "text": "A type representing the parameters used in a typical FACS screen.\n\nnum_genes\nNumber of genes targeted by the screen\ncoverage\nNumber of guides per gene\nrepresentation\nNumber of cells with each guide. representation=10 means that there are 10 times as many cells as guides during transfection. The actual number of cells per guide post-transfection will be less depending on the MOI\nmoi\nThe multiplicity of infection, lambda, of the screen. We model this as a Poisson process during transfection (see Simulation.transfect).\nnote: Note\nWe do not model multiple infections. We assume that the MOI is properly selected and less than half the cells are transfected by any virus, i.e. lambda lt 05 and then select only the cells that have a single transfection occurrence:P(x=1Poisson())For lambda = 025 this works out to being approx 195 of the number of cells (num_genes * coverage * representation).\n\nσ\nThe standard deviation expected for cells during FACS sorting. This should be set according to the biological variance experimentally observed, e.g. in the fluorescence intensity of isogenic cells\nbin_info\nRange of guide phenotypes to collect in each bin\nExample\nIn the following example\np = 0.05\nbin_info = Dict(:bin1 => (0.0, p), :bin2 => (1.0-p, 1.0))\nThe 5th percentile of cells sorted according to their phenotype (fluorescence, size, etc) will be compared to the 95th percentile.\n\nbottleneck_representation\nNumber of cells sorted expressed as an integer multiple of the number of guides\nseq_depth\nSequencing depth as a integer multiple of the number of guides, i.e. seq_depth=10 is equivalent to 10 * num_genes * coverage reads.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.GrowthScreen",
    "page": "Simulation Internals",
    "title": "Simulation.GrowthScreen",
    "category": "Type",
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
    "category": "Function",
    "text": "construct_library(setup, lib)\n\n\nConstructs the guide library for N genes with coverage number of guides per gene. Returns a tuple of guides and their relative frequencies (assigned randomly).\n\n\n\n"
},

{
    "location": "internals.html#Simulation.transfect",
    "page": "Simulation Internals",
    "title": "Simulation.transfect",
    "category": "Function",
    "text": "transfect(setup, lib, guides, guide_freqs_dist)\n\n\nSimulates the transfection the of the library of guides (guides) according to the parameters in Simulation.Library lib and the frequencies represented by the Categorical distribution guide_freqs_dist.\n\nThis function creates a population of cells post-transfection where each cell has single guide RNA present according to the frequencies of the guides in the present in the library. The cells are then grown up so that there are setup.bottleneck_representation * num_guides of them. For FACS screens we assume a simple linear expansion of the cells, while for Growth screens the cells grow according as a function of their phenotype and degree of knockdown.\n\n\n\n"
},

{
    "location": "internals.html#Base.Sort.select",
    "page": "Simulation Internals",
    "title": "Base.Sort.select",
    "category": "Function",
    "text": "select(setup, cells, cell_phenotypes, guides; debug)\n\n\nGiven cells, a vector of integers, and guides, a vector of barcodes, performs simulated FACS sorting of the cells into bins with the given cutoffs. σ refers to the spread of the phenotype during the FACS screen.\n\n\n\nselect(setup, initial_cells, initial_cell_phenotypes, guides; debug)\n\n\nGrowth Screen selection\n\n\n\n"
},

{
    "location": "internals.html#Simulation.sequencing",
    "page": "Simulation Internals",
    "title": "Simulation.sequencing",
    "category": "Function",
    "text": "sequencing(depths, guides, samples)\n\n\nSimulates next-gen sequencing by generating a Categorical distribution based on frequencies of each guide in each bin and then randomly samples from this distributions up to the depth provided in depths for each bin. Returns a dict of DataFrames with the bin names as the keys. This object can then be passed through Simulation.counts_to_freqs followed by Simulation.differences_between_bins.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.counts_to_freqs",
    "page": "Simulation Internals",
    "title": "Simulation.counts_to_freqs",
    "category": "Function",
    "text": "counts_to_freqs(bin_cells, guide_count)\n\n\nConverts raw counts to frequencies by dividing by the total number of reads for each sample.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.differences_between_bins",
    "page": "Simulation Internals",
    "title": "Simulation.differences_between_bins",
    "category": "Function",
    "text": "differences_between_bins(raw_data; first_bin, last_bin)\n\n\nGiven the raw data from Simulation.sequencing returns two DataFrames\n\nguide_data: This DataFrame contains the per-guide level data including the  log2 fold change in the normalized frequencies of each guide between the two  bins.\ngene_data: This DataFrame contains the same information but grouped by  gene. The log2 fold change data from the first DataFrame is used to calculate  the average log2 fold change per gene and a pvalue computed using a  Mann-Whitney U-test as  measure of how consistently shifted the guides are of this gene versus the  population of negative control guides. (see below for more info)\n\nA typical guide_data DataFrame contains the following columns:\n\ngene: the gene ID of that this guide targets\nknockdown: activity of the guide on 0 to 1 scale, where 1 is complete knockout\nbarcodeid: the ID of this specific guide\ntheo_phenotype: expected phenotype of this guide, generally a -1 to 1 scale\nbehavior: whether the target gene displays a linear or sigmoidal response to knockdown\nclass: whether the target gene has a positive, negative, or no phenotype during screening\ninitial_freq: frequency of guide post-transfection (see Simulation.transfect)\ncounts_bin1: the raw number of reads for each guide in the first bin\nfreqs_bin1: the number of reads for each guide divided by the total number of reads in this bin\nrel_freqs_bin1: the frequency of each guide divided by the median frequency of negative control guides\ncounts_bin2: the raw number of reads for each guide in the second bin\nfreqs_bin2: the number of reads for each guide divided by the total number of reads in this bin\nrel_freqs_bin2: the frequency of each guide divided by the median frequency of negative control guides for this bin\nlog2fc_bin2: the log2 fold change in relative guide frequencies between the two bins\n\nA typical gene_data DataFrame contains the following data:\n\ngene: this gene's ID\nbehavior: whether this gene displays a linear or sigmoidal response to knockdown\nclass: whether this gene has a positive, negative, or no phenotype during screening\nmean: the mean log 2 fold change in relative frequencies between the two bins   for all the guides targeting this gene.\npvalue: the -log10 pvalue of the log2 fold changes of all guides targeting   this gene as computed by the non-parametric Mann-Whitney U-test. A measure   of the consistency of the log 2 fold changes[1]\nabsmean: absolute value of mean per-gene\npvalmeanprod: mean multiplied with the pvalue per-gene\n\nFurther reading\n\n[1]: Kampmann M, Bassik MC, Weissman JS. Integrated platform for genome-wide screening and construction of high-density genetic interaction maps in mammalian cells. Proc Natl Acad Sci U S A. 2013;110:E2317–26.\n\n\n\n"
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
    "category": "Type",
    "text": "Wrapper containing all library construction parameters\n\nknockdown_dist\nDistribution of guide knockdown efficiencies\nknockdown_probs\nmax_phenotype_dists\nMaximum phenotype categories mapped to their probability of being selected and the distribution to draw from if they are selected\n\nphenotype_probs\nkd_phenotype_relationships\nKnockdown-phenotype relationships mapped to their probability of being selected and their respective KDPhenotypeRelationship\n\nrelationship_probs\ncas9_behavior\nWhether this library is CRISPRi or CRISPR cutting.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.Barcode",
    "page": "Simulation Internals",
    "title": "Simulation.Barcode",
    "category": "Type",
    "text": "Any entity that is tracked through the pooled experiment. For CRISPR screens, this is equivalent to the sgRNA. This object stores properties relating to the performance of this entity in the screen.\n\ngene\nThe target gene id\nknockdown\nThe knockdown efficiency\ntheo_phenotype\nThe theoretical phenotype exerted\nbehavior\nThe gene behavior, either sigmoidal or linear\nclass\nThe gene activity class, either inactive, up or down\ninitial_freq\nInitial frequency after transfection\n\n\n\n"
},

{
    "location": "internals.html#Simulation.KDPhenotypeRelationship",
    "page": "Simulation Internals",
    "title": "Simulation.KDPhenotypeRelationship",
    "category": "Type",
    "text": "A type representing a relationship between degree of knockdown and effect on phenotype\n\n\n\n"
},

{
    "location": "internals.html#Simulation.Cas9Behavior",
    "page": "Simulation Internals",
    "title": "Simulation.Cas9Behavior",
    "category": "Type",
    "text": "A type representing the behavior of different Cas9s\n\n\n\n"
},

{
    "location": "internals.html#Simulation.CRISPRn",
    "page": "Simulation Internals",
    "title": "Simulation.CRISPRn",
    "category": "Type",
    "text": "type CRISPRn <: Simulation.Cas9Behavior\n\nCRISPR KO behavior is more complex since sgRNA-directed DNA damage repair is stochastic. We assume that 2/3 of repair events at a given locus lead to a frameshift, and that the screen is carried out in diploid cells. The assumption that only bi-allelic frame-shift mutations lead to a phenotype in CRISPRn screens for most sgRNAs is supported by the empirical finding that in-frame deletions mostly do not show strong phenotypes, unless they occur in regions encoding conserved residues or domains[2]\n\n[2]: Horlbeck MA, Gilbert LA, Villalta JE, Adamson B, Pak RA, Chen Y, Fields AP, Park CY, Corn JE, Kampmann M, Weissman JS: Compact and highly active next- generation libraries for CRISPR-mediated gene repression and activation. Elife 2016, 5.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.CRISPRi",
    "page": "Simulation Internals",
    "title": "Simulation.CRISPRi",
    "category": "Type",
    "text": "type CRISPRi <: Simulation.Cas9Behavior\n\nCRISPRi behavior is simply determined by the activity of the guide\n\n\n\n"
},

{
    "location": "internals.html#Miscellaneous-Types-1",
    "page": "Simulation Internals",
    "title": "Miscellaneous Types",
    "category": "section",
    "text": "Simulation.Library\nSimulation.Barcode\nSimulation.KDPhenotypeRelationship\nSimulation.Cas9Behavior\nSimulation.CRISPRn\nSimulation.CRISPRi"
},

{
    "location": "internals.html#Simulation.signal",
    "page": "Simulation Internals",
    "title": "Simulation.signal",
    "category": "Function",
    "text": "signal(bc_counts)\n\n\nComputes the signal in an experiment, where the experimental signal is defined to be the average of signal of true hit genes. That value, the true hit gene signal, is the average ratio of the log2 fold change for a guide targeting a specific gene to the guide's theoretical phenotype.\n\nfrac1N_true times k sum_i=1^N_true sum_j=1^k\nfraclog_2 fc_ijtexttheo phenotype_ij\n\nwhere N_true is the number of true hit genes and k is the number of genes.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.noise",
    "page": "Simulation Internals",
    "title": "Simulation.noise",
    "category": "Function",
    "text": "noise(bc_counts)\n\n\nComputes the noise in an experiment, where noise is defined to be the standard deviation of log2 fold change in the negative controls\n\n\n\n"
},

{
    "location": "internals.html#Simulation.linear",
    "page": "Simulation Internals",
    "title": "Simulation.linear",
    "category": "Function",
    "text": "linear(x, l)\n\n\nGiven value(s) x apply a simple linear function with maximum value l\n\n\n\n"
},

{
    "location": "internals.html#Simulation.sigmoid",
    "page": "Simulation Internals",
    "title": "Simulation.sigmoid",
    "category": "Function",
    "text": "sigmoid(x, l, k, p)\n\n\nGiven value(s) x apply the sigmoidal function with maximum value l, a steepness k, and an inflection point p.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.auprc",
    "page": "Simulation Internals",
    "title": "Simulation.auprc",
    "category": "Function",
    "text": "auprc(scores, classes, pos_labels; rev)\n\n\nComputes the area under the Precision-Recall curve using a lower trapezoidal estimator, which is more accurate for skewed datasets.\n\nK. Boyd, K. H. Eng, and C. D. Page, “Area under the Precision-Recall Curve: Point Estimates and Confidence Intervals,” in Machine Learning and Knowledge Discovery in Databases, H. Blockeel, K. Kersting, S. Nijssen, and F. Železný, Eds. Springer Berlin Heidelberg, 2013, pp. 451–466.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.auroc",
    "page": "Simulation Internals",
    "title": "Simulation.auroc",
    "category": "Function",
    "text": "auroc(scores, classes, pos_labels; rev)\n\n\nOptimized function for computing the area under the receiver operator characteristic curve.\n\n\n\n"
},

{
    "location": "internals.html#Simulation.venn",
    "page": "Simulation Internals",
    "title": "Simulation.venn",
    "category": "Function",
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
