"""
Any entity that is tracked through the pooled experiment. For CRISPR screens,
this is equivalent to the sgRNA. This object stores properties relating to the
performance of this entity in the screen.
"""
type Barcode
    "The target gene id"
    gene::Int64

    "The knockdown efficiency"
    knockdown::Float64

    "The theoretical phenotype exerted"
    theo_phenotype::Float64

    "The observed phenotype"
    obs_phenotype::Float64

    "The gene behavior, either sigmoidal or linear"
    behavior::Symbol

    "The gene activity class, either inactive, up or down"
    class::Symbol

    "Initial frequency after transfection"
    initial_freq::Float64

    function Barcode(gene, knockdown, phenotype, behavior, class)
        if (0 <= knockdown <= 1)
            new(gene, knockdown, phenotype, -Inf, behavior, class, -Inf)
        else
            error("Knockdown must be between 0 and 1, inclusive")
        end
    end
end

as_array(bc::Barcode) = [getfield(bc, fld) for fld in fieldnames(bc)]

"""
A description of screen parameters
"""
abstract ScreenSetup

as_array(ss::ScreenSetup) = [getfield(ss, fld) for fld in fieldnames(ss)]

type FacsScreen <: ScreenSetup
    "Number of target genes"
    num_genes::Int64
    "Number of guides per gene"
    coverage::Int64
    "Number of cells with each guide"
    representation::Int64
    "Multiplicity of infection"
    moi::Float64
    "Std dev expected for cells during facs sorting"
    σ::Float64
    "Range of guide phenotypes to collect in each bin"
    bin_info::Dict{Symbol, Tuple{Float64, Float64}}
    "Number of cells sorted"
    bottleneck_representation::Int64
    "Sequencing depth"
    seq_depth::Int64

    function FacsScreen()
        new(500, 5, 100, 0.25, 1.0, Dict(:bin1 => (0.0, 1/3), :bin2 => (2/3, 1.0)), 1000, 1000)
    end
end

type GrowthScreen <: ScreenSetup
    "Number of target genes"
    num_genes::Int64
    "Number of guides per gene"
    coverage::Int64
    "Number of cells with each guide"
    representation::Int64
    "Multiplicity of infection"
    moi::Float64
    "Sequencing depth"
    seq_depth::Int64
    "For growth screens, how much of a bottleneck is applied"
    bottleneck_representation::Int64
    "For growth screens, how many bottlenecks are applied"
    num_bottlenecks::Int64
    "The σ of the normal noise distribution added to each cell's phenotype"
    noise::Float64

    function GrowthScreen()
        new(500, 5, 100, 0.25, 1000, 1000, 10, 0.01)
    end
end
