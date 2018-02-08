using Compat

"""
Any entity that is tracked through the pooled experiment. For CRISPR screens,
this is equivalent to the sgRNA. This object stores properties relating to the
performance of this entity in the screen.

$(FIELDS)
"""
type Barcode
    "The target gene id"
    gene::Int

    "The knockdown efficiency"
    knockdown::Float64

    "The theoretical phenotype exerted"
    theo_phenotype::Float64

    "The gene behavior, either sigmoidal or linear"
    behavior::Symbol

    "The gene activity class, either inactive, up or down"
    class::Symbol

    "Initial frequency after transfection"
    initial_freq::Float64

    function Barcode(gene, knockdown, phenotype, behavior, class)
        if (0 <= knockdown <= 1)
            new(gene, knockdown, phenotype, behavior, class, -Inf)
        else
            error("Knockdown must be between 0 and 1, inclusive")
        end
    end
end

as_array(bc::Barcode) = [getfield(bc, fld) for fld in fieldnames(bc)]

"""
A description of screen parameters
"""
@compat abstract type ScreenSetup end

as_array(ss::ScreenSetup) = [getfield(ss, fld) for fld in fieldnames(ss)]

"""
A type representing the parameters used in a typical FACS screen.

$(FIELDS)
"""
type FacsScreen <: ScreenSetup
    "Number of genes targeted by the screen"
    num_genes::Int

    "Number of guides per gene"
    coverage::Int

    """Number of cells with each guide. `representation=10` means that there are
    10 times as many cells as guides during transfection. The actual number of
    cells per guide post-transfection will be less depending on the MOI"""
    representation::Int

    """The multiplicity of infection, ``\\lambda``, of the screen. We model this as a Poisson
    process during transfection (see [`Simulation.transfect`](@ref)).

    !!! note

        We **do not** model multiple infections. We assume that the MOI is properly
        selected and less than half the cells are transfected by any virus, i.e.
        ``\\lambda \\lt 0.5`` and then select only the cells that have a single
        transfection occurrence:

        ```math
        P(x = 1; Poisson(λ))
        ```

        For ``\\lambda = 0.25`` this works out to being ``\\approx 19.5\\%`` of the
        number of cells (`num_genes` * `coverage` * `representation`).
    """
    moi::Float64

    """The standard deviation expected for cells during FACS sorting. This should
    be set according to the biological variance experimentally observed, e.g. in
    the fluorescence intensity of isogenic cells"""
    σ::Float64

    """Range of guide phenotypes to collect in each bin

    # Example

    In the following example

    ```julia
    p = 0.05
    bin_info = Dict(:bin1 => (0.0, p), :bin2 => (1.0-p, 1.0))
    ```

    The 5th percentile of cells sorted according to their phenotype (fluorescence,
    size, etc) will be compared to the 95th percentile.
    """
    bin_info::Dict{Symbol, Tuple{Float64, Float64}}

    "Number of cells sorted expressed as an integer multiple of the number of guides"
    bottleneck_representation::Int

    """Sequencing depth as a integer multiple of the number of guides, i.e.
    `seq_depth=10` is equivalent to `10 * num_genes * coverage` reads.
    """
    seq_depth::Int

    function FacsScreen()
        new(500, 5, 100, 0.25, 1.0, Dict(:bin1 => (0.0, 1/3), :bin2 => (2/3, 1.0)), 1000, 1000)
    end
end

"""
A type representing the parameters used in a typical growth-based screen.

$(FIELDS)
"""
type GrowthScreen <: ScreenSetup
    "Number of genes targeted by the screen"
    num_genes::Int

    "Number of guides per gene"
    coverage::Int

    """Number of cells with each guide. `representation=10` means that there are
    10 times as many cells as guides during transfection. The actual number of
    cells per guide post-transfection will be less depending on the MOI"""
    representation::Int

    """The multiplicity of infection, ``\\lambda``, of the screen. We model this as a Poisson
    process during transfection (see [`Simulation.transfect`](@ref)).

    !!! note

        We **do not** model multiple infections. We assume that the MOI is properly
        selected and less than half the cells are transfected by any virus, i.e.
        ``\\lambda \\lt 0.5`` and then select only the cells that have a single
        transfection occurrence:

        ```math
        P(x = 1; Poisson(λ))
        ```

        For ``\\lambda = 0.25`` this works out to being ``\\approx 19.5\\%`` of the
        number of cells (`num_genes` * `coverage` * `representation`).
    """
    moi::Float64

    """Sequencing depth as a integer multiple of the number of guides, i.e.
    `seq_depth=10` is equivalent to `10 * num_genes * coverage` reads.
    """
    seq_depth::Int

    """For growth screens, how much of a bottleneck is applied. This is minimum
    number of cells that is kept when passaging the pool of cells. This is
    expressed as an integer multiple of the number of guides."""
    bottleneck_representation::Int

    """For growth screens, how many bottlenecks are applied. This is the integer
    number of passages during the growth screen."""
    num_bottlenecks::Int

    """Before each passage, the theoretical phenotype of each cell is convolved
    with a normal noise distribution with a standard deviation, σ, dictated by
    this parameter. This should be set based on an expectation of noisiness of
    the subsampling"""
    noise::Float64

    function GrowthScreen()
        new(500, 5, 100, 0.25, 1000, 1000, 10, 0.01)
    end
end
