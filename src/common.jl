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
    
    function Barcode(gene, knockdown, phenotype, behavior, class)
        if (0 <= knockdown <= 1)
            new(gene, knockdown, phenotype, -Inf, behavior, class)
        else
            error("Knockdown must be between 0 and 1, inclusive")
        end
    end
end

as_array(bc::Barcode) = [getfield(bc, fld) for fld in fieldnames(bc)]

type Cell
    barcode_ids::Vector{Int64}
end
