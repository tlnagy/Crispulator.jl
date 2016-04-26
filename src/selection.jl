binprob(p::Float64) = (.3*p+.5, .5-.3*p)

"""
Given `cells`, a vector of integers, and `guides`, a vector of barcodes
performs simulated facs sorting of the cells.
"""
function facs_sort_hinting(cells::Vector{Int64}, guides::Vector{Barcode})
    bins = Array{Int64}[Int64[], Int64[]]
    for cell in cells
        choices = binprob(guides[cell].theo_phenotype, prob_prealloc)
        push!(bins[rand(Categorical([choices[1], choices[2]]))], cell)
    end
    bins
end
