binprob(p::Float64) = (.3*p+.5, .5-.3*p)
binprob!(p::Float64, prealloc::Vector{Float64}) = begin prealloc[1] = .3*p+.5; prealloc[2] = .5-.3*p end

"""
Given `cells`, a vector of integers, and `guides`, a vector of barcodes
performs simulated facs sorting of the cells. Has an optional parameter
`max_guess` that is maximum expect fraction of cells that end up in one
bin for preallocation purposes.
"""
function facs_sort_hinting(cells::Vector{Int64}, guides::Vector{Barcode}; max_guess::Float64 = 0.51)
    prob_prealloc = zeros(2)
    hint = round(Int64, length(cells)*max_guess)
    bins = Array{Int64}[zeros(Int64, hint), zeros(Int64, hint)]
    bin_counts = [0, 0]
    for cell in cells
        binprob!(guides[cell].theo_phenotype, prob_prealloc)
        choice = rand(Categorical(prob_prealloc))
        bin_counts[choice] += 1
        bins[choice][bin_counts[choice]] = cell
    end
    resize!(bins[1], bin_counts[1])
    resize!(bins[2], bin_counts[2])
    bins
end
