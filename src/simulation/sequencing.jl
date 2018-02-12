"""
$(SIGNATURES)

Converts raw counts to frequencies by dividing by the total number of reads
for each sample.
"""
function counts_to_freqs(bin_cells::Associative{Symbol, Vector{Int}}, guide_count::Int)
    results = Dict{Symbol, Vector{Float64}}()
    for (bin, cells) in bin_cells
        counts = StatsBase.counts(cells, 1:guide_count)
        results[bin] = counts ./ sum(counts)
    end
    results
end

"""
$(SIGNATURES)

Simulates next-gen sequencing by generating a Categorical distribution based on
frequencies of each guide in each bin and then randomly samples from this
distributions up to the depth provided in `depths` for each bin. Returns a
dict of DataFrames with the bin names as the keys. This object can then be passed
through [`Simulation.counts_to_freqs`](@ref) followed by
[`Simulation.differences_between_bins`](@ref).
"""
function sequencing(depths::Associative{Symbol, Int}, guides::Vector{Barcode}, samples::Associative{Symbol, Vector{Float64}})
    @assert all(Bool[haskey(depths, key) for key in keys(samples)]) "Supply exactly one sequencing depth per sample"
    sequencing_results = Dict{Symbol, DataFrame}()
    colnames = [fld for fld in fieldnames(Barcode) if fld != :obs_phenotype]
    coltypes = map(x -> fieldtype(Barcode, x), colnames)
    push!(colnames, [:counts, :barcodeid]...)
    push!(coltypes, [Int, Int]...)
    num_guides = length(guides)

    for (bin, freqs) in samples
        reads = rand(Categorical(freqs), depths[bin]*length(guides))

        raw_counts = StatsBase.counts(reads, 1:num_guides)

        data = Any[coltype[] for coltype in coltypes]

        # convert Barcode objects to tabular form
        for id in 1:num_guides
            barcode_info = as_array(guides[id])
            push!(barcode_info, [raw_counts[id], id]...)

            for column in 1:length(barcode_info)
                push!(data[column], barcode_info[column])
            end
        end
        sequencing_results[bin] = DataFrame(data, colnames)
    end
    sequencing_results
end
