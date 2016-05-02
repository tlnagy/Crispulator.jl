using StatsBase

function counts_to_freqs(bin_cells::Dict{Symbol, Vector{Int64}})
    results = Dict{Symbol, Vector{Float64}}()
    for (bin, cells) in bin_cells
        counts = StatsBase.counts(cells, 1:N*coverage)
        results[bin] = counts ./ sum(counts)
    end
    results
end

function sequencing(depths::Dict{Symbol, Int64}, guides::Vector{Barcode}, samples::Dict{Symbol, Vector{Float64}})
    @assert all(Bool[haskey(depths, key) for key in keys(samples)]) "Supply exactly one sequencing depth per sample"
    sequencing_results = Dict{Symbol, DataFrame}()
    colnames = fieldnames(Barcode)
    push!(colnames, [:counts, :barcodeid]...)
    coltypes = map(x -> fieldtype(Barcode, x), fieldnames(Barcode))
    push!(coltypes, [Int64, Int64]...)
    num_guides = length(guides)

    for (bin, freqs) in samples
        reads = rand(Categorical(freqs), depths[bin])

        raw_counts = StatsBase.counts(reads, 1:num_guides)

        bc_fieldnames = fieldnames(Barcode)
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
