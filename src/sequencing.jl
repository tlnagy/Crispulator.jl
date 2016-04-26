function counts_to_freqs(samples::Vector{Int64}...)
    results = Vector{Float64}[]
    for sample in samples
        counts = StatsBase.counts(sample, 1:N*coverage)
        push!(results, counts ./ sum(counts))
    end
    results
end

function sequencing(depths::Vector{Int64}, guides::Vector{Barcode}, samples::Vector{Float64}...)
    @assert (length(depths) == length(samples)) "Supply exactly one sequencing depth per sample"
    sequencing_results = DataFrame[]
    colnames = fieldnames(Barcode)
    push!(colnames, [:counts, :barcodeid]...)
    coltypes = map(x -> fieldtype(Barcode, x), fieldnames(Barcode))
    push!(coltypes, [Int64, Int64]...)

    for i in 1:length(samples)
        reads = rand(Categorical(samples[i]), depths[i])

        raw_counts = StatsBase.counts(reads, 1:N*coverage)

        bc_fieldnames = fieldnames(Barcode)
        data = Any[coltype[] for coltype in coltypes]

        for id in 1:N*coverage
            barcode_info = as_array(guides[id])
            push!(barcode_info, [raw_counts[id], id]...)

            for i in 1:length(barcode_info)
                push!(data[i], barcode_info[i])
            end
        end
        push!(sequencing_results, DataFrame(data, colnames))
    end
    sequencing_results
end
