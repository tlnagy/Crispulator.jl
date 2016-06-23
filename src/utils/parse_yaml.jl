const param_details = Array[[:librarygenome_num_genes, :num_genes, Int64, :both, :screen],
    [:screen_num_bottlenecks, :num_bottlenecks, Int64, :growth, :screen],
    [:screen_type, :screen_type, AbstractString, :both, :screen],
    [:screenrepresentation_selection, :bottleneck_representation, Int64, :both, :screen],
    [:libraryguides_frac_high_quality, :frac_hq, Float64, :both, :library],
    [:screenrepresentation_transfection, :representation, Int64, :both, :screen],
    [:screen_std_noise, :σ, Float64, :facs, :screen],
    [:librarygenome_num_guides_per_gene, :coverage, Int64, :both, :screen],
    [:librarygenome_frac_decreasing_genes, :frac_dec_genes, Float64, :both, :library],
    [:libraryguides_crispr_type, :crisprtype, AbstractString, :both, :library],
    [:libraryguides_mean_high_quality_kd, :mean_hq_kd, Float64, :both, :library],
    [:screen_bin_size, :bin_info, Float64, :facs, :screen],
    [:librarygenome_frac_increasing_genes, :frac_inc_genes, Float64, :both, :library],
    [:screenrepresentation_sequencing, :seq_depth, Int64, :both, :screen]]

function Base.parse(data::Dict{Any, Any})

    input_dict = unravel(data)

    param_table = DataFrame(hcat(param_details...)')
    names!(param_table, [:yaml_name, :sim_name, :type, :relevance, :usage])

    input_mat = [[key, value] for (key, value) in input_dict]
    input_table = DataFrame(hcat(input_mat...)')
    names!(input_table, [:yaml_name, :input])

    combined = join(param_table, input_table, on=:yaml_name)

    println(combined)

    # get screen type
    screentype = symbol(lowercase(combined[combined[:yaml_name] .== :screen_type, :input][1]))

    # get crispr type
    crisprtype = symbol(lowercase(combined[combined[:yaml_name] .== :libraryguides_crispr_type, :input][1]))

    screen = screentype == :facs ? FacsScreen() : GrowthScreen()

    screen_fields = combined[(combined[:usage] .== :screen) &
       BitArray(map(x->x in [screentype, :both], combined[:relevance])), :]
    deleterows!(screen_fields, find(screen_fields[:yaml_name] .== :screen_type)[1])

    for row in eachrow(screen_fields)
        if row[:sim_name] == :bin_info
            p::Float64 = row[:input]
            input_val = Dict{Symbol, Tuple{Float64, Float64}}(:bin1 => (0, p), :bin2 => (1-p, 1))
        else
            input_val = convert(row[:type], row[:input])
        end
        setfield!(screen, row[:sim_name], input_val)
    end

    if lowercase(crisprtype) == "crispri"
        cas9_behavior = CRISPRi()
        knockdown_dist = Dict{Symbol, Tuple{Float64, Sampleable}}(
            :high => (frac_hq, TruncatedNormal(mean_kd, 0.1, 0, 1)),
            :low => (1-frac_hq, TruncatedNormal(0.05, 0.07, 0, 1))
        )
    elseif lowercase(crisprtype) == "crisprn"
        cas9_behavior = CRISPRKO()
        knockdown_dist = Dict{Symbol, Tuple{Float64, Sampleable}}(
            :high => (frac_hq, Delta(1.0)),
            :low => (1-frac_hq, TruncatedNormal(0.05, 0.07, 0, 1))
        )
    else
        error("CRISPR type must be either CRISPRi or CRISPRn")
    end

    max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
        :inactive => (1-frac_inc_genes-frac_dec_genes-0.05, Delta(0.0)),
        :negcontrol => (0.05, Delta(0.0)),
        :increasing => (frac_inc_genes, TruncatedNormal(0.1, 0.1, 0.025, 1)),
        :decreasing => (frac_dec_genes, TruncatedNormal(-0.55, 0.2, -1, -0.1))
    )
    lib = Library(max_phenotype_dists, knockdown_dist, cas9_behavior)

    println(data["screen"])

end

"""
unravel(data::Dict)

Hacky code that unravels the dictionary returned from the YAML
parser and flattens it.
"""
function unravel(data::Dict)
    new_data = Dict{Symbol, Any}()
    for (k1, v1) in data
        if typeof(v1) <: Dict
            for (k2, v2) in unravel(v1)
                new_data[symbol(k1, k2)] = v2
            end
        elseif typeof(v1) <: Array
            for item in v1
                for (k2, v2) in unravel(item)
                    new_data[symbol(k1, k2)] = v2
                end
            end
        else
            new_data[symbol("_", replace(k1, "-", "_"))] = v1
        end
    end
    new_data
end
