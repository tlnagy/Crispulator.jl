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
    [:libraryguides_mean_high_quality_kd, :mean_kd, Float64, :both, :library],
    [:screen_bin_size, :bin_info, Float64, :facs, :screen],
    [:librarygenome_frac_increasing_genes, :frac_inc_genes, Float64, :both, :library],
    [:screenrepresentation_sequencing, :seq_depth, Int64, :both, :screen],
    [:screen_num_runs, :num_runs, Int64, :both, :library]]

# TODO: This is shitty. A reworking of this might be nice at some point
function Base.parse(data::Dict{Any, Any})

    info("Parsing config")

    input_dict = unravel(data)

    param_table = DataFrame(permutedims(hcat(param_details...), [2, 1]))
    names!(param_table, [:yaml_name, :sim_name, :type, :relevance, :usage])

    input_mat = [Any[key, value] for (key, value) in input_dict]
    input_table = DataFrame(permutedims(hcat(input_mat...), [2, 1]))
    names!(input_table, [:yaml_name, :input])

    combined = join(param_table, input_table, on=:yaml_name)

    # get screen type
    screentype = Symbol(lowercase(combined[combined[:yaml_name] .== :screen_type, :input][1]))

    # get crispr type
    crisprtype = Symbol(lowercase(combined[combined[:yaml_name] .== :libraryguides_crispr_type, :input][1]))

    screen = screentype == :facs ? FacsScreen() : GrowthScreen()

    screen_fields = combined[(combined[:usage] .== :screen) .& BitArray(map(x->x in [screentype, :both], combined[:relevance])), :]
    deleterows!(screen_fields, find(screen_fields[:yaml_name] .== :screen_type)[1])

    lib_vals = Dict()
    curr_field, curr_type = :none, Int64

    try
        for row in eachrow(screen_fields)
            curr_field, curr_type = row[:yaml_name], row[:type]
            if row[:sim_name] == :bin_info
                p::Float64 = row[:input]
                input_val = Dict{Symbol, Tuple{Float64, Float64}}(:bin1 => (0, p), :bin2 => (1-p, 1))
            else
                input_val = convert(curr_type, row[:input])
            end
            setfield!(screen, row[:sim_name], input_val)
        end

        for row in eachrow(combined[combined[:usage] .== :library, :])
            curr_field, curr_type = row[:yaml_name], row[:type]
            (row[:sim_name] == :crisprtype) && continue
            lib_vals[row[:sim_name]] = row[:input]::curr_type
        end
    catch e
        if isa(e, TypeError) || isa(e, InexactError)
            println("ERROR: $curr_field must have a type of $curr_type")
            exit(1)
        end
        rethrow(e)
    end

    if crisprtype == :crispri
        cas9_behavior = CRISPRi()
        knockdown_dist = Dict{Symbol, Tuple{Float64, Sampleable}}(
            :high => (lib_vals[:frac_hq], TruncatedNormal(lib_vals[:mean_kd], 0.1, 0, 1)),
            :low => (1-lib_vals[:frac_hq], TruncatedNormal(0.05, 0.07, 0, 1))
        )
    elseif crisprtype == :crisprn
        cas9_behavior = CRISPRn()
        knockdown_dist = Dict{Symbol, Tuple{Float64, Sampleable}}(
            :high => (lib_vals[:frac_hq], Delta(1.0)),
            :low => (1-lib_vals[:frac_hq], TruncatedNormal(0.05, 0.07, 0, 1))
        )
    else
        error("CRISPR type must be either CRISPRi or CRISPRn")
    end

    frac_inc_genes, frac_dec_genes = lib_vals[:frac_inc_genes], lib_vals[:frac_dec_genes]
    max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
        :inactive => (1-frac_inc_genes-frac_dec_genes-0.2, Delta(0.0)),
        :negcontrol => (0.2, Delta(0.0)),
        :increasing => (frac_inc_genes, TruncatedNormal(0.1, 0.1, 0.025, 1)),
        :decreasing => (frac_dec_genes, TruncatedNormal(-0.55, 0.2, -1, -0.1))
    )
    lib = Library(max_phenotype_dists, knockdown_dist, cas9_behavior)

    (screen, lib, lib_vals[:num_runs])
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
                new_data[Symbol(k1, k2)] = v2
            end
        elseif typeof(v1) <: Array && typeof(v1[1]) <: Dict
            for item in v1
                for (k2, v2) in unravel(item)
                    new_data[Symbol(k1, k2)] = v2
                end
            end
        else
            new_data[Symbol("_", replace(k1, "-", "_"))] = v1
        end
    end
    new_data
end
