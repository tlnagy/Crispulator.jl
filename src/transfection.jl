function transfect(guide_freqs_dist::Distributions.Categorical, cell_count::Int64, moi::Float64)
    rand(guide_freqs_dist, round(Int64, pdf(Poisson(moi), 1)*cell_count))
end
