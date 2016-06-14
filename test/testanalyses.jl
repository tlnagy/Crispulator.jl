analyses_path = normpath(joinpath(Base.source_dir(),"..",joinpath("src", "analyses")))

for analysis in readdir(analyses_path)
    include(joinpath(analyses_path, analysis))
    if analysis == "common.jl"
        continue
    else
        tempfile = tempname()
        main(tempfile, debug=true, quiet=true)
        rm(tempfile)
    end
end
