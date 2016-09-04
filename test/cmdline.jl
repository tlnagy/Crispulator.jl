# ensure consistency in the command line interface

include(normpath(joinpath(Base.source_dir(),"..","src", "parsing.jl")))

let s = build_arg_table()

    test_arg(args) = ArgParse.parse_args(args, s)

    @test test_arg(["ls"]) == Dict{AbstractString, Any}("%COMMAND%" => "ls", "ls" => Dict{AbstractString, Any}())
    @test test_arg(["config", "CONFIG", "OUTPUT"]) ==
        Dict{AbstractString,Any}("%COMMAND%"=>"config","config"=>Dict{AbstractString,Any}("config_file"=>"CONFIG","output_dir"=>"OUTPUT"))
    @test test_arg(["exp", "-p", "1", "EXP", "OUTPUT"]) ==
        Dict{AbstractString,Any}("exp"=>Dict{AbstractString,Any}("output_file"=>"OUTPUT","analysis_file"=>"EXP","addprocs"=>1,"debug"=>false),"%COMMAND%"=>"exp")
    @test test_arg(["exp", "EXP", "OUTPUT"]) ==
        Dict{AbstractString,Any}("exp"=>Dict{AbstractString,Any}("output_file"=>"OUTPUT","analysis_file"=>"EXP","addprocs"=>0,"debug"=>false),"%COMMAND%"=>"exp")

end
