println("Using $(nprocs()) threads")
include("load.jl")

function main()
    target_function = symbol(ARGS[1])
    filename = ARGS[2]
    println("Calling $target_function and saving in $filename")
    eval(Expr(:call, target_function, filename))
end

# fire up simulation if run using command line
if !isinteractive()
    main()
end
