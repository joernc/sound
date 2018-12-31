include("sound.jl")
tile_sizes = [64*ones(Int, 1), 64*ones(Int, 4), 64*ones(Int, 16)]
run_model(7001:100000, tile_sizes)
