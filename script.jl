include("sound.jl")
tile_sizes = [64*ones(Int, 4), 64*ones(Int, 4), 64*ones(Int, 2)]
run_model(1001:500000, tile_sizes)
