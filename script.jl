include("rev_sound.jl")
#tile_sizes = [400*ones(Int, 40), 1*ones(Int, 1), 400*ones(Int, 2)]
tile_sizes = [64*ones(Int, 4), 64*ones(Int, 4), 64*ones(Int, 4)]
run_model(4001:10000, tile_sizes)
