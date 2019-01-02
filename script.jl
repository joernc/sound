include("rev_sound.jl")
tile_sizes = tile_sizes = [100*ones(Int, 8), 100*ones(Int, 8), 100*ones(Int, 2)]
run_model(1:200, tile_sizes)
