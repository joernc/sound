include("sound.jl")
tile_sizes = [128*ones(Int, 4), 128*ones(Int, 4), 192*ones(Int, 2)]
run_model(1:5000, tile_sizes)
