include("sound.jl")
tile_sizes = [64*ones(Int, 8), 128*ones(Int, 4), [256]]
run_model(301:30000, tile_sizes)
