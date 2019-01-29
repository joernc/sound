include("slope_sound.jl")
tile_sizes = tile_sizes = [1*ones(Int, 1), 25*ones(Int, 32), 375*ones(Int, 1)]
run_model(250001:350000, tile_sizes)
