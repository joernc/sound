using Printf
using PyPlot

# get index ranges for tiles
tile_range(i, j, tile_sizes) = [sum(tile_sizes[1][1:i-1])+1:sum(tile_sizes[1][1:i]),
				sum(tile_sizes[2][1:j-1])+1:sum(tile_sizes[2][1:j])]

function assemble(steps, tile_sizes)
  ni = length(tile_sizes[1])
  nj = length(tile_sizes[2])
  # assemble and plot results
  for k in steps
    println(k)
    u = Array{Float64}(undef, sum(tile_sizes[1]), sum(tile_sizes[2]))
    for i in 1:ni
      for j in 1:nj
	irange, jrange = tile_range(i, j, tile_sizes)
	open(@sprintf("sound/data/%1d_%1d_%010d", i, j, k), "r") do file
	  u[irange,jrange] = [read(file, Float64) for ii in irange, jj in jrange]
	end
      end
    end
    #m = maximum(abs.(u[.!isnan.(u)]))
    m = .5
    imsave(@sprintf("sound/fig/%010d.png", k), Array(u'), origin="lower", vmin=-m, vmax=m, cmap="RdBu_r")
  end
end
