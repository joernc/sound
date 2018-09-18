using Printf
using PyPlot

# get index ranges for tiles
tile_range(i, j, tile_sizes) = [sum(tile_sizes[1][1:i-1])+1:sum(tile_sizes[1][1:i]),
				sum(tile_sizes[2][1:j-1])+1:sum(tile_sizes[2][1:j])]

# assemble and plot results
function assemble(steps, tile_sizes)
  ni = length(tile_sizes[1])
  nj = length(tile_sizes[2])
  for k in steps
    println(k)
    u = Array{Float64}(undef, sum(tile_sizes[1]), sum(tile_sizes[2]))
    v = Array{Float64}(undef, sum(tile_sizes[1]), sum(tile_sizes[2]))
    ϕ = Array{Float64}(undef, sum(tile_sizes[1]), sum(tile_sizes[2]))
    b = Array{Float64}(undef, sum(tile_sizes[1]), sum(tile_sizes[2]))
    ωz = Array{Float64}(undef, sum(tile_sizes[1]), sum(tile_sizes[2]))
    for i in 1:ni
      for j in 1:nj
	irange, jrange = tile_range(i, j, tile_sizes)
	open(@sprintf("sound/data/u/%1d_%1d_%010d", i, j, k), "r") do file
	  u[irange,jrange] = [read(file, Float64) for ii in irange, jj in jrange]
	end
	open(@sprintf("sound/data/v/%1d_%1d_%010d", i, j, k), "r") do file
	  v[irange,jrange] = [read(file, Float64) for ii in irange, jj in jrange]
	end
	open(@sprintf("sound/data/ϕ/%1d_%1d_%010d", i, j, k), "r") do file
	  ϕ[irange,jrange] = [read(file, Float64) for ii in irange, jj in jrange]
	end
	open(@sprintf("sound/data/b/%1d_%1d_%010d", i, j, k), "r") do file
	  b[irange,jrange] = [read(file, Float64) for ii in irange, jj in jrange]
	end
	open(@sprintf("sound/data/ωz/%1d_%1d_%010d", i, j, k), "r") do file
	  ωz[irange,jrange] = [read(file, Float64) for ii in irange, jj in jrange]
	end
      end
    end
    m = 2.5e-2
    imsave(@sprintf("sound/fig/u/%010d.png", k), Array(u'), origin="lower", vmin=-m, vmax=m, cmap="RdBu_r")
    imsave(@sprintf("sound/fig/v/%010d.png", k), Array(v'), origin="lower", vmin=-m, vmax=m, cmap="RdBu_r")
    imsave(@sprintf("sound/fig/ϕ/%010d.png", k), Array(ϕ'), origin="lower")
    imsave(@sprintf("sound/fig/b/%010d.png", k), Array(b'), origin="lower")
    imsave(@sprintf("sound/fig/ωz/%010d.png", k), Array(ωz'), origin="lower")
#    # print conservation diagnostics
#    println(@sprintf("%16.10e", mass(ϕ, maski, maskx, masky)))
#    println(@sprintf("%16.10e", energy(u, v, ϕ, b, y, maski, maskx, masky)))
#    println(@sprintf("%16.10e", buoyancy(ϕ, b, maski, maskx, masky)))
#    println(@sprintf("%16.10e", buoyancy_variance(ϕ, b, maski, maskx, masky)))
  end
end

# mass
mass(ϕ, maski, maskx, masky) = sum(ϕ[maski]) + sum(ϕ[(maskx.>1).|(masky.>1)])/2

# energy
function energy(u, v, ϕ, b, y, maski, maskx, masky)
  α = ϕ.*b/c^2
  return sum(u[maski].^2)/2 + sum(v[maski].^2)/2 + sum(ϕ[maski].^2)/2c^2 + sum(ϕ[(maskx.>1).|(masky.>1)].^2)/4c^2 - 2sum(α[maski].*y[maski]) - sum(α[masky.>1].*y[masky.>1])
end

# buoyancy
function buoyancy(ϕ, b, maski, maskx, masky)
  maskb = (maskx .> 1) .| (masky .> 1)
  return sum(ϕ[maski].*b[maski]) + sum(ϕ[maskb].*b[maskb])/2
end

# buoyancy variance
function buoyancy_variance(ϕ, b, maski, maskx, masky)
  maskb = (maskx .> 1) .| (masky .> 1)
  return sum(ϕ[maski].*b[maski].^2) + sum(ϕ[maskb].*b[maskb].^2)/2
end