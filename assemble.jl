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
    w = Array{Float64}(undef, sum(tile_sizes[1]), sum(tile_sizes[2]))
    ϕ = Array{Float64}(undef, sum(tile_sizes[1]), sum(tile_sizes[2]))
    b = Array{Float64}(undef, sum(tile_sizes[1]), sum(tile_sizes[2]))
    ωy = Array{Float64}(undef, sum(tile_sizes[1]), sum(tile_sizes[2]))
    for i in 1:ni
      for j in 1:nj
	irange, jrange = tile_range(i, j, tile_sizes)
	open(@sprintf("data/u/%1d_%1d_%010d", i, j, k), "r") do file
	  u[irange,jrange] = [read(file, Float64) for ii in irange, jj in jrange]
	end
	open(@sprintf("data/w/%1d_%1d_%010d", i, j, k), "r") do file
	  w[irange,jrange] = [read(file, Float64) for ii in irange, jj in jrange]
	end
	open(@sprintf("data/ϕ/%1d_%1d_%010d", i, j, k), "r") do file
	  ϕ[irange,jrange] = [read(file, Float64) for ii in irange, jj in jrange]
	end
	open(@sprintf("data/b/%1d_%1d_%010d", i, j, k), "r") do file
	  b[irange,jrange] = [read(file, Float64) for ii in irange, jj in jrange]
	end
	open(@sprintf("data/ωy/%1d_%1d_%010d", i, j, k), "r") do file
	  ωy[irange,jrange] = [read(file, Float64) for ii in irange, jj in jrange]
	end
      end
    end
    m = 1.5e-1
    imsave(@sprintf("fig/u/%010d.png", k), Array(u'), origin="lower", vmin=-m, vmax=m, cmap="RdBu_r")
    imsave(@sprintf("fig/w/%010d.png", k), Array(w'), origin="lower", vmin=-m, vmax=m, cmap="RdBu_r")
    imsave(@sprintf("fig/ϕ/%010d.png", k), Array(ϕ'), origin="lower")
    figure(figsize=(9.6, 4.8))
    PyPlot.axes(aspect=1)
    contour(Array(b'), levels=0 : 1e-8*4000 : 1e-6*4000)
    savefig(@sprintf("fig/b/%010d.png", k), dpi=300)
    close()
    #imsave(@sprintf("fig/b/%010d.png", k), Array(b'), origin="lower")
    imsave(@sprintf("fig/ωy/%010d.png", k), Array(ωy'), origin="lower")
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
