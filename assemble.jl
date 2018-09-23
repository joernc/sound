using Printf
using PyPlot
using HDF5

const Δx = 60e3/1024
const Δy = Δx
const Δz = 1000/512

# background stratification
const N = 1e-3

# slope angle
const θ = 2e-3

pygui(false)
rc("contour", negative_linestyle="solid")

# get index ranges for tiles
tile_range(i, j, k, tile_sizes) = [sum(tile_sizes[1][1:i-1])+1:sum(tile_sizes[1][1:i]),
                                   sum(tile_sizes[2][1:j-1])+1:sum(tile_sizes[2][1:j]),
                                   sum(tile_sizes[3][1:k-1])+1:sum(tile_sizes[3][1:k])]

# assemble and plot results
function assemble(steps, tile_sizes, mu, mv, mw)
  ni = length(tile_sizes[1])
  nj = length(tile_sizes[2])
  nk = length(tile_sizes[2])
  nx = sum(tile_sizes[1])
  ny = sum(tile_sizes[2])
  nz = sum(tile_sizes[3])
  x = [(i-1)*Δx for i = 1:nx, j = 1:ny, k = 1:nz]
  y = [(j-1)*Δx for i = 1:nx, j = 1:ny, k = 1:nz]
  z = [(k-1)*Δz for i = 1:nx, j = 1:ny, k = 1:nz]
  for n in steps
    println(n)
    u = Array{Float64, 3}(undef, nx, ny, nz)
    v = Array{Float64, 3}(undef, nx, ny, nz)
    w = Array{Float64, 3}(undef, nx, ny, nz)
    ϕ = Array{Float64, 3}(undef, nx, ny, nz)
    b = Array{Float64, 3}(undef, nx, ny, nz)
    for k = 1:nk, j = 1:nj, i = 1:ni
      irange, jrange, krange = tile_range(i, j, k, tile_sizes)
      filename = @sprintf("data/%010d_%1d_%1d_%1d.h5", n, i, j, k)
      u[irange,jrange,krange] = h5read(filename, "u")
      v[irange,jrange,krange] = h5read(filename, "v")
      w[irange,jrange,krange] = h5read(filename, "w")
      ϕ[irange,jrange,krange] = h5read(filename, "ϕ")
      b[irange,jrange,krange] = h5read(filename, "b")
    end
    imsave(@sprintf("fig/u/%010d.png", n), Array(u[:,1,:]'), origin="lower", vmin=-mu, vmax=mu, cmap="RdBu_r")
    imsave(@sprintf("fig/v/%010d.png", n), Array(v[:,1,:]'), origin="lower", vmin=-mv, vmax=mv, cmap="RdBu_r")
    imsave(@sprintf("fig/w/%010d.png", n), Array(w[:,1,:]'), origin="lower", vmin=-mw, vmax=mw, cmap="RdBu_r")
    imsave(@sprintf("fig/ϕ/%010d.png", n), Array(ϕ[:,1,:]'), origin="lower")
    imsave(@sprintf("fig/b/%010d.png", n), Array(b[:,1,:]'), origin="lower")
    figure(figsize=(9.6, 4.8))
    PyPlot.axes(aspect=1)
    imshow(Array(u[:,1,:]'), vmin=-mu, vmax=mu, origin="lower", cmap="RdBu_r")
    contour(Array((b[:,1,:] + N^2*(x[:,1,:]*sin(θ) + z[:,1,:]*cos(θ)))'), levels=0:1e-8*1000:2e-6*1000, colors="black",
            linewidths=.75)
    savefig(@sprintf("fig/comb/%010d.svg", n), dpi=300)
    close()
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
