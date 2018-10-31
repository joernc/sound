using Printf
using PyPlot
using HDF5

const Δx = 32768/512
#const Δx = 2048/512
const Δy = Δx
const Δz = 1024/256

#const f = -5.5e-5
const f = -5e-5
const N = 0.#1.3e-3
const θ = 0.#2e-3

const c = 2.

const μ = Δz/Δx

#const datadir = "/central/groups/oceanphysics/data_exp"
#const datadir = "/central/groups/oceanphysics/test"
#const datadir = "/central/groups/oceanphysics/noslope"
#const datadir = "/central/groups/oceanphysics/doubleres"
#const datadir = "/central/groups/oceanphysics/small"
#const datadir = "/central/groups/oceanphysics/salmon"
#const datadir = "/central/groups/oceanphysics/highvisc"
const datadir = "/central/groups/oceanphysics/sztztzsz"

pygui(false)
rc("contour", negative_linestyle="solid")

# get index ranges for tiles
tile_range(i, j, k, tile_sizes) = [sum(tile_sizes[1][1:i-1])+1:sum(tile_sizes[1][1:i]),
                                   sum(tile_sizes[2][1:j-1])+1:sum(tile_sizes[2][1:j]),
                                   sum(tile_sizes[3][1:k-1])+1:sum(tile_sizes[3][1:k])]

# assemble and plot results
function assemble(steps, tile_sizes, mu, mv, mw)
  # number of tiles
  ni = length(tile_sizes[1])
  nj = length(tile_sizes[2])
  nk = length(tile_sizes[3])
  # total number of grid points
  nx = sum(tile_sizes[1])
  ny = sum(tile_sizes[2])
  nz = sum(tile_sizes[3])
  # coordinates
  x = [(i-1)*Δx for i =1:nx, j = 1:ny, k = 1:nz]
  z = [(k-1)*Δz for i =1:nx, j = 1:ny, k = 1:nz]
  xp = x*cos(θ) - z*sin(θ)
  zp = x*sin(θ) + z*cos(θ)
  # assemble masks
  maskx = Array{Int8, 3}(undef, nx, ny, nz)
  masky = Array{Int8, 3}(undef, nx, ny, nz)
  maskz = Array{Int8, 3}(undef, nx, ny, nz)
  maski = BitArray{3}(undef, nx, ny, nz)
  for k = 1:nk, j = 1:nj, i = 1:ni
    irange, jrange, krange = tile_range(i, j, k, tile_sizes)
    filename = @sprintf("%s/masks_%1d_%1d_%1d.h5", datadir, i, j, k)
    maskx[irange,jrange,krange] = h5read(filename, "maskx")
    masky[irange,jrange,krange] = h5read(filename, "masky")
    maskz[irange,jrange,krange] = h5read(filename, "maskz")
    maski[irange,jrange,krange] = convert.(Bool, h5read(filename, "maski"))
  end
  # save mask images
  imsave("fig/maski.png", Array(maski[:,:,160]'), origin="lower")
  imsave("fig/maskx.png", Array(maskx[:,:,160]'), origin="lower")
  imsave("fig/masky.png", Array(masky[:,:,160]'), origin="lower")
  imsave("fig/maskz.png", Array(maskz[:,:,160]'), origin="lower")
  # initialize
  u = Array{Float64,3}(undef, nx, ny, nz)
  v = Array{Float64,3}(undef, nx, ny, nz)
  w = Array{Float64,3}(undef, nx, ny, nz)
  ϕ = Array{Float64,3}(undef, nx, ny, nz)
  α = Array{Float64,3}(undef, nx, ny, nz)
  for n in steps
    println(n)
    for k = 1:nk, j = 1:nj, i = 1:ni
      irange, jrange, krange = tile_range(i, j, k, tile_sizes)
      filename = @sprintf("%s/%010d_%1d_%1d_%1d.h5", datadir, n, i, j, k)
      u[irange,jrange,krange] = h5read(filename, "u")
      v[irange,jrange,krange] = h5read(filename, "v")
      w[irange,jrange,krange] = h5read(filename, "w")
      ϕ[irange,jrange,krange] = h5read(filename, "ϕ")
      α[irange,jrange,krange] = h5read(filename, "α")
    end
#    ω = (circshift(v, (-1, 0, 0)) - circshift(v, (1, 0, 0)))/2Δx - (circshift(u, (0, -1, 0)) - circshift(u, (0, 1, 0)))/2Δy
#    imsave(@sprintf("fig/u/%010d.png", n), Array(u[:,:,256]'), origin="lower", vmin=-mu, vmax=mu, cmap="RdBu_r")
#    imsave(@sprintf("fig/v/%010d.png", n), Array(v[:,:,256]'), origin="lower", vmin=-mv, vmax=mv, cmap="RdBu_r")
    imsave(@sprintf("fig/w/%010d.png", n), Array(w[:,1,:]'), origin="lower") #, vmin=-mw, vmax=mw, cmap="RdBu_r")
#    imsave(@sprintf("fig/ϕ/%010d.png", n), Array(ϕ[:,:,256]'), origin="lower")
#    imsave(@sprintf("fig/b/%010d.png", n), Array(b[:,:,256]'), origin="lower")
#    imsave(@sprintf("fig/ω/%010d.png", n), Array(ω[:,:,256]'), origin="lower", vmin=f, vmax=-f, cmap="RdBu_r")
#    b += N^2*(x*sin(θ) + z*cos(θ))
#    figure(figsize=(9.6, 4.8))
#    PyPlot.axes(aspect=1/μ)
#    pcolormesh(Array(xp[:,1,:]'), Array(zp[:,1,:]'), Array(v[:,416,:]'), vmin=-mv, vmax=mv, cmap="RdBu_r", rasterized=true)
#    contour(Array(xp[:,1,:]'), Array(zp[:,1,:]'), Array(b[:,416,:]'), levels=0:2e-8*1024:2e-6*1024, colors="black", linewidths=.75)
#    savefig(@sprintf("fig/iso/%010d.svg", n), dpi=300)
#    close()
    # print conservation diagnostics
#    println(@sprintf("%16.10e", mass(ϕ, maski, maskx, masky, maskz)))
#    println(@sprintf("%16.10e", energy(u, v, w, ϕ, b, z, maski, maskx, masky, maskz)))
#    println(@sprintf("%16.10e", buoyancy(α, maski, maskx, masky, maskz)))
    println(@sprintf("%16.10e", buoyancy_variance(ϕ, α, maski, maskx, masky, maskz)))
  end
end

# mass
function mass(ϕ, maski, maskx, masky, maskz)
  maskb = (maskx .> 1) .| (masky .> 1) .| (maskz .> 1)
  return sum(ϕ[maski]) + sum(ϕ[maskb])/2
end

# energy
function energy(u, v, w, ϕ, b, z, maski, maskx, masky, maskz)
  α = ϕ.*b/c^2
  maskb = (maskx .> 1) .| (masky .> 1) .| (maskz .> 1)
  return sum(u[maski].^2)/2 + sum(v[maski].^2)/2 + sum(w[maski].^2)/2μ^2 + sum(ϕ[maski].^2)/2c^2 + sum(ϕ[maskb].^2)/4c^2 - sum(α[maski].*z[maski]) - sum(α[maskb].*z[maskb])/2
end

# buoyancy
function buoyancy(α, maski, maskx, masky, maskz)
  maskb = (maskx .> 1) .| (masky .> 1) .| (maskz .> 1)
  return c^2*(sum(α[maski]) + sum(α[maskb])/2)
end

# buoyancy variance
function buoyancy_variance(ϕ, α, maski, maskx, masky, maskz)
  b = c^2*α./ϕ
  maskb = (maskx .> 1) .| (masky .> 1) .| (maskz .> 1)
  return sum(ϕ[maski].*b[maski].^2) + sum(ϕ[maskb].*b[maskb].^2)/2
end
