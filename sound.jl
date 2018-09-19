# Implementation of Rick Salmon's "An Ocean Circulation Model Based on Operator-Splitting, Hamiltonian Brackets, and the Inclusion of Sound Waves" (JPO, 2009)

# next steps:
# - implement aspect ratio trick
# - extend to third dimension
# - allow specification of nonzero BC
# - organize channels better?
# - type annotations for better performance?
# - save as HDF5 files

@everywhere using Printf
@everywhere using HDF5

# grid spacing (domain height 1)
@everywhere const h = 1/512

# sound speed
@everywhere const c = 100.

# inertial period
@everywhere const f = 0.

# forcing frequency
@everywhere const ω = 22.4

# initial stratification
@everywhere const N = 160.

# time step
@everywhere const Δt = h/c

# viscosity/diffusion
#@everywhere const ν = 1.5e-6
@everywhere const ν = 1e-2
@everywhere const α = ν/(c*h)

@everywhere function exchangex!(a, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  # send edges
  put!(chan_send_w, a[2,:])
  put!(chan_send_e, a[end-1,:])
  # receive edges
  a[1,:] = take!(chan_receive_w)
  a[end,:] = take!(chan_receive_e)
end

@everywhere function exchangey!(a, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  # send edges
  put!(chan_send_s, a[:,2])
  put!(chan_send_n, a[:,end-1])
  # receive edges
  a[:,1] = take!(chan_receive_s)
  a[:,end] = take!(chan_receive_n)
end

# sound wave split (x-direction)
@everywhere function Sx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  # tile size
  nx, ny = size(u)
  # fields at initial time
  up = copy(u)
  ϕp = copy(ϕ)
  bp = copy(b)
  # loop over grid points
  for i = 2:nx-1
    for j = 1:ny
      if maskx[i,j] == 1 # interior
	u[i,j] = (up[i-1,j] + up[i+1,j])/2 + (ϕp[i-1,j] - ϕp[i+1,j])/2c
	ϕ[i,j] = (ϕp[i-1,j] + ϕp[i+1,j])/2 + c/2*(up[i-1,j] - up[i+1,j])
	b[i,j] = ϕp[i,j]*bp[i,j]/ϕ[i,j]
      elseif maskx[i,j] == 2 # western boundary
	u[i,j] = 0.
	ϕ[i,j] = ϕp[i+1,j] - c*up[i+1,j]
	b[i,j] = ϕp[i,j]*bp[i,j]/ϕ[i,j]
      elseif maskx[i,j] == 3 # eastern boundary
	u[i,j] = 0.
	ϕ[i,j] = ϕp[i-1,j] + c*up[i-1,j]
	b[i,j] = ϕp[i,j]*bp[i,j]/ϕ[i,j]
      end
    end
  end
  # exchange with neighboring tiles
  exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(ϕ, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(b, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

# sound wave split (y-direction)
@everywhere function Sy!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  # tile size
  nx, ny = size(v)
  # fields at initial time
  vp = copy(v)
  ϕp = copy(ϕ)
  bp = copy(b)
  # loop over grid points
  for i = 1:nx
    for j = 2:ny-1
      if masky[i,j] == 1 # interior
	v[i,j] = (vp[i,j-1] + vp[i,j+1])/2 + (ϕp[i,j-1] - ϕp[i,j+1])/2c
	ϕ[i,j] = (ϕp[i,j-1] + ϕp[i,j+1])/2 + c/2*(vp[i,j-1] - vp[i,j+1])
	b[i,j] = ϕp[i,j]*bp[i,j]/ϕ[i,j]
      elseif masky[i,j] == 2 # southern boundary
	v[i,j] = 0.
	ϕ[i,j] = ϕp[i,j+1] - c*vp[i,j+1]
	b[i,j] = ϕp[i,j]*bp[i,j]/ϕ[i,j]
      elseif masky[i,j] == 3 # northern boundary
	v[i,j] = 0.
	ϕ[i,j] = ϕp[i,j-1] + c*vp[i,j-1]
	b[i,j] = ϕp[i,j]*bp[i,j]/ϕ[i,j]
      end
    end
  end
  # exchange with neighboring tiles
  exchangey!(v, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(ϕ, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(b, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
end

# buoyancy split (x-direction, no gravity)
@everywhere function Tx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  # tile size
  nx, ny = size(u)
  # fields at initial time
  up = copy(u)
  bp = copy(b)
  # fields at midpoint
  ui = Array{Float64, 2}(undef, nx, ny)
  bi = Array{Float64, 2}(undef, nx, ny)
  # loop over grid points
  for i = 2:nx-1
    for j = 1:ny
      if maskx[i,j] == 1 # interior
	ui[i,j] = up[i,j]
	bi[i,j] = bp[i,j] + c/8*((bp[i-1,j] + bp[i,j])*(up[i-1,j] + up[i,j]) - (bp[i,j] + bp[i+1,j])*(up[i,j] + up[i+1,j]))/ϕ[i,j]
      elseif (maskx[i,j] == 2) # western boundary
	ui[i,j] = 0.
	bi[i,j] = bp[i,j] + c/4*(-(bp[i,j] + bp[i+1,j])*up[i+1,j])/ϕ[i,j]
      elseif (maskx[i,j] == 3) # eastern boundary
	ui[i,j] = 0.
	bi[i,j] = bp[i,j] + c/4*((bp[i-1,j] + bp[i,j])*up[i-1,j])/ϕ[i,j]
      end
    end
  end
  # exchange with neighboring tiles
  exchangex!(ui, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(bi, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  # loop over grid points
  for i = 2:nx-1
    for j = 1:ny
      if maskx[i,j] == 1 # interior
	u[i,j] = up[i,j]
	b[i,j] = bp[i,j] + c/4*((bi[i-1,j] + bi[i,j])*(ui[i-1,j] + ui[i,j]) - (bi[i,j] + bi[i+1,j])*(ui[i,j] + ui[i+1,j]))/ϕ[i,j]
      elseif (maskx[i,j] == 2) # west boundary
	u[i,j] = 0.
	b[i,j] = bp[i,j] + c/2*(-(bi[i,j] + bi[i+1,j])*ui[i+1,j])/ϕ[i,j]
      elseif (maskx[i,j] == 3) # east boundary
	u[i,j] = 0.
	b[i,j] = bp[i,j] + c/2*((bi[i-1,j] + bi[i,j])*ui[i-1,j])/ϕ[i,j]
      end
    end
  end
  # exchange with neighboring tiles
  exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(b, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

# buoyancy split (y-direction, gravity)
@everywhere function Ty!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  # tile size
  nx, ny = size(v)
  # fields at initial time
  vp = copy(v)
  bp = copy(b)
  # fields at midpoint
  vi = Array{Float64, 2}(undef, nx, ny)
  bi = Array{Float64, 2}(undef, nx, ny)
  # loop over grid points
  for i = 1:nx
    for j = 2:ny-1
      if masky[i,j] == 1 # interior
	vi[i,j] = vp[i,j] + Δt/8*(bp[i,j-1] + 2bp[i,j] + bp[i,j+1])
	bi[i,j] = bp[i,j] + c/8*((bp[i,j-1] + bp[i,j])*(vp[i,j-1] + vp[i,j]) - (bp[i,j] + bp[i,j+1])*(vp[i,j] + vp[i,j+1]))/ϕ[i,j]
      elseif masky[i,j] == 2 # southern boundary
	vi[i,j] = 0.
	bi[i,j] = bp[i,j] + c/4*(-(bp[i,j] + bp[i,j+1])*vp[i,j+1])/ϕ[i,j]
      elseif masky[i,j] == 3 # northern boundary
	vi[i,j] = 0.
	bi[i,j] = bp[i,j] + c/4*((bp[i,j-1] + bp[i,j])*vp[i,j-1])/ϕ[i,j]
      end
    end
  end
  # exchange with neighboring tiles
  exchangey!(vi, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(bi, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  # loop over grid points
  for i = 1:nx
    for j = 2:ny-1
      if masky[i,j] == 1 # interior
	v[i,j] = vp[i,j] + Δt/4*(bi[i,j-1] + 2bi[i,j] + bi[i,j+1])
	b[i,j] = bp[i,j] + c/4*((bi[i,j-1] + bi[i,j])*(vi[i,j-1] + vi[i,j]) - (bi[i,j] + bi[i,j+1])*(vi[i,j] + vi[i,j+1]))/ϕ[i,j]
      elseif masky[i,j] == 2 # southern boundary
	v[i,j] = 0.
	b[i,j] = bp[i,j] + c/2*(-(bi[i,j] + bi[i,j+1])*vi[i,j+1])/ϕ[i,j]
      elseif masky[i,j] == 3 # northern boundary
	v[i,j] = 0.
	b[i,j] = bp[i,j] + c/2*((bi[i,j-1] + bi[i,j])*vi[i,j-1])/ϕ[i,j]
      end
    end
  end
  # exchange with neighboring tiles
  exchangey!(v, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(b, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
end

# rotation split
@everywhere function Rz!(u, v, ϕ, maski, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
			chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)
  # tile size
  nx, ny = size(u)
  # fields at initial time
  up = copy(u)
  vp = copy(v)
  # loop over grid points
  for i = 2:nx-1
    for j = 2:ny-1
      if maski[i,j] # interior
        ωz = f + (vp[i+1,j] - vp[i-1,j] - up[i,j+1] + up[i,j-1]) / 2h
        γ = Δt*ωz*c^2/ϕ[i,j]
        Cz = (1 - γ.^2/4)/(1 + γ^2/4)
        Sz = γ/(1 + γ^2/4)
        u[i,j] = up[i,j]*Cz + vp[i,j]*Sz
	v[i,j] = vp[i,j]*Cz - up[i,j]*Sz
      end
    end
  end
  # exchange with neighboring tiles
  exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(v, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangey!(u, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(v, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
end

# apply forcing
@everywhere function F!(t, u)
  u .+= ω*cos(ω*t)*Δt
end

# x-diffusion with no-flux boundary conditions (should allow prescription of flux)
@everywhere function Dx!(a, maskx, diri, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  # tile size
  nx, ny = size(a)
  # field at initial time
  ap = copy(a)
  # loop over grid points
  for i = 2:nx-1
    for j = 1:ny
      if maskx[i,j] == 1 # interior
	a[i,j] = (1-2α)*ap[i,j] + α*(ap[i-1,j] + ap[i+1,j])
      elseif maskx[i,j] == 2 # west boundary
	if diri[i,j] # Dirichlet BC
	  a[i,j] = 0.
	else # Neumann BC
	  a[i,j] = (1-2α)*ap[i,j] + 2α*ap[i+1,j]
	end
      elseif maskx[i,j] == 3 # east boundary
	if diri[i,j] # Dirichlet BC
	  a[i,j] = 0.
	else # Neumann BC
	  a[i,j] = (1-2α)*ap[i,j] + 2α*ap[i-1,j]
	end
      end
    end
  end
  # exchange with neighboring tiles
  exchangex!(a, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

# y-diffusion with no-flux boundary conditions (should allow prescription of flux)
@everywhere function Dy!(a, masky, diri, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  # tile size
  nx, ny = size(a)
  # field at initial time
  ap = copy(a)
  # loop over grid points
  for i = 1:nx
    for j = 2:ny-1
      if masky[i,j] == 1 # interior
	a[i,j] = (1-2α)*ap[i,j] + α*(ap[i,j-1] + ap[i,j+1])
      elseif masky[i,j] == 2 # south boundary
	if diri[i,j] # Dirichlet BC
	  a[i,j] = 0.
	else # Neumann BC
	  a[i,j] = (1-2α)*ap[i,j] + 2α*ap[i,j+1]
	end
      elseif masky[i,j] == 3 # north boundary
	if diri[i,j] # Dirichlet BC
	  a[i,j] = 0.
	else # Neumann BC
	  a[i,j] = (1-2α)*ap[i,j] + 2α*ap[i,j-1]
	end
      end
    end
  end
  # exchange with neighboring tiles
  exchangey!(a, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
end

# boundary masks (see Fig. 1 in Salmon for a sketch)
@everywhere function boundary_masks(fluid, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
				    chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)
  nx, ny = size(fluid)
  maskx = zeros(UInt8, nx, ny)
  masky = zeros(UInt8, nx, ny)
  for i = 2:nx-1
    for j = 2:ny-1
      # east–west boundaries
      if fluid[i-1,j-1] & fluid[i,j-1] & fluid[i-1,j] & fluid[i,j] # interior
	maskx[i,j] = 1
      elseif (!fluid[i-1,j-1] | !fluid[i-1,j]) & fluid[i,j-1] & fluid[i,j] # western boundary
	maskx[i,j] = 2
      elseif (!fluid[i,j-1] | !fluid[i,j]) & fluid[i-1,j-1] & fluid[i-1,j] # eastern boundary
	maskx[i,j] = 3
      end
      # north–south boundaries
      if fluid[i-1,j-1] & fluid[i-1,j] & fluid[i,j-1] & fluid[i,j]
	masky[i,j] = 1 # interior
      elseif (!fluid[i-1,j-1] | !fluid[i,j-1]) & fluid[i-1,j] & fluid[i,j]
	masky[i,j] = 2 # southern boundary
      elseif (!fluid[i-1,j] | !fluid[i,j]) & fluid[i-1,j-1] & fluid[i,j-1]
	masky[i,j] = 3 # northern boundary
      end
    end
  end
  # interior points
  maski = (maskx .== 1) .& (masky .== 1)
  # exchange edges
  exchangex!(maski, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangey!(maski, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangex!(maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangey!(maskx, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangex!(masky, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangey!(masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  return maski, maskx, masky
end

# read in fluid mask
@everywhere function read_topo(irange, jrange, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
			       chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)
  # tile size plus two points padding
  nx = length(irange) + 2
  ny = length(jrange) + 2
  # read from file
  full_mask = trues(1024, 512)
  full_mask[256:767,128:383] = h5read("julia_256_512.h5", "m") .> .5
  # add solid boundary at top
  full_mask[:,end] .= false
  # assign to tile
  fluid = Array{Bool}(undef, nx, ny)
  fluid[2:nx-1,2:ny-1] = full_mask[irange,jrange]
  # exchange edges
  exchangex!(fluid, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangey!(fluid, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  # interior, x-, and y-masks
  maski, maskx, masky = boundary_masks(fluid, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
				       chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)
  return maski, maskx, masky
end

# run model on tile
@everywhere function run_tile(i, j, irange, jrange, steps, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
			      chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)
  # read topography mask
  maski, maskx, masky = read_topo(irange, jrange, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
				  chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)
  # initialization
  nx = length(irange) + 2
  ny = length(jrange) + 2
  u = zeros(nx, ny)
  v = zeros(nx, ny)
  y = [(j-1)*h for i = irange[1]-1:irange[end]+1, j = jrange[1]-1:jrange[end]+1]
  ϕ = c^2*ones(nx, ny) + N^2*y.^2/2
  b = N^2*y
  # specify where Dirichlet boundary conditions are to be imposed
  diriu = (maskx .> 1) .| (masky .> 1)
  diriv = (maskx .> 1) .| (masky .> 1)
  dirib = falses(nx, ny)
  # time steps
  for k = 1:steps
    if k % 100 == 1
      println(@sprintf("%6i %9.3e %9.3e", k-1, (k-1)*2Δt, maximum(hypot.(u[2:nx-1,2:ny-1], v[2:nx-1,2:ny-1]))))
      # discard edges
      us = u[2:nx-1,2:ny-1]
      vs = v[2:nx-1,2:ny-1]
      ϕs = ϕ[2:nx-1,2:ny-1]
      bs = b[2:nx-1,2:ny-1]
      # calculate vorticity
      ωz = f .+ (v[3:nx,2:ny-1] - v[1:nx-2,2:ny-1] - u[2:nx-1,3:ny] + u[2:nx-1,1:ny-2])/2h
      # replace missing values with NaNs
      us[((maskx .== 0) .& (masky .== 0))[2:nx-1, 2:ny-1]] .= NaN
      vs[((maskx .== 0) .& (masky .== 0))[2:nx-1, 2:ny-1]] .= NaN
      ϕs[((maskx .== 0) .& (masky .== 0))[2:nx-1, 2:ny-1]] .= NaN
      bs[((maskx .== 0) .& (masky .== 0))[2:nx-1, 2:ny-1]] .= NaN
      ωz[.!maski[2:nx-1,2:ny-1]] .= NaN
      # save data
      open(@sprintf("data/u/%1d_%1d_%010d", i, j, k-1), "w") do file
	write(file, us)
      end
      open(@sprintf("data/v/%1d_%1d_%010d", i, j, k-1), "w") do file
	write(file, vs)
      end
      open(@sprintf("data/ϕ/%1d_%1d_%010d", i, j, k-1), "w") do file
	write(file, ϕs)
      end
      open(@sprintf("data/b/%1d_%1d_%010d", i, j, k-1), "w") do file
	write(file, bs)
      end
      open(@sprintf("data/ωz/%1d_%1d_%010d", i, j, k-1), "w") do file
	write(file, ωz)
      end
    end
    # Strang splitting
    Dx!(u, maskx, diriu, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(v, maskx, diriv, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(b, maskx, dirib, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dy!(u, masky, diriu, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(v, masky, diriv, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(b, masky, dirib, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    F!((2k-1.5)*Δt, u)
    Rz!(u, v, ϕ, maski, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
	chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)
    Tx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Sx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Sy!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Ty!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Ty!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Sy!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Sx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Tx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Rz!(u, v, ϕ, maski, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
        chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)
    F!((2k-.5)*Δt, u)
    Dy!(u, masky, diriu, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(v, masky, diriv, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(b, masky, dirib, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dx!(u, maskx, diriu, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(v, maskx, diriv, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(b, maskx, dirib, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  end
  return
end

# get index ranges for tiles
tile_range(i, j, tile_sizes) = [sum(tile_sizes[1][1:i-1])+1:sum(tile_sizes[1][1:i]),
				sum(tile_sizes[2][1:j-1])+1:sum(tile_sizes[2][1:j])]

# run the model
function run_model(steps, tile_sizes)
  # number of tiles
  ni = length(tile_sizes[1])
  nj = length(tile_sizes[2])
  # set up remote channels for edge communication
  chan_w = [RemoteChannel(()->Channel{Array{Float64, 1}}(1)) for i = 1:ni, j = 1:nj]
  chan_e = [RemoteChannel(()->Channel{Array{Float64, 1}}(1)) for i = 1:ni, j = 1:nj]
  chan_s = [RemoteChannel(()->Channel{Array{Float64, 1}}(1)) for i = 1:ni, j = 1:nj]
  chan_n = [RemoteChannel(()->Channel{Array{Float64, 1}}(1)) for i = 1:ni, j = 1:nj]
  # spawn work on tiles
  a = Array{Future}(undef, ni, nj)
  for i in 1:ni
    for j in 1:nj
      irange, jrange = tile_range(i, j, tile_sizes)
      p = workers()[mod1((i-1)*nj+j, nprocs()-1)]
      println(p, ": ", i, " ", j, " ", irange, " ", jrange)
      a[i,j] = @spawnat p run_tile(i, j, irange, jrange, steps, chan_w[i,j], chan_e[i,j], chan_s[i,j], chan_n[i,j],
				   chan_e[mod1(i-1,ni),j], chan_w[mod1(i+1,ni),j], chan_n[i,mod1(j-1,nj)], chan_s[i,mod1(j+1,nj)])
    end
  end
  # wait for results
  wait.(a)
end
