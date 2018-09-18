# Implementation of Rick Salmon's "An Ocean Circulation Model Based on Operator-Splitting, Hamiltonian Brackets, and the Inclusion of Sound Waves" (JPO, 2009)

# next steps:
# - extend to third dimension
# - allow specification of nonzero BC
# - organize channels better?
# - type annotations for better performance?

@everywhere using Printf
@everywhere using HDF5

# grid spacing (domain height 1)
@everywhere const h = 1/256

# sound speed
@everywhere const c = 1.

# inertial period
@everywhere const f = 0.

# forcing frequency
@everywhere const ω = .01

# forcing amplitude
@everywhere const u0 = .01

# time step
@everywhere const Δt = h/c

# viscosity/diffusion
@everywhere const ν = 1e-5
@everywhere const γ = ν/(c*h)

@everywhere function exchangex!(a, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  # send edges
  put!(channel_send_west, a[2,:])
  put!(channel_send_east, a[end-1,:])
  # receive edges
  a[1,:] = take!(channel_receive_west)
  a[end,:] = take!(channel_receive_east)
  return
end

@everywhere function exchangey!(a, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  # send edges
  put!(channel_send_south, a[:,2])
  put!(channel_send_north, a[:,end-1])
  # receive edges
  a[:,1] = take!(channel_receive_south)
  a[:,end] = take!(channel_receive_north)
  return
end

# sound wave split (x-direction)
@everywhere function Sx(u, ϕ, b, maskx, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  nx, ny = size(u)
  up = Array{Float64, 2}(undef, nx, ny)
  ϕp = Array{Float64, 2}(undef, nx, ny)
  bp = Array{Float64, 2}(undef, nx, ny)
  for i = 2:nx-1
    for j = 1:ny
      if maskx[i,j] == 1 # interior
	up[i,j] = (u[i-1,j] + u[i+1,j])/2 + (ϕ[i-1,j] - ϕ[i+1,j])/2c
	ϕp[i,j] = (ϕ[i-1,j] + ϕ[i+1,j])/2 + c/2*(u[i-1,j] - u[i+1,j])
	bp[i,j] = ϕ[i,j]*b[i,j]/ϕp[i,j]
      elseif maskx[i,j] == 2 # western boundary
	up[i,j] = 0.
	ϕp[i,j] = ϕ[i+1,j] - c*u[i+1,j]
	bp[i,j] = ϕ[i,j]*b[i,j]/ϕp[i,j]
      elseif maskx[i,j] == 3 # eastern boundary
	up[i,j] = 0.
	ϕp[i,j] = ϕ[i-1,j] + c*u[i-1,j]
	bp[i,j] = ϕ[i,j]*b[i,j]/ϕp[i,j]
      else
	up[i,j] = u[i,j]
	ϕp[i,j] = ϕ[i,j]
	bp[i,j] = b[i,j]
      end
    end
  end
  exchangex!(up, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangex!(ϕp, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangex!(bp, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  return up, ϕp, bp
end

# sound wave split (y-direction)
@everywhere function Sy(v, ϕ, b, masky, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  nx, ny = size(v)
  vp = Array{Float64, 2}(undef, nx, ny)
  ϕp = Array{Float64, 2}(undef, nx, ny)
  bp = Array{Float64, 2}(undef, nx, ny)
  for i = 1:nx
    for j = 2:ny-1
      if masky[i,j] == 1 # interior
	vp[i,j] = (v[i,j-1] + v[i,j+1])/2 + (ϕ[i,j-1] - ϕ[i,j+1])/2c
	ϕp[i,j] = (ϕ[i,j-1] + ϕ[i,j+1])/2 + c/2*(v[i,j-1] - v[i,j+1])
	bp[i,j] = ϕ[i,j]*b[i,j]/ϕp[i,j]
      elseif masky[i,j] == 2 # south boundary
	vp[i,j] = 0.
	ϕp[i,j] = ϕ[i,j+1] - c*v[i,j+1]
	bp[i,j] = ϕ[i,j]*b[i,j]/ϕp[i,j]
      elseif masky[i,j] == 3 # north boundary
	vp[i,j] = 0.
	ϕp[i,j] = ϕ[i,j-1] + c*v[i,j-1]
	bp[i,j] = ϕ[i,j]*b[i,j]/ϕp[i,j]
      else
	vp[i,j] = v[i,j]
	ϕp[i,j] = ϕ[i,j]
	bp[i,j] = b[i,j]
      end
    end
  end
  exchangey!(vp, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  exchangey!(ϕp, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  exchangey!(bp, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  return vp, ϕp, bp
end

# buoyancy split (x-direction, no gravity)
@everywhere function Tx(u, ϕ, b, maskx, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  nx, ny = size(u)
  ui = Array{Float64, 2}(undef, nx, ny)
  bi = Array{Float64, 2}(undef, nx, ny)
  for i = 2:nx-1
    for j = 1:ny
      if maskx[i,j] == 1 # interior
	ui[i,j] = u[i,j]
	α = ϕ[i,j]*b[i,j]/c^2 + ((b[i-1,j] + b[i,j])*(u[i-1,j] + u[i,j]) - (b[i,j] + b[i+1,j])*(u[i,j] + u[i+1,j]))/8c
	bi[i,j] = c^2*α/ϕ[i,j]
      elseif (maskx[i,j] == 2) # western boundary
	ui[i,j] = 0.
	α = ϕ[i,j]*b[i,j]/c^2 + (-(b[i,j] + b[i+1,j])*u[i+1,j])/4c
	bi[i,j] = c^2*α/ϕ[i,j]
      elseif (maskx[i,j] == 3) # eastern boundary
	ui[i,j] = 0.
	α = ϕ[i,j]*b[i,j]/c^2 + ((b[i-1,j] + b[i,j])*u[i-1,j])/4c
	bi[i,j] = c^2*α/ϕ[i,j]
      else
	ui[i,j] = u[i,j]
	bi[i,j] = b[i,j]
      end
    end
  end
  exchangex!(ui, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangex!(bi, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  up = Array{Float64, 2}(undef, nx, ny)
  bp = Array{Float64, 2}(undef, nx, ny)
  for i = 2:nx-1
    for j = 1:ny
      if maskx[i,j] == 1 # interior
	up[i,j] = u[i,j]
	α = ϕ[i,j]*b[i,j]/c^2 + ((bi[i-1,j] + bi[i,j])*(ui[i-1,j] + ui[i,j]) - (bi[i,j] + bi[i+1,j])*(ui[i,j] + ui[i+1,j]))/4c
	bp[i,j] = c^2*α/ϕ[i,j]
      elseif (maskx[i,j] == 2) # west boundary
	up[i,j] = 0.
	α = ϕ[i,j]*b[i,j]/c^2 + (-(bi[i,j] + bi[i+1,j])*ui[i+1,j])/2c
	bp[i,j] = c^2*α/ϕ[i,j]
      elseif (maskx[i,j] == 3) # east boundary
	up[i,j] = 0.
	α = ϕ[i,j]*b[i,j]/c^2 + ((bi[i-1,j] + bi[i,j])*ui[i-1,j])/2c
	bp[i,j] = c^2*α/ϕ[i,j]
      else
	up[i,j] = ui[i,j]
	bp[i,j] = bi[i,j]
      end
    end
  end
  exchangex!(up, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangex!(bp, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  return up, bp
end

# buoyancy split (y-direction, gravity)
@everywhere function Ty(v, ϕ, b, masky, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  nx, ny = size(v)
  vi = Array{Float64, 2}(undef, nx, ny)
  bi = Array{Float64, 2}(undef, nx, ny)
  for i = 1:nx
    for j = 2:ny-1
      if masky[i,j] == 1 # interior
	vi[i,j] = v[i,j] + Δt/8*(b[i,j-1] + 2b[i,j] + b[i,j+1])
	α = ϕ[i,j]*b[i,j]/c^2 + ((b[i,j-1] + b[i,j])*(v[i,j-1] + v[i,j]) - (b[i,j] + b[i,j+1])*(v[i,j] + v[i,j+1]))/8c
	bi[i,j] = c^2*α/ϕ[i,j]
      elseif masky[i,j] == 2 # southern boundary
	vi[i,j] = 0.
	α = ϕ[i,j]*b[i,j]/c^2 + (-(b[i,j] + b[i,j+1])*v[i,j+1])/4c
	bi[i,j] = c^2*α/ϕ[i,j]
      elseif masky[i,j] == 3 # northern boundary
	vi[i,j] = 0.
	α = ϕ[i,j]*b[i,j]/c^2 + ((b[i,j-1] + b[i,j])*v[i,j-1])/4c
	bi[i,j] = c^2*α/ϕ[i,j]
      else
	vi[i,j] = v[i,j]
	bi[i,j] = b[i,j]
      end
    end
  end
  exchangey!(vi, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  exchangey!(bi, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  vp = Array{Float64, 2}(undef, nx, ny)
  bp = Array{Float64, 2}(undef, nx, ny)
  for i = 1:nx
    for j = 2:ny-1
      if masky[i,j] == 1 # interior
	vp[i,j] = v[i,j] + Δt/4*(bi[i,j-1] + 2bi[i,j] + bi[i,j+1])
	α = ϕ[i,j]*b[i,j]/c^2 + ((bi[i,j-1] + bi[i,j])*(vi[i,j-1] + vi[i,j]) - (bi[i,j] + bi[i,j+1])*(vi[i,j] + vi[i,j+1]))/4c
	bp[i,j] = c^2*α/ϕ[i,j]
      elseif masky[i,j] == 2 # southern boundary
	vp[i,j] = 0.
	α = ϕ[i,j]*b[i,j]/c^2 + (-(bi[i,j] + bi[i,j+1])*vi[i,j+1])/2c
	bp[i,j] = c^2*α/ϕ[i,j]
      elseif masky[i,j] == 3 # northern boundary
	vp[i,j] = 0.
	α = ϕ[i,j]*b[i,j]/c^2 + ((bi[i,j-1] + bi[i,j])*vi[i,j-1])/2c
	bp[i,j] = c^2*α/ϕ[i,j]
      else
	vp[i,j] = v[i,j]
	bp[i,j] = b[i,j]
      end
    end
  end
  exchangey!(vp, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  exchangey!(bp, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  return vp, bp
end

# rotation split
@everywhere function Rz(u, v, ϕ, maski, channel_send_west, channel_send_east, channel_send_south, channel_send_north,
			channel_receive_west, channel_receive_east, channel_receive_south, channel_receive_north)
  nx, ny = size(u)
  up = Array{Float64, 2}(undef, nx, ny)
  vp = Array{Float64, 2}(undef, nx, ny)
  for i = 2:nx-1
    for j = 2:ny-1
      if maski[i,j] # interior
        ωz = f + (v[i+1,j] - v[i-1,j] - u[i,j+1] + u[i,j-1]) / 2h
        γ = Δt*ωz*c^2/ϕ[i,j]
        Cz = (1 - γ.^2/4)/(1 + γ^2/4)
        Sz = γ/(1 + γ^2/4)
        up[i,j] = u[i,j]*Cz + v[i,j]*Sz
	vp[i,j] = v[i,j]*Cz - u[i,j]*Sz
      else # boundary or solid
	up[i,j] = u[i,j]
	vp[i,j] = v[i,j]
      end
    end
  end
  exchangex!(up, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangex!(vp, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangey!(up, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  exchangey!(vp, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  return up, vp
end

# apply forcing
@everywhere function forcing!(t, u, maski)
  u[maski] .+= u0*ω*cos(ω*t)*Δt
end

# x-diffusion with no-flux boundary conditions (should allow prescription of flux)
@everywhere function Dx(a, maskx, diri, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  nx, ny = size(a)
  ap = Array{Float64, 2}(undef, nx, ny)
  for i = 2:nx-1
    for j = 1:ny
      if maskx[i,j] == 1 # interior
	ap[i,j] = (1-2γ)*a[i,j] + γ*(a[i-1,j] + a[i+1,j])
      elseif maskx[i,j] == 2 # west boundary
	if diri[i,j]
	  ap[i,j] = 0.
	else
	  ap[i,j] = (1-2γ)*a[i,j] + 2γ*a[i+1,j]
	end
      elseif maskx[i,j] == 3 # east boundary
	if diri[i,j]
	  ap[i,j] = 0.
	else
	  ap[i,j] = (1-2γ)*a[i,j] + 2γ*a[i-1,j]
	end
      else
	ap[i,j] = a[i,j]
      end
    end
  end
  exchangex!(ap, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  return ap
end

# y-diffusion with no-flux boundary conditions (should allow prescription of flux)
@everywhere function Dy(a, masky, diri, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  nx, ny = size(a)
  ap = Array{Float64, 2}(undef, nx, ny)
  for i = 1:nx
    for j = 2:ny-1
      if masky[i,j] == 1 # interior
	ap[i,j] = (1-2γ)*a[i,j] + γ*(a[i,j-1] + a[i,j+1])
      elseif masky[i,j] == 2 # south boundary
	if diri[i,j]
	  ap[i,j] = 0.
	else
	  ap[i,j] = (1-2γ)*a[i,j] + 2γ*a[i,j+1]
	end
      elseif masky[i,j] == 3 # north boundary
	if diri[i,j]
	  ap[i,j] = 0.
	else
	  ap[i,j] = (1-2γ)*a[i,j] + 2γ*a[i,j-1]
	end
      else
	ap[i,j] = a[i,j]
      end
    end
  end
  exchangey!(ap, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  return ap
end

# boundary masks (see Fig. 1 in Salmon for a sketch)
@everywhere function boundary_masks(fluid, channel_send_west, channel_send_east, channel_send_south, channel_send_north,
				    channel_receive_west, channel_receive_east, channel_receive_south, channel_receive_north)
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
  exchangex!(maski, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangey!(maski, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  exchangex!(maskx, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangey!(maskx, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  exchangex!(masky, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangey!(masky, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  return maski, maskx, masky
end

# read in fluid mask
@everywhere function read_topo(irange, jrange, channel_send_west, channel_send_east, channel_send_south, channel_send_north,
			       channel_receive_west, channel_receive_east, channel_receive_south, channel_receive_north)
  # tile size plus two points padding
  nx = length(irange) + 2
  ny = length(jrange) + 2
  # read from file
  full_mask = h5read("julia_256_512.h5", "m") .> .5
  # add solid boundary at top
  full_mask[:,end] .= false
  # assign to tile
  fluid = Array{Bool}(undef, nx, ny)
  fluid[2:nx-1,2:ny-1] = full_mask[irange,jrange]
  # exchange edges
  exchangex!(fluid, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangey!(fluid, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  # interior, x-, and y-masks
  maski, maskx, masky = boundary_masks(fluid, channel_send_west, channel_send_east, channel_send_south, channel_send_north,
				       channel_receive_west, channel_receive_east, channel_receive_south, channel_receive_north)
  return maski, maskx, masky
end

# run model on tile
@everywhere function run_tile(i, j, irange, jrange, steps, channel_send_west, channel_send_east, channel_send_south,
			      channel_send_north, channel_receive_west, channel_receive_east, channel_receive_south,
			      channel_receive_north)
  # read topography mask
  maski, maskx, masky = read_topo(irange, jrange, channel_send_west, channel_send_east, channel_send_south, channel_send_north,
				  channel_receive_west, channel_receive_east, channel_receive_south, channel_receive_north)
  # initialization
  nx = length(irange) + 2
  ny = length(jrange) + 2
  u = zeros(nx, ny)
  v = zeros(nx, ny)
  y = [(j-1)*h for i = irange[1]-1:irange[end]+1, j = jrange[1]-1:jrange[end]+1]
  N = .02
  ϕ = c^2*ones(nx, ny) + N^2*y.^2/2
  b = N^2*y
  # specify where Dirichlet boundary conditions are to be imposed
  diriu = (maskx .> 1) .| (masky .> 1)
  diriv = (maskx .> 1) .| (masky .> 1)
  dirib = falses(nx, ny)
  # time steps
  for k = 1:steps
    if k % 100 == 1
      println(@sprintf("%6i %9.3e %9.3e", k-1, (k-1)*Δt, maximum(hypot.(u[2:nx-1,2:ny-1], v[2:nx-1,2:ny-1]))))
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
      open(@sprintf("sound/data/u/%1d_%1d_%010d", i, j, k-1), "w") do file
	write(file, us)
      end
      open(@sprintf("sound/data/v/%1d_%1d_%010d", i, j, k-1), "w") do file
	write(file, vs)
      end
      open(@sprintf("sound/data/ϕ/%1d_%1d_%010d", i, j, k-1), "w") do file
	write(file, ϕs)
      end
      open(@sprintf("sound/data/b/%1d_%1d_%010d", i, j, k-1), "w") do file
	write(file, bs)
      end
      open(@sprintf("sound/data/ωz/%1d_%1d_%010d", i, j, k-1), "w") do file
	write(file, ωz)
      end
    end
    # Dx steps
    u = Dx(u, maskx, diriu, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    v = Dx(v, maskx, diriv, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    b = Dx(b, maskx, dirib, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    # Dy steps
    u = Dy(u, masky, diriu, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    v = Dy(v, masky, diriv, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    b = Dy(b, masky, dirib, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    # forcing
    forcing!((k+.5)*Δt, u, maski)
    # Rz step
    u, v = Rz(u, v, ϕ, maski, channel_send_west, channel_send_east, channel_send_south, channel_send_north, channel_receive_west,
        channel_receive_east, channel_receive_south, channel_receive_north)
    # Tx split
    u, b = Tx(u, ϕ, b, maskx, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    # Sx step
    u, ϕ, b = Sx(u, ϕ, b, maskx, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    # Sy step
    v, ϕ, b = Sy(v, ϕ, b, masky, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    # Ty split
    v, b = Ty(v, ϕ, b, masky, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    # Ty split
    v, b = Ty(v, ϕ, b, masky, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    # Sy step
    v, ϕ, b = Sy(v, ϕ, b, masky, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    # Sx step
    u, ϕ, b = Sx(u, ϕ, b, maskx, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    # Tx split
    u, b = Tx(u, ϕ, b, maskx, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    # Rz step
    u, v = Rz(u, v, ϕ, maski, channel_send_west, channel_send_east, channel_send_south, channel_send_north, channel_receive_west,
	      channel_receive_east, channel_receive_south, channel_receive_north)
    # forcing
    forcing!((k+1.5)*Δt, u, maski)
    # Dy steps
    u = Dy(u, masky, diriu, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    v = Dy(v, masky, diriv, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    b = Dy(b, masky, dirib, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    # Dx steps
    u = Dx(u, maskx, diriu, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    v = Dx(v, maskx, diriv, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    b = Dx(b, maskx, dirib, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
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
  channels_west = [RemoteChannel(()->Channel{Array{Union{Missing, Float64}}}(1)) for i = 1:ni, j = 1:nj]
  channels_east = [RemoteChannel(()->Channel{Array{Union{Missing, Float64}}}(1)) for i = 1:ni, j = 1:nj]
  channels_south = [RemoteChannel(()->Channel{Array{Union{Missing, Float64}}}(1)) for i = 1:ni, j = 1:nj]
  channels_north = [RemoteChannel(()->Channel{Array{Union{Missing, Float64}}}(1)) for i = 1:ni, j = 1:nj]
  # spawn work on tiles
  a = Array{Future}(undef, ni, nj)
  for i in 1:ni
    for j in 1:nj
      irange, jrange = tile_range(i, j, tile_sizes)
      p = workers()[mod1((i-1)*nj+j, nprocs()-1)]
      println(p, ": ", i, " ", j, " ", irange, " ", jrange)
      a[i,j] = @spawnat p run_tile(i, j, irange, jrange, steps, channels_west[i,j], channels_east[i,j], channels_south[i,j],
				   channels_north[i,j], channels_east[mod1(i-1,ni),j], channels_west[mod1(i+1,ni),j],
				   channels_north[i,mod1(j-1,nj)], channels_south[i,mod1(j+1,nj)])
    end
  end
  # wait for result
  for i in 1:ni
    for j in 1:nj
      wait(a[i,j])
    end
  end
end
