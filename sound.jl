# Implementation of Rick Salmon's "An Ocean Circulation Model Based on Operator-Splitting, Hamiltonian Brackets, and the Inclusion of Sound Waves" (JPO, 2009)

# next steps:
# - extend to third dimension
# - organize channels better
# - type annotations for better performance?
# - allow specification of flux BC

@everywhere using Printf
@everywhere using HDF5

# grid spacing (domain width 2π)
@everywhere const h = 2π/1024

# sound speed
@everywhere const c = 1.

# inertial period
@everywhere const f = 0.

# time step
@everywhere const Δt = h/c

# viscosity
@everywhere const ν = 2e-5

# parameter for diffusion
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

#"""Sound wave split (x-direction)"""
@everywhere function Sx(u, ϕ, maskx, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  nx, ny = size(u)
  up = Array{Float64, 2}(undef, nx, ny)
  ϕp = Array{Float64, 2}(undef, nx, ny)
  for i = 2:nx-1
    for j = 1:ny
      if maskx[i,j] == 1 # interior
	up[i,j] = (u[i-1,j] + u[i+1,j])/2 + (ϕ[i-1,j] - ϕ[i+1,j])/2c
	ϕp[i,j] = (ϕ[i-1,j] + ϕ[i+1,j])/2 + c/2*(u[i-1,j] - u[i+1,j])
      elseif maskx[i,j] == 2 # western boundary
	up[i,j] = 0.
	ϕp[i,j] = ϕ[i+1,j] - c*u[i+1,j]
      elseif maskx[i,j] == 3 # eastern boundary
	up[i,j] = 0.
	ϕp[i,j] = ϕ[i-1,j] + c*u[i-1,j]
      else
	up[i,j] = u[i,j]
	ϕp[i,j] = ϕ[i,j]
      end
    end
  end
  exchangex!(up, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangex!(ϕp, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  return up, ϕp
end

#"""Sound wave split (y-direction)"""
@everywhere function Sy(v, ϕ, masky, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  nx, ny = size(v)
  vp = Array{Float64, 2}(undef, nx, ny)
  ϕp = Array{Float64, 2}(undef, nx, ny)
  for i = 1:nx
    for j = 2:ny-1
      if masky[i,j] == 1 # interior
	vp[i,j] = (v[i,j-1] + v[i,j+1])/2 + (ϕ[i,j-1] - ϕ[i,j+1])/2c
	ϕp[i,j] = (ϕ[i,j-1] + ϕ[i,j+1])/2 + c/2*(v[i,j-1] - v[i,j+1])
      elseif masky[i,j] == 2 # south boundary
	vp[i,j] = 0.
	ϕp[i,j] = ϕ[i,j+1] - c*v[i,j+1]
      elseif masky[i,j] == 3 # north boundary
	vp[i,j] = 0.
	ϕp[i,j] = ϕ[i,j-1] + c*v[i,j-1]
      else
	vp[i,j] = v[i,j]
	ϕp[i,j] = ϕ[i,j]
      end
    end
  end
  exchangey!(vp, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  exchangey!(ϕp, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  return vp, ϕp
end

@everywhere function Tx(u, ϕ, b, maskx, diri, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  nx, ny = size(u)
  up = Array{Float64, 2}(undef, nx, ny)
  bp = Array{Float64, 2}(undef, nx, ny)
  for i = 2:nx-1
    for j = 1:ny
      if maskx[i,j] == 1 # interior
	up[i,j] = u[i,j]
	α = ϕ[i,j]*b[i,j]/c^2 + Δt/4h*((b[i-1,j] + b[i,j])*(u[i-1,j] + u[i,j]) - (b[i,j] + b[i+1,j])*(u[i,j] + u[i+1,j]))
	bp[i,j] = c^2*α/ϕ[i,j]
      elseif (maskx[i,j] == 2) & !diri[i,j] # west boundary
	up[i,j] = 0.
	α = ϕ[i,j]*b[i,j]/c^2 + Δt/2h*(-(b[i,j] + b[i+1,j])*u[i+1,j])
	bp[i,j] = c^2*α/ϕ[i,j]
      elseif (maskx[i,j] == 3) & !diri[i,j] # east boundary
	up[i,j] = 0.
	α = ϕ[i,j]*b[i,j]/c^2 + Δt/2h*((b[i-1,j] + b[i,j])*u[i-1,j])
	bp[i,j] = c^2*α/ϕ[i,j]
      else
	up[i,j] = u[i,j]
	bp[i,j] = b[i,j]
      end
    end
  end
  exchangex!(up, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangex!(bp, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  return up, bp
end

@everywhere function Ty(v, ϕ, b, masky, diri, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  nx, ny = size(v)
  vp = Array{Float64, 2}(undef, nx, ny)
  bp = Array{Float64, 2}(undef, nx, ny)
  for i = 1:nx
    for j = 2:ny-1
      if masky[i,j] == 1 # interior
	vp[i,j] = v[i,j] + Δt/4*(b[i,j-1] + 2b[i,j] + b[i,j+1])
	α = ϕ[i,j]*b[i,j]/c^2 + Δt/4h*((b[i,j-1] + b[i,j])*(v[i,j-1] + v[i,j]) - (b[i,j] + b[i,j+1])*(v[i,j] + v[i,j+1]))
	bp[i,j] = c^2*α/ϕ[i,j]
      elseif (masky[i,j] == 2) & !diri[i,j] # south boundary
	vp[i,j] = 0.
	α = ϕ[i,j]*b[i,j]/c^2 + Δt/2h*(-(b[i,j] + b[i,j+1])*v[i,j+1])
	bp[i,j] = c^2*α/ϕ[i,j]
      elseif (masky[i,j] == 3) & !diri[i,j] # north boundary
	vp[i,j] = 0.
	α = ϕ[i,j]*b[i,j]/c^2 + Δt/2h*((b[i,j-1] + b[i,j])*v[i,j-1])
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

#"""Rotation split"""
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

#"""Apply forcing"""
@everywhere function forcing!(t, u, maski)
  u[maski] .+= .01*.01*cos(.01t)*Δt
end

# x-diffusion with no-flux boundary conditions (should allow prescription of flux)
@everywhere function Dx(a, maskx, diri, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  nx, ny = size(a)
  ap = Array{Float64, 2}(undef, nx, ny)
  for i = 2:nx-1
    for j = 1:ny
      if maskx[i,j] == 1 # interior
	ap[i,j] = (1-2γ)*a[i,j] + γ*(a[i-1,j] + a[i+1,j])
      elseif (maskx[i,j] == 2) & !diri[i,j] # west boundary
	ap[i,j] = (1-2γ)*a[i,j] + 2γ*a[i+1,j]
      elseif (maskx[i,j] == 3) & !diri[i,j] # east boundary
	ap[i,j] = (1-2γ)*a[i,j] + 2γ*a[i-1,j]
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
      elseif (masky[i,j] == 2) & !diri[i,j] # south boundary
	ap[i,j] = (1-2γ)*a[i,j] + 2γ*a[i,j+1]
      elseif (masky[i,j] == 3) & !diri[i,j] # north boundary
	ap[i,j] = (1-2γ)*a[i,j] + 2γ*a[i,j-1]
      else
	ap[i,j] = a[i,j]
      end
    end
  end
  exchangey!(ap, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  return ap
end

#"""Boundary masks (see Fig. 1 in Salmon for a sketch)"""
@everywhere function boundary_masks(fluid, channel_send_west, channel_send_east, channel_send_south, channel_send_north,
				    channel_receive_west, channel_receive_east, channel_receive_south, channel_receive_north)
  nx, ny = size(fluid)
  maskx = zeros(UInt8, nx, ny)
  masky = zeros(UInt8, nx, ny)
  for i = 2:nx-1
    for j = 2:ny-1
      # east–west boundaries
      if fluid[i-1,j-1] & fluid[i,j-1] & fluid[i-1,j] & fluid[i,j]
	maskx[i,j] = 1 # interior
      elseif (!fluid[i-1,j-1] | !fluid[i-1,j]) & fluid[i,j-1] & fluid[i,j]
	maskx[i,j] = 2 # western boundary
      elseif (!fluid[i,j-1] | !fluid[i,j]) & fluid[i-1,j-1] & fluid[i-1,j]
	maskx[i,j] = 3 # eastern boundary
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

@everywhere function read_topo(irange, jrange, channel_send_west, channel_send_east, channel_send_south, channel_send_north,
			       channel_receive_west, channel_receive_east, channel_receive_south, channel_receive_north)
  # tile size plus two points padding
  nx = length(irange) + 2
  ny = length(jrange) + 2
  # read from file
  full_mask = h5read("julia_512_1024.h5", "m") .> .5
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
  ϕ = c^2*ones(nx, ny)
  b = .0005*randn(nx, ny)
  b[:,2] .= .005/2π
  b[:,end-1] .= -.005/2π
  # specify where Dirichlet boundary conditions are to be imposed
  diriu = (maskx .> 1) .| (masky .> 1)
  diriv = (maskx .> 1) .| (masky .> 1)
  dirib = falses(nx, ny)
  dirib[:,2] .= true
  dirib[:,end-1] .= true
  # time steps
  for k = 1:steps
    if k % 100 == 1
      println(@sprintf("%6i %8.3e", k-1, maximum(hypot.(u[2:nx-1,2:ny-1], v[2:nx-1,2:ny-1]))))
      # vorticity
#      ωz = f .+ (v[3:nx,2:ny-1] - v[1:nx-2,2:ny-1] - u[2:nx-1,3:ny] + u[2:nx-1,1:ny-2])/2h
#      ωz[.!maski[2:nx-1,2:ny-1]] .= NaN
      # buoyancy
      bs = b[2:nx-1,2:ny-1]
      bs[((maskx .== 0) .& (masky .== 0))[2:nx-1, 2:ny-1]] .= NaN
      # save data
      open(@sprintf("sound/data/%1d_%1d_%010d", i, j, k-1), "w") do file
	write(file, bs)#ωz)
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
#    # forcing
#    forcing!((2k+.5)*Δt, u, maski)
    # Rz step
    u, v = Rz(u, v, ϕ, maski, channel_send_west, channel_send_east, channel_send_south, channel_send_north, channel_receive_west,
        channel_receive_east, channel_receive_south, channel_receive_north)
    # Sx step
    u, ϕ = Sx(u, ϕ, maskx, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    # Sy step
    v, ϕ = Sy(v, ϕ, masky, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    # Tx split
    u, b = Tx(u, ϕ, b, maskx, dirib, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    # Ty split
    v, b = Ty(v, ϕ, b, masky, dirib, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    # Ty split
    v, b = Ty(v, ϕ, b, masky, dirib, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    # Tx split
    u, b = Tx(u, ϕ, b, maskx, dirib, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    # Sy step
    v, ϕ = Sy(v, ϕ, masky, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    # Sx step
    u, ϕ = Sx(u, ϕ, maskx, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    # Rz step
    u, v = Rz(u, v, ϕ, maski, channel_send_west, channel_send_east, channel_send_south, channel_send_north, channel_receive_west,
	      channel_receive_east, channel_receive_south, channel_receive_north)
#    # forcing
#    forcing!((2k+1.5)*Δt, u, maski)
    # Dy steps
    u = Dy(u, masky, diriu, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    v = Dy(v, masky, diriv, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    b = Dy(b, masky, dirib, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    # Dx steps
    u = Dx(u, maskx, diriu, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    v = Dx(v, maskx, diriv, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    b = Dx(b, maskx, dirib, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
#    # print conservation diagnostics
#    println(@sprintf("%16.10e", mass(ϕ[2:nx-1,2:ny-1], maski[2:nx-1,2:ny-1], maskx[2:nx-1,2:ny-1], masky[2:nx-1,2:ny-1])))
#    println(@sprintf("%16.10e", energy(u[2:nx-1,2:ny-1], v[2:nx-1,2:ny-1], ϕ[2:nx-1,2:ny-1], b[2:nx-1,2:ny-1],
#				       [j*h for i = 2:nx-1, j = 2:ny-1], maski[2:nx-1,2:ny-1], maskx[2:nx-1,2:ny-1],
#				       masky[2:nx-1,2:ny-1])))
#    println(@sprintf("%16.10e", buoyancy(ϕ[2:nx-1,2:ny-1], b[2:nx-1,2:ny-1], maski[2:nx-1,2:ny-1], maskx[2:nx-1,2:ny-1],
#					 masky[2:nx-1,2:ny-1])))
#    println(@sprintf("%16.10e", buoyancy_variance(ϕ[2:nx-1,2:ny-1], b[2:nx-1,2:ny-1], maski[2:nx-1,2:ny-1], maskx[2:nx-1,2:ny-1],
#						  masky[2:nx-1,2:ny-1])))
  end
  return
end

# get index ranges for tiles
tile_range(i, j, tile_sizes) = [sum(tile_sizes[1][1:i-1])+1:sum(tile_sizes[1][1:i]),
				sum(tile_sizes[2][1:j-1])+1:sum(tile_sizes[2][1:j])]

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

#"""Mass"""
@everywhere mass(ϕ, maski, maskx, masky) = sum(ϕ[maski]) + sum(ϕ[(maskx.>1).|(masky.>1)])/2

#"""Energy"""
@everywhere function energy(u, v, ϕ, b, y, maski, maskx, masky)
  α = ϕ.*b/c^2
  return sum(u[maski].^2)/2 + sum(v[maski].^2)/2 + sum(ϕ[maski].^2)/2c^2 + sum(ϕ[(maskx.>1).|(masky.>1)].^2)/4c^2 - 2sum(α[maski].*y[maski]) - sum(α[masky.>1].*y[masky.>1])
end

# buoyancy
@everywhere buoyancy(ϕ, b, maski, maskx, masky) = sum(ϕ[maski].*b[maski]) + sum(ϕ[(maskx.>1).|(masky.>1)].*b[(maskx.>1).|(masky.>1)])/2

# buoyancy variance
@everywhere buoyancy_variance(ϕ, b, maski, maskx, masky) = sum(ϕ[maski].*b[maski].^2) + sum(ϕ[(maskx.>1).|(masky.>1)].*b[(maskx.>1).|(masky.>1)].^2)/2
