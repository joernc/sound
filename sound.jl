# Implementation of Rick Salmon's "An Ocean Circulation Model Based on Operator-Splitting, Hamiltonian Brackets, and the Inclusion of Sound Waves" (JPO, 2009)

# next steps:
# - add buoyancy
# - extend to third dimension
# - organize channels better

@everywhere using Printf
@everywhere using HDF5
@everywhere using PyPlot

# grid spacing (domain size: 2pi)
@everywhere const h = 2pi/512

# sound speed
@everywhere const c = 1.

# inertial period
@everywhere const f = 0.

# time step
@everywhere const dt = h/c

# viscosity
@everywhere const nu = 2e-5

# parameter for diffusion
@everywhere const alp = nu/(c*h)

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
@everywhere function Sx(u, p, maskx, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangex!(u, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangex!(p, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  nx, ny = size(u)
  up = Array{Float64, 2}(undef, nx, ny)
  pp = Array{Float64, 2}(undef, nx, ny)
  for i = 2:nx-1
    for j = 2:ny-1
      if maskx[i,j] == 1 # interior
	up[i,j] = (u[i-1,j] + u[i+1,j])/2 + (p[i-1,j] - p[i+1,j])/2c
	pp[i,j] = (p[i-1,j] + p[i+1,j])/2 + c/2*(u[i-1,j] - u[i+1,j])
      elseif maskx[i,j] == 2 # western boundary
	up[i,j] = 0.
	pp[i,j] = p[i+1,j] - c*u[i+1,j]
      elseif maskx[i,j] == 3 # eastern boundary
	up[i,j] = 0.
	pp[i,j] = p[i-1,j] + c*u[i-1,j]
      else
	up[i,j] = u[i,j]
	pp[i,j] = p[i,j]
      end
    end
  end
  return up, pp
end

#"""Sound wave split (y-direction)"""
@everywhere function Sy(v, p, masky, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  exchangey!(v, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  exchangey!(p, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  nx, ny = size(v)
  vp = Array{Float64, 2}(undef, nx, ny)
  pp = Array{Float64, 2}(undef, nx, ny)
  for i = 2:nx-1
    for j = 2:ny-1
      if masky[i,j] == 1 # interior
	vp[i,j] = (v[i,j-1] + v[i,j+1])/2 + (p[i,j-1] - p[i,j+1])/2c
	pp[i,j] = (p[i,j-1] + p[i,j+1])/2 + c/2*(v[i,j-1] - v[i,j+1])
      elseif masky[i,j] == 2 # south boundary
	vp[i,j] = 0.
	pp[i,j] = p[i,j+1] - c*v[i,j+1]
      elseif masky[i,j] == 3 # north boundary
	vp[i,j] = 0.
	pp[i,j] = p[i,j-1] + c*v[i,j-1]
      else
	vp[i,j] = v[i,j]
	pp[i,j] = p[i,j]
      end
    end
  end
  return vp, pp
end

#"""Rotation split"""
@everywhere function Rz(u, v, p, maski, channel_send_west, channel_send_east, channel_send_south, channel_send_north, channel_receive_west, channel_receive_east, channel_receive_south, channel_receive_north)
  exchangex!(u, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangex!(v, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangey!(u, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  exchangey!(v, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  nx, ny = size(u)
  up = Array{Float64, 2}(undef, nx, ny)
  vp = Array{Float64, 2}(undef, nx, ny)
  for i = 2:nx-1
    for j = 2:ny-1
      if maski[i,j] # interior
        omgz = f + (v[i+1,j] - v[i-1,j] - u[i,j+1] + u[i,j-1]) / 2h
        gam = dt*omgz*c^2/p[i,j]
        Cz = (1 - gam.^2/4)/(1 + gam^2/4)
        Sz = gam/(1 + gam^2/4)
        up[i,j] = u[i,j]*Cz + v[i,j]*Sz
	vp[i,j] = v[i,j]*Cz - u[i,j]*Sz
      else # boundary or solid
	up[i,j] = u[i,j]
	vp[i,j] = v[i,j]
      end
    end
  end
  return up, vp
end

#"""Apply forcing"""
@everywhere function forcing!(t, u, maski)
  u[maski] .+= .01*.01*cos(.01t)*dt
end

#"Diffusion split (2D)"
@everywhere function Dxy(a, maski, channel_send_west, channel_send_east, channel_send_south, channel_send_north,
    channel_receive_west, channel_receive_east, channel_receive_south, channel_receive_north)
  exchangex!(a, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangey!(a, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  nx, ny = size(a)
  ap = Array{Float64, 2}(undef, nx, ny)
  for i = 2:nx-1
    for j = 2:ny-1
      if maski[i,j]
	ap[i,j] = (1-4alp)*a[i,j] + alp*(a[i-1,j] + a[i+1,j] + a[i,j-1] + a[i,j+1])
      else
	ap[i,j] = a[i,j]
      end
    end
  end
  return ap
end

#"""Boundary masks (see Fig. 1 in Salmon for a sketch)"""
@everywhere function boundary_masks(fluid)
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
  return maski, maskx, masky
end

@everywhere function run_tile(i, j, irange, jrange, steps, channel_send_west, channel_send_east, channel_send_south,
    channel_send_north, channel_receive_west, channel_receive_east, channel_receive_south, channel_receive_north)
  # topography mask
  nx = length(irange) + 2
  ny = length(jrange) + 2
  fluid = Array{Bool}(undef, nx, ny)
  fluid[2:nx-1,2:ny-1] = (h5read("julia_256_512.h5", "m") .> .5)[irange,jrange]
  exchangex!(fluid, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
  exchangey!(fluid, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
  maski, maskx, masky = boundary_masks(fluid)
  # initialization
  u = zeros(nx, ny)
  v = zeros(nx, ny)
  p = c^2*ones(nx, ny)
  u[(maskx.==0).&(masky.==0)] .= 0
  v[(maskx.==0).&(masky.==0)] .= 0
  p[(maskx.==0).&(masky.==0)] .= 0
  # time steps
  for k = 1:steps
    println(k)
    if k % 100 == 1
      # save data
      imsave(@sprintf("sound/fig/%1d_%1d_%010d.png", i, j, k-1), Array(u[2:nx-1,2:ny-1]'), origin="lower")
      open(@sprintf("sound/data/%1d_%1d_%010d", i, j, k-1), "w") do file
	write(file, u[2:nx-1,2:ny-1])
      end
    end
    # Dxy steps
    u = Dxy(u, maski, channel_send_west, channel_send_east, channel_send_south, channel_send_north, channel_receive_west,
        channel_receive_east, channel_receive_south, channel_receive_north)
    v = Dxy(v, maski, channel_send_west, channel_send_east, channel_send_south, channel_send_north, channel_receive_west,
        channel_receive_east, channel_receive_south, channel_receive_north)
    # Sx step
    u, p = Sx(u, p, maskx, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    # Sy step
    v, p = Sy(v, p, masky, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    # Rz step
    u, v = Rz(u, v, p, maski, channel_send_west, channel_send_east, channel_send_south, channel_send_north, channel_receive_west,
        channel_receive_east, channel_receive_south, channel_receive_north)
    # forcing
    forcing!((2k+.5)*dt, u, maski)
    # forcing
    forcing!((2k+1.5)*dt, u, maski)
    # Rz step
    u, v = Rz(u, v, p, maski, channel_send_west, channel_send_east, channel_send_south, channel_send_north, channel_receive_west,
        channel_receive_east, channel_receive_south, channel_receive_north)
    # Sy step
    v, p = Sy(v, p, masky, channel_send_south, channel_send_north, channel_receive_south, channel_receive_north)
    # Sx step
    u, p = Sx(u, p, maskx, channel_send_west, channel_send_east, channel_receive_west, channel_receive_east)
    # Dxy steps
    u = Dxy(u, maski, channel_send_west, channel_send_east, channel_send_south, channel_send_north, channel_receive_west,
        channel_receive_east, channel_receive_south, channel_receive_north)
    v = Dxy(v, maski, channel_send_west, channel_send_east, channel_send_south, channel_send_north, channel_receive_west,
        channel_receive_east, channel_receive_south, channel_receive_north)
  end
  return
end

# get index ranges for tiles
tile_range(i, j, tile_sizes) = sum(tile_sizes[1][1:i-1])+1:sum(tile_sizes[1][1:i]),
    sum(tile_sizes[2][1:j-1])+1:sum(tile_sizes[2][1:j])

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
      p = workers()[mod1(((i-1)*ni+j), nprocs()-1)]
      println(p, ": ", i, " ", j, " ", irange, " ", jrange)
      a[i,j] = @spawnat p run_tile(i, j, irange, jrange, steps, channels_west[i,j], channels_east[i,j], channels_south[i,j],
          channels_north[i,j], channels_east[mod1(i-1,ni),j], channels_west[mod1(i+1,ni),j], channels_north[i,mod1(j-1,nj)],
          channels_south[i,mod1(j+1,nj)])
    end
  end
  # wait for retult
  for i in 1:ni
    for j in 1:nj
      wait(a[i,j])
    end
  end
  # assemble and plot results
  for k in 1:steps
    if k % 100 == 1
      println(k-1)
      u = Array{Float64}(undef, sum(tile_sizes[1]), sum(tile_sizes[2]))
      for i in 1:ni
	for j in 1:nj
	  irange, jrange = tile_range(i, j, tile_sizes)
	  open(@sprintf("sound/data/%1d_%1d_%010d", i, j, k-1), "r") do file
	    u[irange,jrange] = [read(file, Float64) for ii in irange, jj in jrange]
	  end
	end
      end
      imsave(@sprintf("sound/fig/%010d.png", k-1), Array(u'), origin="lower")
    end
  end
end

#"""Mass"""
#mass(p, m, mx, my) = sum(p[m]) + sum(p[(mx.>1).|(my.>1)])/2
#
#"""Energy"""
#energy(u, v, p, m, mx, my) = sum(u[m].^2)/2 + sum(v[m].^2)/2 + sum(p[m].^2)/2c^2 + sum(p[(mx.>1).|(my.>1)].^2)/4c^2
