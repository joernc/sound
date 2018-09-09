# Implementation of Rick Salmon's "An Ocean Circulation Model Based on Operator-Splitting, Hamiltonian Brackets, and the Inclusion of Sound Waves" (JPO, 2009)

# next steps:
# - add buoyancy
# - extend to third dimension
# - cut up into tiles to allow parallelization

using Printf
using HDF5
using PyPlot

# grid points
const nx = 990 + 14 # add twice the overlap
const ny = 422 + 1 # add 1 for solid boundary

# overlap
const ox = 7
const oy = 0

# grid spacing (domain size: 2pi)
const h = 2pi/nx

# sound speed
const c = 1.

# inertial period
const f = 0.

# time step
const dt = h/c

# viscosity
const nu = 2e-5

# parameter for diffusion
const alp = nu/(c*h)

"""Sound wave split (x-direction)"""
function Sx(u, p, mx)
  up = Array{Float64,2}(undef, nx, ny)
  pp = Array{Float64,2}(undef, nx, ny)
  for i = 1:nx
    for j = 1:ny
      if mx[i,j] == 1 # interior
	up[i,j] = (u[i-1,j] + u[i+1,j])/2 + (p[i-1,j] - p[i+1,j])/2c
	pp[i,j] = (p[i-1,j] + p[i+1,j])/2 + c/2*(u[i-1,j] - u[i+1,j])
      elseif mx[i,j] == 2 # west boundary
	up[i,j] = 0.
	pp[i,j] = p[i+1,j] - c*u[i+1,j]
      elseif mx[i,j] == 3 # east boundary
	up[i,j] = 0.
	pp[i,j] = p[i-1,j] + c*u[i-1,j]
      else # solid
	up[i,j] = u[i,j]
	pp[i,j] = p[i,j]
      end
    end
  end
  return up, pp
end

"""Sound wave split (y-direction)"""
function Sy(v, p, my)
  vp = Array{Float64,2}(undef, nx, ny)
  pp = Array{Float64,2}(undef, nx, ny)
  for i = 1:nx
    for j = 1:ny
      if my[i,j] == 1 # interior
	vp[i,j] = (v[i,j-1] + v[i,j+1])/2 + (p[i,j-1] - p[i,j+1])/2c
	pp[i,j] = (p[i,j-1] + p[i,j+1])/2 + c/2*(v[i,j-1] - v[i,j+1])
      elseif my[i,j] == 2 # south boundary
	vp[i,j] = 0.
	pp[i,j] = p[i,j+1] - c*v[i,j+1]
      elseif my[i,j] == 3 # north boundary
	vp[i,j] = 0.
	pp[i,j] = p[i,j-1] + c*v[i,j-1]
      else # solid
	vp[i,j] = v[i,j]
	pp[i,j] = p[i,j]
      end
    end
  end
  return vp, pp
end

"Diffusion split (2D)"
function Dxy(a, m)
  ap = Array{Float64,2}(undef, nx, ny)
  for i = 1:nx
    for j = 1:ny
      if m[i,j]
#	ap[i,j] = (a[i,j] + alp*(a[i-1,j] + a[i+1,j] + a[i,j-1] + a[i,j+1]))/(1 + 4alp)
	ap[i,j] = (1-4alp)*a[i,j] + alp*(a[i-1,j] + a[i+1,j] + a[i,j-1] + a[i,j+1])
      else
	ap[i,j] = a[i,j]
      end
    end
  end
  return ap
end

"""Rotation split"""
function Rz(u, v, p, m)
  up = Array{Float64,2}(undef, nx, ny)
  vp = Array{Float64,2}(undef, nx, ny)
  for i = 1:nx
    for j = 1:ny
      if m[i,j] # interior
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

"""Make periodic (in x-dicrection)"""
function periodicize_x!(a)
  a[1:ox,:] = a[nx-2ox+1:nx-ox,:]
  a[nx-ox+1:nx,:] = a[ox+1:2ox,:]
end

"""Make periodic (in y-dicrection)"""
function periodicize_y!(a)
  a[:,1:oy] = a[:,ny-2oy+1:ny-oy]
  a[:,ny-oy+1:ny] = a[:,oy+1:2oy]
end

"""Apply forcing"""
function forcing!(t, u, m)
  u[m] .+= .01*.01*cos(.01t)*dt
end

"""Time step (Strang splitting)"""
function time_step(t, u, v, p, m, mx, my)
  # apply splits
  u = Dxy(u, m)
  v = Dxy(v, m)
  forcing!(t+dt/2, u, m)
  u, v = Rz(u, v, p, m)
  u, p = Sx(u, p, mx)
  v, p = Sy(v, p, my)
  v, p = Sy(v, p, my)
  u, p = Sx(u, p, mx)
  u, v = Rz(u, v, p, m)
  forcing!(t+3dt/2, u, m)
  u = Dxy(u, m)
  v = Dxy(v, m)
  # periodicize
  periodicize_x!(u)
  periodicize_x!(v)
  periodicize_x!(p)
  periodicize_y!(u)
  periodicize_y!(v)
  periodicize_y!(p)
  return u, v, p
end

"""Mass"""
mass(p, m, mx, my) = sum(p[m]) + sum(p[(mx.>1).|(my.>1)])/2

"""Energy"""
energy(u, v, p, m, mx, my) = sum(u[m].^2)/2 + sum(v[m].^2)/2 + sum(p[m].^2)/2c^2 + sum(p[(mx.>1).|(my.>1)].^2)/4c^2

"""Boundary masks (see Fig. 1 in Salmon for a sketch)"""
function boundary_masks(s)
  # pad s with overlap
  s = [s; s[1:2ox,:]]
  s = [s s[:,1:2oy]]
  # east–west boundaries
  mx = zeros(UInt8, nx, ny)
  for j = 2:ny-1
    if s[1,j-1] & s[1,j]
      mx[1,j] = 2 # western boundary
    end
    if s[nx-1,j-1] & s[nx-1,j]
      mx[nx,j] = 3 # eastern boundary
    end
  end
  for i = 2:nx-1
    for j = 2:ny-1
      if s[i-1,j-1] & s[i,j-1] & s[i-1,j] & s[i,j]
	mx[i,j] = 1 # interior
      elseif (!s[i-1,j-1] | !s[i-1,j]) & s[i,j-1] & s[i,j]
	mx[i,j] = 2 # western boundary
      elseif (!s[i,j-1] | !s[i,j]) & s[i-1,j-1] & s[i-1,j]
	mx[i,j] = 3 # eastern boundary
      end
    end
  end
  # north–south boundaries
  my = zeros(UInt8, nx, ny)
  for i = 2:nx-1
    if s[i-1,1] & s[i,1]
      my[i,1] = 2 # southern boundary
    end
    if s[i-1,ny-1] & s[i,ny-1]
      my[i,ny] = 3 # northern boundary
    end
  end
  for i = 2:nx-1
    for j = 2:ny-1
      if s[i-1,j-1] & s[i-1,j] & s[i,j-1] & s[i,j]
	my[i,j] = 1 # interior
      elseif (!s[i-1,j-1] | !s[i,j-1]) & s[i-1,j] & s[i,j]
	my[i,j] = 2 # southern boundary
      elseif (!s[i-1,j] | !s[i,j]) & s[i-1,j-1] & s[i,j-1]
	my[i,j] = 3 # northern boundary
      end
    end
  end
  # interior points
  m = (mx .== 1) .& (my .== 1)
  return m, mx, my
end

# save every se time steps
global se = 200

function run_model()
  # topography mask
  s = h5read("caltech_990.h5", "m") .< .5
  m, mx, my = boundary_masks(s)
  # initial conditions
  u = zeros(nx, ny)
  v = zeros(nx, ny)
  p = c^2*ones(nx, ny)
  # time stepping
  for i = 1:100000
    # print diagnostics
    idx = CartesianIndices((1+ox:nx-ox,1+oy:ny-oy))
    @printf("%6d %6d %12.10e %12.10e %12.10e\n", i-1, (i-1)*se, (i-1)*se*2dt, mass(p[idx], m[idx], mx[idx], my[idx]),
        energy(u[idx], v[idx], p[idx], m[idx], mx[idx], my[idx]))
    # save plots
    idx = CartesianIndices((2+ox:nx-1-ox,2+oy:ny-1-oy))
    omgz = f .+ (v[idx.+CartesianIndex(1,0)] - v[idx.-CartesianIndex(1,0)] - u[idx.+CartesianIndex(0,1)] + u[idx.-CartesianIndex(0,1)]) / 2h
    omgz[.!m[idx]] .= 0.
    mm = maximum(abs.(omgz))
    imsave(@sprintf("fig_caltech/vort/%010d.png", i-1), Array(omgz'), vmin=-mm, vmax=mm, origin="lower")
    idx = CartesianIndices((1+ox:nx-ox,1+oy:ny-oy))
    imsave(@sprintf("fig_caltech/u/%010d.png", i-1), Array(u[idx]'), vmin=-.01, vmax=.01, origin="lower")
    imsave(@sprintf("fig_caltech/v/%010d.png", i-1), Array(v[idx]'), vmin=-.01, vmax=.01, origin="lower")
    imsave(@sprintf("fig_caltech/p/%010d.png", i-1), Array(p[idx]'), origin="lower")
    # step se time steps
    for j = 1:se
      u, v, p = time_step(((i-1)*se+j-1)*2dt, u, v, p, m, mx, my)
    end

  end
  return u, v, p
end

# run the model
u, v, p = run_model()
