@everywhere using Distributed
@everywhere using Printf
@everywhere using PyPlot
@everywhere using HDF5

const cw = [RemoteChannel(()->Channel{Array{Bool}}(1)) for i = 1:2, j = 1:2]
const ce = [RemoteChannel(()->Channel{Array{Bool}}(1)) for i = 1:2, j = 1:2]
const cs = [RemoteChannel(()->Channel{Array{Bool}}(1)) for i = 1:2, j = 1:2]
const cn = [RemoteChannel(()->Channel{Array{Bool}}(1)) for i = 1:2, j = 1:2]

@everywhere function step!(a)
  # perform step of Conway's game of life
  n, m = size(a)
  # number of neighbor cells that are alive
  N = [sum(a[i-1:i+1,j-1:j+1]) - a[i,j] for i = 2:n-1, j = 2:m-1]
  for i = 2:n-1
    for j = 2:m-1
      if a[i,j] & ((N[i-1,j-1] < 2) | (N[i-1,j-1] > 3)) # death by under- or overpopulation
	a[i,j] = false
      elseif !a[i,j] & (N[i-1,j-1] == 3) # reproduction
	a[i,j] = true
      end
    end
  end
end

@everywhere function run_tile(tile, irange, jrange, steps, cmw, cme, cms, cmn, crw, cre, crs, crn)
  # initialize
  nx = length(irange)
  ny = length(jrange)
  a = Array{Bool}(undef, nx+2, ny+2)
  a[2:nx+1,2:ny+1] = (h5read("caltech_990.h5", "m") .> .5)[irange,jrange]
  for i = 1:steps
    println(i)
    # send edges
    println("writing for W edge")
    put!(cmw, a[2,:])
    println("done")
    println("writing for E edge")
    put!(cme, a[nx+1,:])
    println("done")
    println("writing for S edge")
    put!(cms, a[:,2])
    println("done")
    println("writing for N edge")
    put!(cmn, a[:,ny+1])
    println("done")
    # receive edges
    println("reading for W edge")
    a[1,:] = take!(crw)
    println("done")
    println("reading for E edge")
    a[nx+2,:] = take!(cre)
    println("done")
    println("reading for S edge")
    a[:,1] = take!(crs)
    println("done")
    println("reading for N edge")
    a[:,ny+2] = take!(crn)
    println("done")
    # save image
    imsave(@sprintf("gol/%1d_%04d.png", tile, i), Array(a[2:nx+1,2:ny+1]'), origin="lower")
    # perform step
    step!(a)
  end
  return a
end

function merge_images(steps)
  @distributed for i = 1:steps
    println(i)
    f1 = @sprintf("gol/1_%04d.png", i)
    f2 = @sprintf("gol/2_%04d.png", i)
    f3 = @sprintf("gol/3_%04d.png", i)
    f4 = @sprintf("gol/4_%04d.png", i)
    out = @sprintf("gol/merged/%04d.png", i)
    run(`montage $f3 $f4 $f1 $f2 -tile 2x2 -geometry +0+0 $out`)
  end
end

function run_gol(steps)
  a1 = @spawnat 2 run_tile(1,   1:445,   1:211, steps, cw[1,1], ce[1,1], cs[1,1], cn[1,1], ce[2,1], cw[2,1], cn[1,2], cs[1,2])
  a2 = @spawnat 3 run_tile(2, 446:990,   1:211, steps, cw[2,1], ce[2,1], cs[2,1], cn[2,1], ce[1,1], cw[1,1], cn[2,2], cs[2,2])
  a3 = @spawnat 4 run_tile(3,   1:445, 212:422, steps, cw[1,2], ce[1,2], cs[1,2], cn[1,2], ce[2,2], cw[2,2], cn[1,1], cs[1,1])
  a4 = @spawnat 5 run_tile(4, 446:990, 212:422, steps, cw[2,2], ce[2,2], cs[2,2], cn[2,2], ce[1,2], cw[1,2], cn[2,1], cs[2,1])
  wait(a1)
  wait(a2)
  wait(a3)
  wait(a4)
  merge_images(steps)
end
