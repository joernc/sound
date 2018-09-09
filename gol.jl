# run, e.g., as
#   run_gol(16, [repeat([198], 5), repeat([211], 2)])

@everywhere using Distributed
@everywhere using Printf
@everywhere using PyPlot
@everywhere using HDF5

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

@everywhere function run_tile(i, j, irange, jrange, steps, cmw, cme, cms, cmn, crw, cre, crs, crn)
  # initialize
  nx = length(irange)
  ny = length(jrange)
  a = Array{Bool}(undef, nx+2, ny+2)
  a[2:nx+1,2:ny+1] = (h5read("caltech_990.h5", "m") .> .5)[irange,jrange]
  for k = 1:steps
    println(k)
    # send edges
    put!(cmw, a[2,:])
    put!(cme, a[nx+1,:])
    put!(cms, a[:,2])
    put!(cmn, a[:,ny+1])
    # receive edges
    a[1,:] = take!(crw)
    a[nx+2,:] = take!(cre)
    a[:,1] = take!(crs)
    a[:,ny+2] = take!(crn)
    # save data
    open(@sprintf("gol/data/%1d_%1d_%04d", i, j, k), "w") do file
      write(file, a[2:nx+1,2:ny+1])
    end
    # perform step
    step!(a)
  end
  return a
end

# get index ranges for tiles
tile_range(i, j, tile_sizes) = sum(tile_sizes[1][1:i-1])+1:sum(tile_sizes[1][1:i]), sum(tile_sizes[2][1:j-1])+1:sum(tile_sizes[2][1:j])

function run_gol(steps, tile_sizes)
  # number of tiles
  ni = length(tile_sizes[1])
  nj = length(tile_sizes[2])
  # set up remote channels for edge communication
  cw = [RemoteChannel(()->Channel{Array{Bool}}(1)) for i = 1:ni, j = 1:nj]
  ce = [RemoteChannel(()->Channel{Array{Bool}}(1)) for i = 1:ni, j = 1:nj]
  cs = [RemoteChannel(()->Channel{Array{Bool}}(1)) for i = 1:ni, j = 1:nj]
  cn = [RemoteChannel(()->Channel{Array{Bool}}(1)) for i = 1:ni, j = 1:nj]
  # spawn work on tiles
  a = Array{Future}(undef, ni, nj)
  for i in 1:ni
    for j in 1:nj
      irange, jrange = tile_range(i, j, tile_sizes)
      p = workers()[mod1(((i-1)*ni+j), nprocs()-1)]
      println(p, ": ", i, " ", j, " ", irange, " ", jrange)
      a[i,j] = @spawnat p run_tile(i, j, irange, jrange, steps,
          cw[i,j], ce[i,j], cs[i,j], cn[i,j], ce[mod1(i-1,ni),j], cw[mod1(i+1,ni),j], cn[i,mod1(j-1,nj)], cs[i,mod1(j+1,nj)])
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
    println(k)
    a = Array{Bool}(undef, sum(tile_sizes[1]), sum(tile_sizes[2]))
    for i in 1:ni
      for j in 1:nj
	irange, jrange = tile_range(i, j, tile_sizes)
	open(@sprintf("gol/data/%1d_%1d_%04d", i, j, k), "r") do file
	  a[irange,jrange] = read(file)
	end
      end
    end
    imsave(@sprintf("gol/fig/%04d.png", k), Array(a'), origin="lower")
  end
end
