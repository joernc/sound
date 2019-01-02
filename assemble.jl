using HDF5

# get index ranges for tiles
tile_range(i, j, k, tile_sizes) = [sum(tile_sizes[1][1:i-1])+1:sum(tile_sizes[1][1:i]),
                                   sum(tile_sizes[2][1:j-1])+1:sum(tile_sizes[2][1:j]),
                                   sum(tile_sizes[3][1:k-1])+1:sum(tile_sizes[3][1:k])]

function assemble(experiment, step, variable, tile_sizes)
  # number of tiles
  ni = length(tile_sizes[1])
  nj = length(tile_sizes[2])
  nk = length(tile_sizes[3])
  # total number of grid points
  nx = sum(tile_sizes[1])
  ny = sum(tile_sizes[2])
  nz = sum(tile_sizes[3])
  # assemble masks
  a = Array{Float64, 3}(undef, nx, ny, nz)
  for k = 1:nk, j = 1:nj, i = 1:ni
    irange, jrange, krange = tile_range(i, j, k, tile_sizes)
    filename = @sprintf("/central/groups/oceanphysics/%s/%010d_%1d_%1d_%1d.h5", experiment, step, i, j, k)
    a[irange,jrange,krange] = h5read(filename, variable)
  end
  return a
end
