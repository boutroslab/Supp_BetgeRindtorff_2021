import numpy as np
import h5py
import imageio

#rsync -rav rindtorf@b110-sc2cn01://data/valentini/promise/PROMISE/data-10x-4t-c-16z/hdf5projection/D018T01P907L08/D018T01P907L08_E_13_contrastProjections.h5 ~/Downloads

# path = '/data/valentini/promise/PROMISE/data-10x-4t-c-16z/hdf5projection/D018T01P907L08/D018T01P907L08_E_13_contrastProjections.h5'
path = 'Downloads/D018T01P907L08_E_13_contrastProjections.h5'
prefix = 'Downloads/'
fld = [0,1,2,3]
ch = [0,1,2]


hf = h5py.File(path, 'r')
print(hf.keys())
set = hf['images']

for i in fld:
  for j in ch:
    imageio.imwrite('{}export_fld{}_ch{}.tif'.format(prefix, i+1, j+1), set[i,j])

