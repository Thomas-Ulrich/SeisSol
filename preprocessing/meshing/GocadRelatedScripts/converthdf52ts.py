#convert SeisSol hdf5 to ts (Gocad)

# parsing python arguments
import argparse
import h5py
print(h5py.__version__)
import numpy as np
parser = argparse.ArgumentParser(description='convert hdf5 to ts')
parser.add_argument('hdf5_filename', help='hdf5 filename')
parser.add_argument('ts_filename', nargs='?',help='output filname (if not used = hdf5_basename.ts)',default='')
args = parser.parse_args()

#Read Hdf5
h5f = h5py.File(args.hdf5_filename,'r')
print(h5f.keys())
connect = h5f['connect'][:,:] 
xyz = h5f['geometry'][:,:]
#dimensions
nElements = connect.shape[0]
nPoints =  xyz.shape[0]
h5f.close()


if args.ts_filename == '':
   args.ts_filename = args.hdf5_filename[0:-3]+'.ts'

fout = open (args.ts_filename, 'w')

mystring = "GOCAD TSURF 1\n\
HEADER {\n\
name:s1\n\
}\n\
TRIANGLES\n"

for i in range(nPoints):
   mystring = mystring + 'VRTX %d %e %e %e\n' %(i, xyz[i,0], xyz[i,1], xyz[i,2])

fout.write(mystring)

for i in range(nElements):
      fout.write('TRGL %d %d %d\n' %(connect[i,0], connect[i,1],connect[i,2]))
fout.write('END')
fout.close()
print('done writing %s' %args.ts_filename)
