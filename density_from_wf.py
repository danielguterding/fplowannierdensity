import os
import sys
import numpy as np
try:
  import tricubic #please get the pytricubic library from https://github.com/danielguterding/pytricubic
except ImportError:
  print "Error loading pytricubic library. Please get it from  https://github.com/danielguterding/pytricubic."
  sys.exit(0)

def read_dx_file_and_get_values(filename, nx, ny, nz):
  #read input file
  infilehandle = open(filename, 'r')
  lines = infilehandle.readlines()
  infilehandle.close()
  #convert entries to linear array
  wfvalues = []
  for l in lines:
    splitline = l.strip().split()
    for v in splitline:
      wfvalues.append(float(v))
  #reshape array from (i,j,k) to (k,j,i)
  wfvalues = np.array(wfvalues, dtype=float).reshape((nx,ny,nz)).swapaxes(0,2) #now we have values given in order so that x coordinate varies first
  return wfvalues

def get_density(samplegrid, wfvalues, atompos_rel, nx, ny, nz, samplex, sampley, samplez):
  shift = +0.5*np.ones(3, dtype=float)
  sampledens = np.zeros(len(samplegrid), dtype=float)
  for w,p in zip(wfvalues, atompos_rel):
    #initialize interpolator for this partial density
    ip = tricubic.tricubic(list(w), [nx, ny, nz])
    for i,vec in enumerate(samplegrid):
      transformedvec = vec + shift + p[::-1]
      transformedvec = np.remainder(transformedvec, 1) #transform coordinates and fold back to range 0 ... 1
      transformedvec = list(np.multiply(transformedvec, np.array([nx-1, ny-1, nz-1], dtype=float))) #multiply by n-1, because tricubic interpolator has fixed spacing of one between data points
      sampledens[i] += (ip.ip(transformedvec))**2 #square to get density from WF
  return sampledens      
  
def get_coordinate_grid_xz_plane(nx, nz, y):
  coordinates = []
  for x in np.linspace(0,1,num=nx,endpoint=True):
    for z in np.linspace(0,1,num=nz,endpoint=True):
      vec = np.array([x, y, z], dtype=float)
      coordinates.append(vec)
  return coordinates

def write_cut_xz(filename, coordinates, density, samplex, samplez):
  outfilehandle = open(filename, 'w')
  outfilehandle.write('#x y z density\n')
  for i in range(samplex):
    for j in range(samplez):
      k = i*samplez + j
      c = coordinates[k]
      d = density[k]
      outstr = '%f %f %f %f\n' % (c[0], c[1], c[2], d)
      outfilehandle.write(outstr)
    outfilehandle.write('\n')
  outfilehandle.close()

def main():
  filelist = ['As_4pz/4p24GPa/wfdata028', 'As_4pz/4p24GPa/wfdata031'] #list of wannier function inputs
  outfilename = 'denscutxz_4p24GPa'
  nx, ny, nz = [51, 51, 51] #number of points in original grid
  samplex, sampley, samplez = [300, 100, 300] #desired resolution of sample grid
  ypos = 0.5 #y-coordinate at which to evaluate xz-cut
  atompos_rel = [np.array([0.5, 0.5, 0.37678]), np.array([0.5, 0.5, 0.62322])] #positions corresponding to input files
  
  wfvalues = []
  for f in filelist:
    wfvalues.append(read_dx_file_and_get_values(f, nx, ny, nz))
  
  samplegrid = get_coordinate_grid_xz_plane(samplex, samplez, ypos)
  density = get_density(samplegrid, wfvalues, atompos_rel, nx, ny, nz, samplex, sampley, samplez)
  write_cut_xz(outfilename, samplegrid, density, samplex, samplez)
  
  return 0
  
main()