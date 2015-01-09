import sys
import numpy as na
from itertools import izip, compress, chain
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
#from matplotlib.markers import TICKLEFT, TICKRIGHT, TICKUP, TICKDOWN
from matplotlib.collections import PathCollection
from matplotlib.path import Path
from ambre.utils import SAM

CONTOUR_RES = 500

Forward = """121513427
121516044
121519217
121522597
121525278
121528631
121532137
121535273
121538152
121541268
121543470
121545998
121548610
121551131
121554451
121557690
121560189
121563255
121566017
121568802
121571125
121573861
121577105
121579657
121582062
121585306
121588748
121591920
121595214
121598522
121601508
121604208
121607596
121611660"""

Reverse = """116313181
116315037
116318961
116321683
116324497
116327819
116330881
116333948
116337099
116339862"""

Forward = map(int, Forward.split('\n'))
Reverse = map(int, Reverse.split('\n'))

class PE(object):
  def __init__(self,x,y):
    self.x = x
    self.y = y
    self.rev_orient_flag = SAM.SAM_FLAGS_H['rev_strand_of_query']
  def get_coord(self):
    return int(self.x.pos), int(self.y.pos)
  def _orient(self, flags):
    return SAM.get_flag(self.rev_orient_flag, int(flags))
  def get_orient(self):
    return self._orient(self.x.flags), self._orient(self.y.flags)

def swap(x,y):
  if int(x)>int(y):
    return True
  return False

def contour_plot(min_x,max_x,min_y,max_y, coords):
  # make these smaller to increase the resolution
  dx, dy = CONTOUR_RES,CONTOUR_RES

  # generate 2 2d grids for the x & y bounds
  y, x = na.mgrid[slice(min_y, max_y + dy, dy),
                slice(min_x, max_x + dx, dx)]

  z = na.zeros(y.shape, dtype=int)
  z = z[:-1, :-1]

  for p1,p2 in coords:
    z[(p2-min_y)/CONTOUR_RES, (p1-min_x)/CONTOUR_RES] += 1

  # x and y are bounds, so z should be the value *inside* those bounds.
  # Therefore, remove the last value from the z array.
  levels = MaxNLocator(nbins=100).bin_boundaries(z.min(), z.max())

  # pick the desired colormap, sensible levels, and define a normalization
  # instance which takes data values and translates those into levels.
  cmap = plt.get_cmap('Greens')
  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

  fig = plt.figure()
  fig.add_subplot(111,  aspect='equal')
  #fig.add_subplot(111,  aspect='equal')
  #plt.axis([x.min(), x.max(), y.min(), y.max()])

  # contours are *point* based plots, so convert our bound into point
  # centers
  plt.contourf(x[:-1, :-1] + dx / 2.,
             y[:-1, :-1] + dy / 2., z, levels=levels,
             cmap=cmap)
  plt.colorbar()
  #plt.show()

def plot_pe(sam_fpath):
  p = SAM.SAMParser()
  reads = list(p.parse(sam_fpath))
  pes = []
  for r1,r2 in izip(reads[::2],reads[1::2]):
    if swap(r1.pos, r2.pos):
      r1,r2 = r2,r1
    pes.append(PE(r1,r2))
  
  coords = [pe.get_coord() for pe in pes]
  x,y = zip(*coords)
  minx,maxx=min(x),max(x)
  miny,maxy=min(y),max(y)
  contour_plot(minx,maxx,miny,maxy,coords)

  n_r,n_f = len(Reverse), len(Forward)
  plt.vlines(Reverse, ymin=miny,ymax=maxy, colors='r', alpha=0.1)
  plt.hlines(Reverse, xmin=miny,xmax=maxy, colors='r', alpha=0.1)
  plt.vlines(Forward, ymin=minx,ymax=maxx, colors='b', alpha=0.1)
  plt.hlines(Forward, xmin=minx,xmax=maxx, colors='b', alpha=0.1)

  read_length = 300
  xorients,yorients = zip(*[pe.get_orient() for pe in pes])
  xorients,yorients = na.array(xorients,dtype=bool), na.array(yorients,dtype=bool)

  points = na.array(coords)
  direct_x = points.copy()
  direct_x[:,0] += read_length
  paths_x = [Path(na.array([i,j])) for i,j in izip(points,direct_x)]


  direct_y = points.copy()
  direct_y[:,1] += read_length
  paths_y = [Path(na.array([i,j])) for i,j in izip(points,direct_y)]
  
  rev_c = PathCollection(list(chain(compress(paths_x,xorients),compress(paths_y,yorients))),
                         edgecolors='r',linewidths=2)
  for_c = PathCollection(list(chain(compress(paths_x,na.logical_not(xorients)),compress(paths_y,na.logical_not(yorients)))),edgecolors='b',linewidths=2)

  ax = plt.gca()
  ax.add_collection(rev_c)
  ax.add_collection(for_c)
  ax.set_xlabel('met')
  ax.set_ylabel('ptprz1')

  """
  x,y = zip(*[pe.get_coord() for pe in compress(pes, xorients)])
  plt.scatter(x,y,marker=TICKLEFT,c='r')
  x,y = zip(*[pe.get_coord() for pe in compress(pes, na.logical_not(xorients))])
  plt.scatter(x,y,marker=TICKRIGHT,c='b')

  x,y = zip(*[pe.get_coord() for pe in compress(pes, yorients)])
  plt.scatter(x,y,marker=TICKDOWN,c='r')
  x,y = zip(*[pe.get_coord() for pe in compress(pes, na.logical_not(yorients))])
  plt.scatter(x,y,marker=TICKUP,c='b')
  """

  plt.show()

if __name__ == '__main__':
  plot_pe(sys.argv[1])
