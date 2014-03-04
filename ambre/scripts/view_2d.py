import sys
from matplotlib import pyplot as plt
from ambre.utils import SAM

class PE(object):
  def __init__(self,x,y):
    self.x = x
    self.y = y

def swap(x,y):
  if x>y:
    return y,x
  return x,y

def plot_pe(sam_fpath):
  with open(sam_fpath) as f:
    positions = [int(line.split('\t', 4)[-2]) for line in f]
  
  points = [swap(a,b) for a,b in zip(positions[::2],positions[1::2])]

  x,y = zip(*points)
  plt.scatter(x,y)
  plt.show()

if __name__ == '__main__':
  plot_pe(sys.argv[1])
