from ambre.analyze.algo import OverlapBPs
import unittest

class TestOverlap(unittest.TestCase):
  def setUp(self):
    pass
  def test_overlapping_triangles_small(self):
    import numpy as na
    na.random.seed(1)
    
    n = 2000
    data = na.zeros((n,3), dtype=na.int)
    data[:,0] = na.arange(1,2*n+1,2)
    data[:,1] = na.arange(4,2*n+4,2)
    data[:,1] = 0
    data[:,2] = 8
    
    w = OverlapBPs()
    w.overlapping_triangles(data)
    
    print "OverlapCount:", sum([len(u) for u in w.overlaps.values()])
    print w.get_connected_components()
    
    from matplotlib import pyplot as plt
    x,y = [],[]
    for i in xrange(n):
      x.append((data[i,0],data[i,0],data[i,0]-data[i,2]))
      y.append((data[i,1],data[i,1]-data[i,2],data[i,1]))
    call_list = ['-b']*(len(x)*3)
    call_list[::3] = x
    call_list[1::3] = y
    plt.fill(*call_list, 
              alpha=.9, facecolor='b', edgecolor='k', lw=1.2)
    plt.show()
    
  def test_overlapping_triangles(self):
    import numpy as na
    na.random.seed(1)
    
    n = 2000
    data = na.zeros((n,3), dtype=na.int)
    data[:,0] = na.random.random_integers(-1000000, 1000000,n)
    data[:,1] = na.random.random_integers(-1000000, 1000000,n)
    data[:,2] = na.random.random_integers(0,10000,n)
    
    w = OverlapBPs()
    w.overlapping_triangles(data)
    
    print "OverlapCount:", sum([len(u) for u in w.overlaps.values()])
    print w.get_connected_components()

if __name__ == '__main__':
  suite = unittest.TestLoader().loadTestsFromTestCase(TestOverlap)
  unittest.TextTestRunner(verbosity=2).run(suite)
