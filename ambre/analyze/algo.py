'''
Created on May 23, 2012

@author: root
'''
import sys
import numpy as na
import bisect
from collections import defaultdict
from operator import itemgetter
debug_log = sys.stderr

class OverlapBPs(object):
  def __init__(self):
    self.overlaps = None
    self.a = None
    self.b = None
    self.d = None
  def _add_intersection_q(self, a_idx, b_idx):
    if a_idx<=b_idx:
      self.overlaps[a_idx].add(b_idx)
    else:
      self.overlaps[b_idx].add(a_idx)
  def _add_intersection(self, idxs):
    for j in xrange(len(idxs)):
      for i in xrange(1,j):
        if idxs[j]<=idxs[i]:
          self.overlaps[idxs[j]].add(idxs[i])
        else:
          self.overlaps[idxs[i]].add(idxs[j])
  def _get_next_sweep_line(self, c, sorted_c_idx):
    prev_c = None
    prev_c_list = []
    prev_a_list = []
    for c_idx in sorted_c_idx:
      c_v = c[c_idx]
      if not prev_c is None and prev_c!=c_v:
        yield prev_c, prev_c_list, prev_a_list
        prev_c_list = []
        prev_a_list = []
      
      if len(prev_c_list)==0:
        prev_c_list.append(c_idx)
        prev_a_list.append(self.a[c_idx]-self.d[c_idx])
      else:
        i = bisect.bisect(prev_a_list, self.a[c_idx]-self.d[c_idx])
        prev_c_list.insert(i, c_idx)
        prev_a_list.insert(i, self.a[c_idx]-self.d[c_idx])
        
      prev_c = c_v  
    yield prev_c, prev_c_list, prev_a_list
  
  def _add_overlaps_on_same_sweep_line(self, c_list, a_list):
    if len(c_list)<2:
      return 
    # Assumes elements on sweep are ordered. 
    alive = []
    for i,s,t in zip(c_list, a_list, self.a[c_list]):
      p = bisect.bisect(alive, (s,None))
      alive = alive[p+1:]
      if len(alive)>0:
        for u,v in zip(*alive):
          self._add_intersection_q(i, v)
      alive.append((t,i))
    
  def _add_overlaps_from_new_sweep_line(self, prev_c_list, prev_a_list, c_list, a_list):
    
    new_sweep = zip(c_list, a_list, self.a[c_list])
    sweep = zip(prev_a_list, self.a[prev_c_list], prev_c_list)
    if len(sweep)==0 or len(new_sweep)==0:
      return
    for i,s,t in new_sweep:
      u_u = bisect.bisect(sweep, (s,None,None))
      u_v = bisect.bisect(sweep, (t,self.max_a,None))
      for u,v,j in sweep[u_u:u_v]:
        self._add_intersection_q(i, j)
       
  def overlapping_triangles(self, element_array):
    #elements are in an (n,3)-array of the form (a \in R ,b \in R,d \in N).
    # Generate a sweep line schedule to iterate over all event starts and ends
    # if two triangles share a segment, [i,j] for i>=j ,of a sweep line then
    # store the elements as shared. 
    
    # All elements are sorted by a-axis.
    
    self.overlaps = defaultdict(set)
    self.a = element_array[:,0]
    self.b = element_array[:,1]
    self.d = element_array[:,2]
     
    self.max_a = na.max(self.a)+1
    # No rotation as it may introduce floating point errors
    c = self.a+self.b-(self.d)
    sorted_c_idx = na.argsort(c)    
    print >>debug_log, "Elements", element_array.shape[0]
    print >>debug_log, "SweepLines", len(set(list(c[sorted_c_idx])))
    
    prev_c_v, prev_c_list, prev_a_list, d_list = None,None, None, None
    for c_v, c_list, a_list in self._get_next_sweep_line(c, sorted_c_idx):
      self._add_overlaps_on_same_sweep_line(c_list, a_list)

      if not prev_c_v is None:
        # update distance from sweep line for previous
        c_v_diff = c_v-prev_c_v
        sorted_list = sorted([(u_c, u_a-(dist-c_v_diff), dist-c_v_diff) for u_c, u_a, dist in zip(prev_c_list, [self.a[prev_c_idx] for prev_c_idx in prev_c_list], d_list) if dist-c_v_diff>=0], key=itemgetter(1))
          
        if len(sorted_list)>0:
          prev_c_list, prev_a_list, d_list = zip(*sorted_list)
          prev_c_list, prev_a_list, d_list = list(prev_c_list), list(prev_a_list), list(d_list) 
          self._add_overlaps_from_new_sweep_line(prev_c_list, prev_a_list, c_list, a_list)
        else:
        	prev_c_list, prev_a_list, d_list = [],[],[]
        for c_idx, a in zip(c_list, a_list):  
          i = bisect.bisect(prev_a_list, a)
          prev_c_list.insert(i, c_idx)
          prev_a_list.insert(i, a)
          d_list.insert(i, self.d[c_idx])
          
      else:
        prev_c_list = c_list
        prev_a_list = a_list
        d_list = list(self.d[c_list])
      prev_c_v = c_v    
  
  
  def get_connected_components(self):
    
    n_overlaps = self.overlaps.copy()
    for i,v in self.overlaps.iteritems():
      for j in v:
        n_overlaps[j].add(i)
    l = zip(sorted(n_overlaps.keys(), reverse=True), [True]*len(n_overlaps))
    
    ccs = []
    processed_keys = set()
    while len(l)>0:
      i,v = l.pop()
      if i in processed_keys:
        continue
      processed_keys.add(i)
      
      l.extend(zip(n_overlaps[i], [False]*len(n_overlaps[i])))
      if v:
        ccs.append([i,])
      else:
        ccs[-1].append(i)
    return ccs
