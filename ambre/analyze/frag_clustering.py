'''
#  ambre.analyze.frag_clustering.py
#
#  Copyright March 2013 by Anand D. Patel
#
#  This program is free software; you may redistribute it and/or modify its
#  under the terms of the GNU General Public License as published by the Free
#  Software Foundation; either version 2 of the License or
#  any later version.
#
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#  for more details.
#
#  For any other inquiries send an email to: Anand D. Patel (adp002@ucsd.edu)
'''
import numpy as na
import sys
from operator import itemgetter
from collections import defaultdict

debug_log = sys.stderr

class BreakPoint(object):
  def __init__(self, support=None):
    self.name = ''
    self.lbp_contig = ''
    self.rbp_contig = ''
    if support is None:
      self.support = []
    else:
      self.support = list(support)
    self.topo_idx = None
  def set_name(self, name):
    self.name = name
    
  def add_support(self, frag_idx, frag_bp_idx, lbp, rbp, d):
    self.support.append((frag_idx, frag_bp_idx, lbp, rbp, d))
  def get_pair_mode(self):
    frag_idx, frag_bp_idx, lbp, rbp, d = zip(*self.support)
    counter = defaultdict(int)
    for p in zip(lbp, rbp):
      counter[p] += 1
    s_counts = sorted(counter.items(), key=itemgetter(1), reverse=True)
    mlbp, mrbp = s_counts[0][0]
    return mlbp, mrbp, self._get_mode(d)
  def _get_mode(self, values):
    v = na.sort(values)
    uniq_v = na.unique(v)
    bins = uniq_v.searchsorted(v)  
    s = zip(na.bincount(bins), uniq_v)
    
    s.sort(reverse=True)
    if len(s)>1 and s[0][0]==s[1][0]:
      return (s[0][1]+s[1][1])*0.5 
    return s[0][1]
  def get_min(self):
    frag_idx, frag_bp_idx, lbp, rbp, d = zip(*self.support)
    return na.min(lbp), na.min(rbp), na.min(d)
  def get_mode(self):
    frag_idx, frag_bp_idx, lbp, rbp, d = zip(*self.support)
    return self._get_mode(lbp), self._get_mode(rbp), self._get_mode(d)
  def get_mean(self):
    frag_idx, frag_bp_idx, lbp, rbp, d = zip(*self.support)
    return na.mean(lbp), na.mean(rbp), na.mean(d)
  def get_tri_vertices(self):
    frag_idx, frag_bp_idx, lbp, rbp, d = zip(*self.support)
    a_lbp, a_rbp, a_d = na.array(lbp), na.array(rbp), na.array(d)
    x = zip(lbp,a_lbp-a_d,lbp)
    y = zip(rbp, rbp, a_rbp-a_d)
    return x,y
  def set_contigs(self, lbp_contig, rbp_contig):
    self.lbp_contig = lbp_contig
    self.rbp_contig = rbp_contig
  def to_file(self, f):
    self.support.sort()
    mlbp,mrbp,md = self.get_pair_mode()
    print >>f, "#%s\t%s\t%s\t%.1f\t%.1f\t%.1f"%(self.name, self.lbp_contig, self.rbp_contig, mlbp,mrbp,md)
    for support in self.support:
      print >>f, "%d\t%d\t%d\t%d\t%d"%support
    print >>f, '='
  def from_file(self, f):
    while True:
      line = f.next()
      if line[0]=='#':
        tokens = line.strip()[1:].split('\t')
        self.name, self.lbp_contig, self.rbp_contig = tokens[:3]
        continue
      elif line[0]=='=':
        return
      self.support.append(map(int,line.strip().split('\t')))

    
class AssemblyTAlign(object):
  def __init__(self, talign_fpath):
    self.frags_dict = defaultdict(list)
    with open(talign_fpath, 'rb') as f:
      for line in f:
        if line.startswith('#'):
          continue
        tokens = line.strip().split('\t')
        self.frags_dict[tokens[0]].append(map(int, tokens[2:6]))
      
    self.frags = self.frags_dict.keys()
    
    try:
      self.frags.sort(key=lambda x:x.split('/')[1])
    except KeyError:
      print >>debug_log, "Error: assumes PacBio read names"
      sys.exit(1)
    
  def _get_overlap_right_triangle(self, x1,y1,d1,c1, x2,y2,d2,c2):
    assert c2>=c1
    if d1 < c2-c1 or y2-d2 > y1 or x2-d2 > x1:
      return None
    x12,y12 = min(x1,x2), min(y1,y2)
    d12 = x12+y12-c2
    return x12,y12,d12
    
  def condense_fragments_clustering(self, max_d=None):
    from ambre.analyze.algo import OverlapBPs
    
    cluster = OverlapBPs()
    
    triangles = []
    frag_idx = []
    frag_bp_pos = []
    for idx, frag in enumerate(self.frags):
      if len(self.frags_dict[frag])<2:
        continue
      lbp, rbp, lfp, rfp = zip(*self.frags_dict[frag])
      lbp, rbp, lfp, rfp = na.array(lbp), na.array(rbp), na.array(lfp), na.array(rfp)
      
      orient = na.ones(lbp.size, dtype=na.int)
      orient[lbp>rbp] = -1
      
      # appends +/- as break orientations
      #bp_sequence = zip(rbp[:-1]*orient[:-1]*-1, lbp[1:]*orient[1:], lfp[1:]-rfp[:-1])
      bp_sequence = zip(rbp[:-1]*orient[:-1]*-1, lbp[1:]*orient[1:], lfp[1:]-rfp[:-1]-1)
      
      if not max_d is None:
        bp_sequence = [(a,b,d) for a,b,d in bp_sequence if d<max_d]
      
      triangles.extend(bp_sequence)
      frag_bp_pos.extend(range(len(bp_sequence)))
      frag_idx.extend([idx]*len(bp_sequence))
      
    data = na.array(triangles)

    unique_map = defaultdict(list)
    # get a dataset with unique a,b,d
    # No point in finding overlaps between repeating triangles.
    unique_data = []
    unique_data_c = []
    
    l = sorted(range(len(triangles)), key=triangles.__getitem__)
    
    for data_idx in l:
      if len(unique_data)!=0 and (data[data_idx,:]==unique_data[-1]).all():
        unique_data_c[-1]+=1
      else:
        unique_data.append(data[data_idx,:])
        unique_data_c.append(1)
        
      unique_map[len(unique_data)-1].append(data_idx) 
    
    unique_data = na.array(unique_data)
    print >>debug_log, "Data: ", data.shape
    print >>debug_log, "UniqueData: ", unique_data.shape
    cluster.overlapping_triangles(unique_data)
    ccs = cluster.get_connected_components()
    
    breakpoints = [BreakPoint() for i in xrange(len(ccs))]
    print >>debug_log, "ClusterCount: %d"%len(ccs)
    for bp, bp_unique_idxs in zip(breakpoints, ccs):
      for bp_unique_idx in bp_unique_idxs:
        a,b,d = unique_data[bp_unique_idx,:]
        for bp_idx in unique_map[bp_unique_idx]:
          bp.add_support(frag_idx[bp_idx], frag_bp_pos[bp_idx], a,b,d)
    breakpoints.sort(reverse=True, key=lambda bp:len(bp.support))
    return breakpoints

  def output_breakpoints(self, breakpoints, breakpoint_fpath=None, min_cluster_size=25):
    if breakpoint_fpath is None:
      breakpoint_file = sys.stdout
    else:
      breakpoint_file = open(breakpoint_fpath, 'wb') 
    
    for bp in breakpoints:
      if len(bp.support)<min_cluster_size:
        break
      
      a2,b2,d2 = bp.get_mode()
      
      print >>debug_log, "Size:%d, %d, %d, %d"%(len(bp.support), a2,b2,d2)
      
      bp.to_file(breakpoint_file)
    
    if not breakpoint_fpath is None:
      breakpoint_file.close()
      
    return breakpoints

  
