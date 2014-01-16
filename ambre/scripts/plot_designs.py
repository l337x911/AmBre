'''
#  ambre.design.plot_designs.py
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
import sys, os
from collections import defaultdict
import numpy as na
import bisect
import locale


def cell_line_from_fpath(fpath):
  return os.path.basename(fpath).split('.')[0]

def number_formatter(x):
  x = '%d' % x
  return ','.join([x[:(len(x) % 3)]] + [''.join(x) for x in zip(x[-3::-3], x[-2::-3], x[-1::-3])][::-1])

    
class PlotAmBreSolution(object):
  def __init__(self):
    self.samples_for = None
    self.samples_rev = None
    self.inputs = []
    self.input_names = []
    self.d = 6500
  def parse_inputs(self, fpath):
    regions = []
    with open(fpath, 'rb') as f:
      for p in f:
        tokens = p.strip().split('\t')
        #print tokens
        contig, norm_pos, length, orient = tokens[0], int(tokens[1]), int(tokens[2]), tokens[3] == 'forward'
        regions.append((contig, norm_pos, length, orient))
    self.inputs.append(regions)
    self.input_names.append(os.path.basename(fpath))
    #print self.primers
  def parse_breakpoints(self, fpath):
    self.samples_for = defaultdict(list)
    self.samples_rev = defaultdict(list)
    with open(fpath, 'rb') as f:
      for s in f:
        if s[0] == '#':
          continue
        
        tokens = s.strip().split('\t')
        sample_name = tokens[0]
        pos = map(int, tokens[2:])
        if pos[4] > 0:
          self.samples_for[sample_name].append(pos[4])
        else:
          self.samples_for[sample_name].append(pos[0])
          self.samples_for[sample_name].append(pos[1])
        if pos[5] > 0:
          self.samples_rev[sample_name].append(pos[5])
        else:
          self.samples_rev[sample_name].append(pos[2])
          self.samples_rev[sample_name].append(pos[3])

  def load(self, input_fpaths, breakpoints_fpath):
    for input_fpath in input_fpaths:
      self.parse_inputs(input_fpath)
    self.parse_breakpoints(breakpoints_fpath)
        
  def plot_same_contig(self, show=True):
  
    from matplotlib import pyplot as plt
    
    points = [pos for name in self.samples_for.keys() for pos in self.samples_for[name]]
    points += [pos for name in self.samples_rev.keys() for pos in self.samples_rev[name]]

    pos = []
    for regions in self.inputs:
      c, p, l, o = zip(*regions)
      pos.extend(pos)
      pos.extend(na.array(p)+na.array(l))

    points += list(pos)
    points.sort()
    rounding = 10000
    plot_start = points[0] - rounding
    plot_end = points[-1] + rounding
    
    plot_start = na.floor(plot_start / rounding) * rounding
    plot_end = na.ceil(plot_end / rounding) * rounding
    
    sample_names = list(set(self.samples_for.keys()).union(set(self.samples_rev.keys())))
    sample_names.sort()
    sample_names = sample_names[:2]+sample_names[3:]
    sample_names.append('Detroit562')    
    

    for name in sample_names:
      try:
        self.samples_for[name].sort()
      except:
        pass
      try:
        self.samples_rev[name].sort()
      except:
        pass

    plt.figure()

    ori_color = {True:'b', False:'r'}
    for idx, regions in enumerate(self.inputs):
      idx += 1
      for c,p,l,o in regions:
	plt.plot((p,p+l),(-idx,-idx), color=ori_color[o], lw=5)
    
    breakpoints = []
    for idx, sample_name in enumerate(sample_names):
      idx = len(sample_names) - idx - 1
      print idx, sample_name
      plt.plot((plot_start, plot_end), (idx, idx), color='k', alpha=0.5, lw=2)      
      plt.plot((self.samples_for[sample_name][0], self.samples_rev[sample_name][-1]), (idx, idx), linestyle='--', color='g', lw=3)
    
      if len(self.samples_for[sample_name]) == 2:
        plt.plot(self.samples_for[sample_name], (idx, idx), color='g', marker='o', lw=5, markersize=7)
      elif len(self.samples_for[sample_name]) == 1:
        breakpoints.append((self.samples_for[sample_name][0], idx))

      if len(self.samples_rev[sample_name]) == 2:
        plt.plot(self.samples_rev[sample_name], (idx, idx), color='g', marker='o', lw=5, markersize=7)
      elif len(self.samples_rev[sample_name]) == 1:
        breakpoints.append((self.samples_rev[sample_name][0], idx))
    if len(breakpoints) > 0:
      bp_x, bp_y = zip(*breakpoints)
      plt.scatter(bp_x, bp_y, s=30, c='g', marker='o')

    plt.xlim(plot_start, plot_end)
    plt.ylim(-1-len(self.inputs), len(sample_names))
    plt.xticks(na.arange(plot_start, plot_end + 2, 20000), [number_formatter(x / 1000) for x in na.arange(plot_start, plot_end + 2, 20000)])
    plt.yticks(na.arange(-len(self.inputs), len(sample_names)), self.input_names[::-1] + [name for name in sample_names[::-1]])
    
    plt.xlabel('Genomic Position in Kb')
    plt.ylabel('Sample')
    plt.title('Regions')
    if show:
      plt.show()

if __name__ == '__main__':
  import argparse
  
  parser = argparse.ArgumentParser(description="""Displays expected amplicon lengths given estimated breakpoints. 
  If simple deletions, then plots primers against estimated breakpoints for simple deletions.""")
  parser.add_argument('breakpoints', type=str, nargs=1, help="""TAB-delimited deletions file. A row is of the form:
  <sample_id> <fasta_seqid> <left bp range start>  <left bp range end> <right bp range start>  <right bp range end> <left bp>  <right bp>
where left and right bps can be empty (-) """) 
  parser.add_argument('inputs', type=str, nargs='+', help='AmBre input regions')

  args = parser.parse_args()
  
  w = PlotAmBreSolution()
  w.load(input_fpaths=args.inputs, breakpoints_fpath=args.breakpoints[0])
  try:
    import matplotlib
    w.plot_same_contig()
  except ImportError:
    print "No matplotlib"
    pass
  
