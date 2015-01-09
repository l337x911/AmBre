'''
#  ambre.design.sa_analysis.py
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

from matplotlib import pyplot as plt, cm
import numpy as na
import sys, os
from collections import defaultdict

def get_parameters(fpath):
  # slope,intercept, repeat
  tokens = fpath.split('.')
  
  return map(int, tokens[-3:])

class SAStats(object):
  def __init__(self, fpath=None):
    self.fpath = fpath
    self.runtime = None
    self.avg_runtime_per_iter = None
    self.iterations = []
    self.min_cost = []
    self.curr_cost = []
    
    if not fpath is None and os.path.isfile(fpath):
      self.load()
  def load(self):
    self.runtime = 0
    self.avg_runtime_per_iter = 0
    
    
    with open(self.fpath, 'rb') as f:
      rows = []
      for line in f:
        if line.startswith('Iterations'):
          tokens = line.strip().split('\t')
          rows.append((int(tokens[1]), float(tokens[2]), int(tokens[3]), int(tokens[4])))
        if line.startswith('TerminationClock'):
          self.runtime = float(line.strip().split(':')[-1])

    if len(rows) == 0:
      return
    self.iterations, self.avg_runtime_per_iter, self.min_cost, self.curr_cost = zip(*rows)
    self.iterations, self.avg_runtime_per_iter, self.min_cost, self.curr_cost = na.array(self.iterations), na.array(self.avg_runtime_per_iter), na.array(self.min_cost), na.array(self.curr_cost)  
    self.avg_runtime_per_iter = na.average(self.avg_runtime_per_iter)

  def average_stats(self, stats_objects):
    
    self.avg_runtime_per_iter = 0
    self.runtime = 0
    
    # Assumes same iteration lengths
    self.iterations = stats_objects[0].iterations
    self.min_cost = na.zeros(len(stats_objects[0].min_cost))
    self.curr_cost = na.zeros(len(stats_objects[0].curr_cost))
    
    for stat_obj in stats_objects:
      #assert not ((stat_obj.iterations - self.iterations) > 0).any()
      self.min_cost += stat_obj.min_cost
      self.curr_cost += stat_obj.curr_cost
      self.avg_runtime_per_iter += stat_obj.avg_runtime_per_iter
      self.runtime += stat_obj.runtime
      
    self.min_cost /= len(stats_objects)
    self.curr_cost /= len(stats_objects)
    self.avg_runtime_per_iter /= len(stats_objects)
    self.runtime /= len(stats_objects)

def get_colordict(val_set, val_cm):
  ms = list(val_set)
  ms.sort()
  colors_m = {}
  for i, m in enumerate(ms):
    colors_m[m] = val_cm(int(64 + ((256 - 64) * (i)) / len(ms)))
  return colors_m

def simulated_annealing_curves(sa_runs):
  
  fig = plt.figure()
  
  
  ms, bs = set(), set()
  repeat_params = defaultdict(list)
  for fpath in sa_runs:
    m, b, i = get_parameters(fpath)
    ms.add(m)
    bs.add(b)
    repeat_params[(m, b)].append(SAStats(fpath))
  
  m_colormap = cm.get_cmap('jet')
  b_colormap = cm.get_cmap('Greys')
  
  colors_m, colors_b = get_colordict(ms, m_colormap), get_colordict(bs, b_colormap)
  
  axMinCost = fig.add_subplot(131)
  axCurrCost = fig.add_subplot(132)
  axLegend = fig.add_subplot(133, adjustable='box', aspect=0.2)
  counter = 0
  
  ks = repeat_params.keys()
  ks.sort()
  
  l_repeat_params = [(k, repeat_params[k]) for k in ks]
  for (m, b), stats_objs in l_repeat_params:
    
    avg_stat = SAStats()
    avg_stat.average_stats(stats_objs)
    print m, b, avg_stat.runtime, avg_stat.avg_runtime_per_iter
    
    mlw, blw = 1.5, 5
    
    alpha = 0.8
    
    axMinCost.plot(avg_stat.iterations, avg_stat.min_cost, color=colors_b[b], lw=blw, alpha=alpha)
    axMinCost.plot(avg_stat.iterations, avg_stat.min_cost, color=colors_m[m], lw=mlw, alpha=alpha)
    
    axCurrCost.plot(avg_stat.iterations, avg_stat.curr_cost, color=colors_b[b], lw=blw, alpha=alpha)
    axCurrCost.plot(avg_stat.iterations, avg_stat.curr_cost, color=colors_m[m], lw=mlw, alpha=alpha)
    
    label = "m=10^%d b=10^%d" % (m, b)
    idx = len(repeat_params) - counter
    axLegend.plot((0, 1), (idx, idx), color=colors_b[b], lw=blw, alpha=alpha)
    axLegend.plot((0, 1), (idx, idx), color=colors_m[m], lw=mlw, alpha=alpha)
    axLegend.text(0.5, idx + .2, label, fontsize=12, horizontalalignment='center')
    counter += 1
          
  axLegend.set_ylim(0, len(repeat_params) + 1)
  axLegend.set_xticks([])
  axLegend.set_yticks([])
  
  
  axMinCost.set_xlabel('Iterations')
  axCurrCost.set_xlabel('Iterations')
  axMinCost.set_ylabel('Cost')
  axCurrCost.set_ylabel('Cost')
  axMinCost.set_title('Minimum State Cost')
  axCurrCost.set_title('Current State Cost')
  axMinCost.set_xscale('log')
  axCurrCost.set_xscale('log')
       
  x1, x2 = axMinCost.get_xlim()
  x2 *= 2
  axMinCost.set_xlim((x1, x2))
  axCurrCost.set_xlim((x1, x2))

  plt.show()
  
if __name__ == '__main__':
  import os
  import argparse
  
  parser = argparse.ArgumentParser(description='Plotting simulated annealing for 2-parameter schedules from a single AmBre run')
  parser.add_argument('temptag', type=str, nargs=1, default=None, help='Prefix for the run id')
  args = parser.parse_args()
  dir_basename = os.path.dirname(args.temptag[0])
  basename = "%s.sa"%os.path.basename(args.temptag[0])
  sa_runs = [os.path.join(dir_basename, f) for f in os.listdir(dir_basename) if f!=basename and f.startswith(basename)]
  simulated_annealing_curves(sa_runs)
