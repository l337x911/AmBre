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
    self.primers = None
    self.d = 6500
  def parse_primers(self, fpath):
    self.primers = []
    with open(fpath, 'rb') as f:
      # Quench header    
      f.next()
      for p in f:
        tokens = p.strip().split('\t')
        #print tokens
        contig, norm_pos, orient, sequence, penalty = tokens[0], int(tokens[1]), tokens[2] == 'True', tokens[3], float(tokens[4])
        self.primers.append((contig, norm_pos, orient))
    self.primers.sort()
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
  def load(self, primers_fpath, breakpoints_fpath):
    self.parse_primers(primers_fpath)
    self.parse_breakpoints(breakpoints_fpath)
        
  def check_positions(self):
    # assumes that breakpoints are on the same contig
    
    c, primer_pos, orient = zip(*self.primers)

    primer_pos, orient = na.array(primer_pos), na.array(orient, bool)    
    
    coverage = {True:'left', False:'right'}

    valid_count = defaultdict(int)
    sample_count = defaultdict(int)
    #print "Position Checks", primer_pos, orient
    
    samples = [(name, pos_for, True) for name in self.samples_for.keys() for pos_for in self.samples_for[name]]
    samples += [(name, pos_rev, False) for name in self.samples_rev.keys() for pos_rev in self.samples_rev[name]]

    sample_dist = (defaultdict(list), defaultdict(list),)
    for name, bp_pos, bp_orient in sorted(samples):
      m = orient == bp_orient
      #print "Sample", bp_orient, bp_pos
      
      x = primer_pos[m] - bp_pos

      if bp_orient:
        x *= -1

      #print x
      primer_indices = na.logical_and(x > 0, x < self.d)
      x_valid = primer_pos[m][primer_indices]
      if len(x_valid) > 0:
        valid_count[name] += 1
      
      x_invalid = list(x[na.logical_and(x > 0, x < 2 * self.d)])
      x_invalid.sort()
      try:
        dist_to_closest_primer = x_invalid[0]
      except:
        dist_to_closest_primer = 2 * self.d
      
      if bp_orient:
        sample_dist[0][name].append(dist_to_closest_primer)
      else:
        sample_dist[1][name].append(dist_to_closest_primer)
        
      
      sample_count[name] += 1
      print name, coverage[bp_orient], bp_pos, ",".join(["%d" % i for i in x_valid]), list(na.nonzero(primer_indices)[0] + 1)
      
    for i in sample_count.keys():
      # Checks if sample's breakpoints are covered by some primer.
      # note that only the boundaries are checked.
      print "%s\t%02d/%02d" % (i, valid_count[i], sample_count[i])
    
    print "Sample\tMinAmpliconSize\tMaxAmpliconSize"
    for name in sample_dist[0].keys():
      for_l, rev_l = sample_dist[0][name], sample_dist[1][name]
      for_l.sort()
      rev_l.sort()
      if len(for_l) == 0 or len(rev_l) == 0:
        print "%s\t0\t0"
      else:
        print "%s\t%06d\t%06d" % (name, for_l[0] + rev_l[0], for_l[-1] + rev_l[-1])
    
    # Primer Design Coverage 
    for_primers = primer_pos[orient]
    rev_primers = primer_pos[na.logical_not(orient)]
    
    f_m = na.min(for_primers)
    r_m = na.min(rev_primers)
    for_basepairs = na.zeros((na.max(for_primers) - f_m + self.d + 1), dtype=na.bool)
    rev_basepairs = na.zeros((na.max(rev_primers) - r_m + self.d + 1), dtype=na.bool) 
    for i in for_primers:
      i -= f_m 
      for_basepairs[i:i + self.d] = True
    r_m -= self.d
    for j in rev_primers:
      j -= r_m
      rev_basepairs[j - self.d:j] = True

    
    print "ForCov %d/%d ; RevCov %d/%d" % (na.sum(for_basepairs), len(for_basepairs), na.sum(rev_basepairs), len(rev_basepairs))
    print "ForCov %.4f ; RevCov %.4f" % (na.sum(for_basepairs) / float(len(for_basepairs)), na.sum(rev_basepairs) / float(len(rev_basepairs)))
    print "ForCost %d ; RevCost %d" % (len(for_basepairs) - na.sum(for_basepairs), len(rev_basepairs) - na.sum(rev_basepairs))
    
  def plot_same_contig(self, show=True):
  
    from matplotlib import pyplot as plt
    
    points = [pos for name in self.samples_for.keys() for pos in self.samples_for[name]]
    points += [pos for name in self.samples_rev.keys() for pos in self.samples_rev[name]]
    c, pos, orient = zip(*self.primers)
    pos = na.array(pos)
    orient = na.array(orient)

    points += list(pos)
    points.sort()
    rounding = 10000
    plot_start = points[0] - rounding
    plot_end = points[-1] + rounding
    
    plot_start = na.floor(plot_start / rounding) * rounding
    plot_end = na.ceil(plot_end / rounding) * rounding
    
    sample_names = list(set(self.samples_for.keys()).union(set(self.samples_rev.keys())))
    sample_names.sort()

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
    #print pos, orient
    for_pos = pos[orient > 0]
    rev_pos = pos[orient < 1]
    for_heights = na.ones(len(for_pos)) * (len(sample_names) + 1)
    rev_heights = na.ones(len(rev_pos)) * (len(sample_names) + 1)

    plt.bar(for_pos, for_heights, width=self.d, bottom= -1, color='blue', alpha=0.2)
    plt.bar(rev_pos - self.d, rev_heights, width=self.d, bottom= -1, color='blue', alpha=0.2)

    plt.bar(for_pos, for_heights, width=30, bottom= -1, color='blue', alpha=0.9)
    plt.bar(rev_pos - 30, rev_heights, width=30, bottom= -1, color='blue', alpha=0.9)
    #print pos[orient>0]
    
    breakpoints = []
    for idx, sample_name in enumerate(sample_names):
      idx = len(sample_names) - idx - 1
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
    plt.ylim(-1, len(sample_names))
    plt.xticks(na.arange(plot_start, plot_end + 2, 20000), [number_formatter(x / 1000) for x in na.arange(plot_start, plot_end + 2, 20000)])
    plt.yticks(na.arange(len(sample_names)), [name for name in sample_names[::-1]])
    plt.xlabel('Genomic Position in Kb')
    plt.ylabel('Sample')
    plt.title('Regions')
    if show:
      plt.show()

if __name__ == '__main__':
  import argparse
  
  parser = argparse.ArgumentParser(description="""Displays expected amplicon lengths given estimated breakpoints. 
  If simple deletions, then plots primers against estimated breakpoints for simple deletions.""")
  parser.add_argument('primers', type=str, nargs=1, help='Single primer design (from workflow output). First line is header, followed by primers.')
  parser.add_argument('breakpoints', type=str, nargs=1, help="""TAB-delimited deletions file. A row is of the form:
  <sample_id> <fasta_seqid> <left bp range start>  <left bp range end> <right bp range start>  <right bp range end> <left bp>  <right bp>
where left and right bps can be empty (-) """)
   
  args = parser.parse_args()
  
  w = PlotAmBreSolution()
  w.load(primers_fpath=args.primers[0], breakpoints_fpath=args.breakpoints[0])
  w.check_positions()
  try:
    import matplotlib
    w.plot_same_contig()
  except:
    pass
  
