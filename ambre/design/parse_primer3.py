'''
#  ambre.design.parse_primer3.py
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

from collections import defaultdict

def parse_info(r):
  primer_idx_pos = r.find('_')
  if primer_idx_pos == -1:
    primer_idx_pos = r.find('=')
  primer_idx = int(r[:primer_idx_pos])
  info = r[primer_idx_pos+1:].split('=')
  # starting position
  if len(info)==1:
    info.insert(0,'')
    info[1] = info[1].split(',')
    info[1] = (int(info[1][0]), int(info[1][1]))
  
  try:
    info[1] = float(info[1])
  except:
    # Let the value be a string
    pass

  
  return primer_idx, info[0], info[1]

def get_primer_dict(primer_file):
  
  primer_data = [i.strip() for i in primer_file if i.startswith('PRIMER')]
  left_primers, right_primers = defaultdict(dict), defaultdict(dict)
  
  for x in primer_data:
    if x.find('RETURNED')!=-1:
      continue
    if x.startswith('PRIMER_LEFT_'):
      idx, info_n, info_v = parse_info(x[12:])
      left_primers[idx][info_n] = info_v
    if x.startswith('PRIMER_RIGHT_'):
      idx, info_n, info_v = parse_info(x[13:])
      right_primers[idx][info_n] = info_v
  return left_primers, right_primers

def get_primer_sublist(fpath):
  import bisect
  l_primers, r_primers = get_primer_dict(fpath)
  pos_l = [(d[''][0],d['SEQUENCE']) for d in l_primers.itervalues()]
  pos_r = [(d[''][0],d['SEQUENCE']) for d in r_primers.itervalues()]
  ls_idx, le_idx = bisect.bisect_left(pos_l, (181600-80000,'')), bisect.bisect_right(pos_l, (181600,''))
  rs_idx, re_idx = bisect.bisect_left(pos_r, (181600+12000,'')), bisect.bisect_right(pos_r, (181600+12000+80000,''))
  pos_l.sort()
  pos_r.sort()
  ls_idx, le_idx = bisect.bisect_left(pos_l, (181600-80000,'')), bisect.bisect_right(pos_l, (181600,''))
  rs_idx, re_idx = bisect.bisect_left(pos_r, (181600+12000,'')), bisect.bisect_right(pos_r, (181600+12000+80000,''))
  return pos_l[ls_idx:le_idx+1], pos_r[rs_idx:re_idx+1]

def get_primer_subset(left_fpath,right_fpath, primer_designs_fpath):
  
  # double the runtime, but who gives a fuck when your computer is fast
  left_primers, x = get_primer_dict(open(left_fpath, 'rb'))
  x, right_primers = get_primer_dict(open(right_fpath, 'rb'))

  left_primer_by_position = dict([(info_dict[''][0], info_dict) for (idx, info_dict) in left_primers.iteritems()])
  right_primer_by_position = dict([(info_dict[''][0], info_dict) for (idx, info_dict) in right_primers.iteritems()])
  RIGHT_OFFSET = 104000
  x = []
  with open(primer_designs_fpath, 'rb') as f:
    for i,line in enumerate(f):
      tokens = line.strip().split('\t')
      l,r = [], []
      left_positions = map(int, tokens[2].split(','))
      left_positions.sort()
      for j in left_positions:
        l.append(left_primer_by_position[j])
      
      right_positions = map(int, tokens[3].split(','))
      right_positions.sort()
      for j in right_positions:
        r.append(right_primer_by_position[j-RIGHT_OFFSET])
      x.append((l,r))
  return x

def print_primer3_input(x):
  for l,r in  x:
    print "Solution"
    for i in l:
      print "SEQUENCE_PRIMER=%s"%i
    for i in r:
      print "SEQUENCE_PRIMER_REVCOMP=%s"%i

def plot_primer3_penalty_distr(pl_dict, pr_dict, penalty_bins=None):
  import numpy as na
  from matplotlib import pyplot as plt

  if penalty_bins is None:
    penalty_bins = na.arange(0.2,1.7, 0.2)
    
  cm = plt.get_cmap('hot')
  bin_width = 2000
  
  colors = [cm(i) for i in na.linspace(0,1,penalty_bins.size)]
  
  # pos and penalty
  x = [(i[''][0], i['PENALTY']) for k,i in pl_dict.iteritems() if i.has_key('')] + [(i[''][0], i['PENALTY']) for k,i in pr_dict.iteritems() if i.has_key('')]
  x.sort()
  
  bins = na.arange(na.around(x[0][0],-5), na.around(x[-1][0],-5), bin_width)
  
  pos,pen = map(na.array, zip(*x))

  prev = 0 
  for idx, i in enumerate(penalty_bins):
    p = pos[pen<i]  
    print i, len(p)      
    h,b = na.histogram(p, bins=bins)
    plt.bar(left=b[:-1], height=h-prev, width=bin_width, bottom=prev, color=colors[idx], label="%.2f"%i)
    prev = na.copy(h)
  plt.xlim(b[0],b[-1])
  plt.title('Distribution of Primer3 primers separated by penalty')
  plt.legend()
  plt.show()
  
  
def get_penalty_sum():
  x = get_primer_subset('/home/vdeshpan/PAMP/data/primer3_chr9del_left.out', 
            '/home/vdeshpan/PAMP/data/primer3_chr9del_right.out', 
            '/home/anand/Projects/PAMP/best_pamp.out')
  OUT_DIR = '/home/anand/Projects/PAMP/best_pamp_sln_'
  for sln_idx, (l,r) in enumerate(x):
    
    penalties = sum([i['PENALTY'] for i in l])+sum([i['PENALTY'] for i in r])
    with open("%s%d.fa"%(OUT_DIR, sln_idx), 'wb') as f:
      
      t = "solution_%d %.2f"%(sln_idx, penalties)
      for idx, i in enumerate(l):
        print >>f, ">LP_%s%s"%(idx,t)
        print >>f, i['SEQUENCE']
      for idx, i in enumerate(r):
        print >>f, ">RP_%s%s"%(idx, t)
        print >>f, i['SEQUENCE']    

def plot_penalty():
  with open('/home/anand/Projects/PAMP/pamp_w/pamp_A6/a4.primer3.out', 'rb') as f:
    pl, pr = get_primer_dict(f)
  plot_primer3_penalty_distr(pl, pr)

if __name__ == '__main__':
  plot_penalty()
  #get_penalty_sum()
    