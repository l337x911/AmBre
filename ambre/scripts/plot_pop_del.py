from collections import defaultdict
from matplotlib import gridspec, pyplot as plt, ticker

import numpy as na

all_panel = '/media/T01/data/1000genomes/1k.ALL.panel'
wgs_lce = '/media/T01/data/1000genomes/ALL.wgs.merged_5_del_call_sets_bps.20101123.sv_dels.low_coverage.sites.lce.vcf'

cmap = plt.get_cmap('Set1')

class Population(object):
  def __init__(self, name, c, idx, size):
    self.name = name
    self.c = cmap(c)
    self.i = idx
    self.size = size

def uniqify_deletions(dels):
  pct_overlap = 0.80
  srtd_dels = sorted(dels, key=lambda x:(x[0], x[0]+x[1],x[2]))
  print "Sortded dels", [(x[0],x[0]+x[1]) for x in srtd_dels]
  new_dels = [srtd_dels[0],]
  for i in xrange(1,len(srtd_dels)):
    c_p,c_l,c_s = srtd_dels[i]
    p_p,p_l,p_s = new_dels[-1]
    print "#", c_p,c_l,p_p,p_l
    print float(p_p+p_l-c_p)/max(c_l,p_l)
    if float(p_p+p_l-c_p)/max(c_l,p_l) > pct_overlap:
      t_p = min(p_p, c_p)
      t_l = max(p_p+p_l, c_p+c_l)-t_p
      t_s = p_s
      t_s.update(c_s)
      new_dels[-1] = (t_p, t_l, t_s)
    else:
      print "New del"
      new_dels.append(srtd_dels[i])
  return new_dels

sample_pop = {}
pops_c = defaultdict(int)

with open(all_panel, 'rb') as f:
  for line in f:
    tokens = line.strip().split('\t')
    sample_pop[tokens[0]] = tokens[1]
    pops_c[tokens[1]] += 1

pops = sorted(pops_c.keys())
pops = [Population(n, c, i, pops_c[n]) for i, n,c  in zip(range(len(pops)), pops, na.linspace(0,0.8,num=len(pops)))]
y_pops = dict([(p.name, p) for p in pops])

y_shift = na.linspace(-0.4,0.4,num=len(pops))

fig = plt.figure()

gs = gridspec.GridSpec(1,3,width_ratios=[10,2,1])
gs.update(left=0.05, right=0.95, wspace=0.02)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax3 = plt.subplot(gs[2])

def lw_func(x):
  return x/10.

y_idx = 1
delsfound = []

with open(wgs_lce, 'rb') as f:
  for line in f:
    if line.startswith('#'): continue
    tokens = line.strip().split('\t')
    pos1 = int(tokens[1])
    l = len(tokens[3])
    if l<500: continue
    samples = [s for s in tokens[7].split(';') if s.startswith('SAMPLES=')][0][8:].split(',')
    delsfound.append((pos1,l,set(samples)))

delsfound = uniqify_deletions(delsfound)

for pos1,l,samples in delsfound:
    counts = defaultdict(int)
    for s in samples:
      try:
        counts[sample_pop[s]] += 1
      except KeyError:
        print "%s is in VCF but not in ALL panel"%s

    print "chr1:%d-%d"%(pos1, pos1+l), sum(counts.values())
    for n,i in counts.iteritems():
      p = y_pops[n]
      ax1.plot((pos1,pos1+l),(y_idx+y_shift[p.i], y_idx+y_shift[p.i]), c=p.c, lw=lw_func(i), marker=None, solid_capstyle='butt')
      if i>0:
        print "\t%s\t%d\t%d\t%.2f"%(n,i,p.size, float(i)/p.size)
    y_idx += 1

text_orient = 'left'
for y, p in zip(y_shift, pops):
  ax2.plot((0,1),(y,y), c=p.c, label=p.name, lw=lw_func(pops_c[p.name]))
  ax2.text(.1,y, p.name, horizontalalignment=text_orient, verticalalignment='center')
 
step = 5
for y, i in zip(y_shift, na.linspace(step, max(pops_c.values())+step, num=len(pops))):
  ax3.plot((0,.4),(y,y), c='k', lw=lw_func(i))
  ax3.text(.5,y, "%d"%i, horizontalalignment=text_orient, verticalalignment='center')

hg_formatter = ticker.ScalarFormatter()
hg_formatter.set_scientific(False)
ax1.xaxis.set_major_formatter(hg_formatter)
ax2.axis('off')
ax3.axis('off')

ax2.set_ylim(-0.5,0.5)
ax3.set_ylim(-0.5,0.5)
ax1.set_ylabel('Deletions')
ax1.set_xlabel('CHR1 pos (grch37)')

ax2.set_title('Populations\n (# in line width)')
ax3.set_title('Number of individ\n in populations')
plt.show()

