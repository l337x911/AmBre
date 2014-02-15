import sys
from collections import defaultdict
from itertools import izip
import pandas as pd

# assume mrna

regions_fpath = '/home/anand/Projects/ambre/gbm/gbm-rna-screen_cDNA04.txt'
region_links_fpath = '/home/anand/Projects/ambre/gbm/gbm-rna-screen_cDNA04_links.txt'
fa_fpath = '/home/anand/Projects/ambre/gbm/gbm-rna-screen_cDNA04.b4.fa' 
region_count = 19
soi = 1

def get_id(fa):
  tokens = fa[1:10].split('-')
  assert tokens[0]=='sln' 
  return int(tokens[1]), int(tokens[2])

def get_pos_orient(fa):
  tokens = fa[10:].split('_')
  isForward = tokens[-1].startswith('True')
  isReverse = tokens[-1].startswith('False')
  assert isForward or isReverse 
  return int(tokens[-2]), isForward 

def unique_regions(primer_list):
  unique_primers = dict()
  for primer_idx, pos, orient, name, seq in primer_list:
    if primer_idx in unique_primers.viewkeys():
      old_pos, old_orient, old_name, old_seq = unique_primers[primer_idx]
      assert orient==old_orient
      if orient and pos<old_pos: continue
      if (not orient) and pos>old_pos: continue

    unique_primers[primer_idx] = pos, orient, name, seq
  return unique_primers
  #k,v = unique_primers.items()        
  #return sorted(zip(k,*zip(*v)))

def _get_length(regions, uniq_sln, i):
  c,a,b,o = regions[i]
  try:
    pos,orient,name, seq = uniq_sln[i]
  except KeyError:
    return b-a

  assert o==orient
  if orient:
    return b-pos
  else:
    return pos-a

def get_product_lengths(regions_fpath, region_links_fpath, uniq_sln):
  regions = []
  links = []

  with open(regions_fpath, 'rb') as f:
    for line in f:
      if line.startswith('#'): continue
      tokens = line.strip().split('\t')
      contig, a, l, orient = (tokens[0], int(tokens[1]), int(tokens[2]), tokens[3]=='forward')
      regions.append((contig, a, a+l, orient))

  with open(region_links_fpath, 'rb') as f:
    for line in f:
      if line.startswith('#'): continue
      tokens = line.strip().split('\t')
      idxi,idxj,name = int(tokens[0]),int(tokens[1]),tokens[2]
      links.append((idxi,idxj,name))
  ret = []
  for i,j,n in links:
    l_i,l_j = _get_length(regions, uniq_sln, i), _get_length(regions, uniq_sln, j)
    ret.append((n, l_i+l_j))
  return ret    

with open(fa_fpath, 'rb') as f:
  fas = [line.strip() for line in f]

t = zip(fas[::2], fas[1::2])

slns = defaultdict(set)
slns_full = defaultdict(list)
for name, seq in izip(fas[::2], fas[1::2]):
  sln_idx, primer_idx = get_id(name)
  pos, orient = get_pos_orient(name)
  slns[sln_idx].add(primer_idx)

  slns_full[sln_idx].append((primer_idx, pos, orient, name, seq))

print "# Slns with a primer in each region"
for sln_idx in sorted(slns.iterkeys()):
  region_with_primers = slns[sln_idx]
  if len(region_with_primers) != region_count: continue
  print "%02d: "%sln_idx, sorted(list(region_with_primers))

uniq_sln = unique_regions(slns_full[soi])
print "# Product lengths from %d"%soi
for n,l in get_product_lengths(regions_fpath, region_links_fpath, uniq_sln):
  print "%04d\t%s"%(l,n)

print "# Primers in %d"%soi
for i in sorted(uniq_sln.keys()):
  pos,orient,name, seq = uniq_sln[i]
  print "%s\n%s"%(name,seq)

