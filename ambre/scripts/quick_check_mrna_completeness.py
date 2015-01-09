import sys
from ambre.utils import reference
from collections import defaultdict
from itertools import izip, count
import pandas as pd
import re
import numpy as na

# assume mrna
#regions_fpath = '/home/anand/Projects/ambre/gbm/gbm-rna-screen_cDNA04.1.txt'
#regions_fpath = '/home/anand/Projects/ambre/gbm/gbm-rna-screen_cDNA06.txt'
#regions_fpath = '/home/anand/Projects/ambre/na12878p/hic_het_01.txt'
#region_links_fpath = '/home/anand/Projects/ambre/na12878p/hic_het_01_links.txt'
#region_links_fpath = '/home/anand/Projects/ambre/gbm/gbm-rna-screen_cDNA04_links.txt'
#region_links_fpath = '/home/anand/Projects/ambre/gbm/gbm-rna-screen_cDNA06_links.txt'
#fa_fpath = '/home/anand/Projects/ambre/gbm/gbm-rna-screen_cDNA04.b4.fa' 
#fa_fpath = '/home/anand/Projects/ambre/na12878p/hic_het_01.a01.clean.fa' 
#ref_fasta_fpath = '/media/T02/data/hg/hg18.fa'

class NucleicAcidsLinks(object):
  def __init__(self, ref, regions_fpath, region_links_fpath, fa_fpath, seq_extension, soi=1):
    self.regions_fpath = regions_fpath
    self.region_links_fpath = region_links_fpath
    self.fa_fpath = fa_fpath
    self.seq_extension = seq_extension
    self.soi = soi
    with open(regions_fpath, 'rb') as f:
      self.region_count = sum([1 for l in f if not l.startswith('#')])

    self.ref = reference.Reference(ref_fasta_fpath)

  def get_id(self, fa):
    tokens = fa[1:10].split('-')
    assert tokens[0]=='sln' 
    return int(tokens[1]), int(tokens[2])

  def get_pos_orient(self, fa):
    tokens = fa[10:].split('_')
    isForward = tokens[-1].startswith('True')
    isReverse = tokens[-1].startswith('False')
    assert isForward or isReverse 
    return int(tokens[-2]), isForward 

  def unique_regions(self, primer_list):
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

  def _get_length(self, regions, uniq_sln, i):
    c,a,b,o = regions[i]
    try:
      pos,orient,name, seq = uniq_sln[i]
    except KeyError:
      raise
    return b-a
    assert o==orient
    if orient:
      return b-pos
    else:
      return len(seq)+pos-a

def _get_seq(regions, uniq_sln, i, extend=0):
  c,a,b,o = regions[i]
  try:
    pos,orient,name, seq = uniq_sln[i]
  except KeyError:
    raise
    return b-a

  m = re.search("(gi\|.*\|)[^|]*$",name)
  contig = m.group(1)

  assert o==orient
  if orient:
    #return b-pos
    reg_seq = ref.get_contig(contig)
    assert reg_seq[pos:pos+len(seq)]==seq
    reg_seq = reg_seq[pos-extend:b]  
    try:
      assert seq in reg_seq
    except:
      print reg_seq
      print seq
      raise
    p = reg_seq.find(seq)
    reg_seq = "{0}>{1}>{2}".format(reg_seq[:p],reg_seq[p:p+len(seq)],reg_seq[p+len(seq):])
  else:
    #return pos-a
    reg_seq = ref.get_contig(contig)[a:pos+len(seq)+extend]
    s_seq = seq.translate(reference.BASE_COMPLEMENT)[::-1] 
    assert s_seq in reg_seq
    p = reg_seq.find(s_seq)
    reg_seq = "{0}<{1}<{2}".format(reg_seq[:p],reg_seq[p:p+len(seq)],reg_seq[p+len(seq):])
  return reg_seq

def validate_primers_in_regions(regions, uniq_sln):
  for idx, region in enumerate(regions):
    c,a,b,o = region
    pos,orient,name,seq = uniq_sln[idx]
    assert o==orient
    if pos<a and pos>b:
      print "Invalid primer: Outside exon", c,a,b,pos,orient,name
    elif pos>b and orient:
      print "Invalid primer: Forward past exon", c,a,b,pos,orient,name
    elif pos<a and not orient:
      print "Invalid primer: Reverse before exon", c,a,b,pos,orient,name
    elif orient:
      print "Valid forward primer:", b-pos, b-(pos+len(seq)), name
    elif not orient:
      print "Valid reverse primer:", pos-a, pos+len(seq)-a, name
    
      
# ASSUMES 1 primer per REGION! and they are stored in order!

def get_products(regions_fpath, region_links_fpath, uniq_sln, extend=0):
  regions = []
  links = []

  with open(regions_fpath, 'rb') as f:
    lines = [line.strip() for line in f if not line.startswith('##')]

  for l1,l2 in izip(lines[::2], lines[1::2]):
    assert l1.startswith('#')
    tokens = l2.split('\t')
    tokens2 = l1.split(' ')

    contig, a, b, orient = (tokens[0], int(tokens2[3]), int(tokens2[4]), tokens[3]=='forward')
    regions.append((contig, a, b, orient))

  with open(region_links_fpath, 'rb') as f:
    for line in f:
      if line.startswith('#'): continue
      tokens = line.strip().split('\t')
      idxi,idxj,name = int(tokens[0]),int(tokens[1]),tokens[2]
      links.append((idxi,idxj,name))

  validate_primers_in_regions(regions, uniq_sln)
  ret = []
  for i,j,n in links:
    l_i,l_j = _get_length(regions, uniq_sln, i), _get_length(regions, uniq_sln, j)
    s_i,s_j = _get_seq(regions, uniq_sln, i, extend=extend), _get_seq(regions, uniq_sln, j, extend=extend)
    ret.append((n, l_i+l_j, s_i+s_j))
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

region_histogram = na.zeros(region_count, dtype=na.float)
print "# Slns with a primer in each region %d"%(max(map(len,slns.values())))
for sln_idx in sorted(slns.iterkeys()):
  region_with_primers = slns[sln_idx]
  region_histogram[list(region_with_primers)] += 1
  if len(region_with_primers) != region_count: continue
  print "%02d: "%sln_idx, sorted(list(region_with_primers))

print ' | '.join(["%d,%.2f"%(i,v) for i,v in enumerate(region_histogram/float(len(slns)))])

uniq_sln = unique_regions(slns_full[soi])
print "# Products from %d"%soi
#for n,l,s in get_products(regions_fpath, region_links_fpath, uniq_sln, extend=seq_extension):
#  try:
#    assert len(s.replace('>','').replace('<',''))==l+2*seq_extension
#  except:
#    raise
#  print "%04d\t%s\t%s"%(l,n,s)

print "# Primers in %d"%soi
for i in sorted(uniq_sln.keys()):
  pos,orient,name, seq = uniq_sln[i]
  print "%s\n%s"%(name,seq)

