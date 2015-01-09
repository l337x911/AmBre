import os
import sys
from ambre.utils import reference
from collections import defaultdict
from itertools import izip, count
from operator import attrgetter
import pandas as pd
import re
import numpy as na

class Primer(object):
  def __init__(self, sln_idx=None, region_idx=None, contig=None, pos=None, orient=None, name=None, seq=None):
    self.sln_idx = sln_idx
    self.region_idx = region_idx
    self.contig = contig
    self.pos = pos
    self.orient = orient
    self.name = name
    self.seq = seq
  def print_pos(self):
    assert not self.pos is None
    print self.contig, self.pos
  def info(self):
    print self.name, self.seq, self.orient

class Region(object):
  def __init__(self, contig=None, pos=None, length=None, orient=None):
    self.contig = contig
    self.pos = pos
    self.length = length
    self.orient = orient
  def primer_in(self, p):
    return p.orient==self.orient and p.contig==self.contig and p.pos >= self.pos and p.pos<self.pos+self.length 
  def synthesize_length(self, p, seq_ext=0):
    assert p.orient==self.orient
    assert p.contig==self.contig
    if self.orient:
      return (self.pos+self.length)-p.pos
    else:
      return p.pos-self.pos+seq_ext
  def min_dist_from_region(self, r):
    assert r.contig==self.contig
    x,y = na.meshgrid((r.pos,r.pos+r.length),(self.pos,self.pos+self.length))
    return na.min(na.abs(x-y))  


class NucleicAcidsLinks(object):
  def __init__(self, ref_fasta_fpath, regions_fpath, region_links_fpath, fa_fpath, seq_extension):
    self.regions_fpath = regions_fpath
    self.region_links_fpath = region_links_fpath
    self.fa_fpath = fa_fpath
    self.seq_ext = seq_extension
    self.regions = []
    with open(regions_fpath, 'rb') as f:
      for l in f:
        if l.startswith('#'): continue
        tokens = l.strip().split('\t')
       
        self.regions.append(Region(contig=tokens[0], pos=int(tokens[1]), length=int(tokens[2]), orient='forward'==tokens[3]))

    self.links = []
    with open(self.region_links_fpath, 'rb') as f:
      for l in f:
        if l.startswith('#'): continue
        tokens = l.strip().split('\t')
        self.links.append((int(tokens[0]), int(tokens[1]), tokens[2]))

    self.ref = reference.Reference(ref_fasta_fpath)

    self.slns = defaultdict(set)
    self.slns_full = defaultdict(list)
  def get_id(self, fa):
    tokens = fa[1:10].split('-')
    assert tokens[0]=='sln' 
    return int(tokens[1]), int(tokens[2])

  def get_pos_orient(self, fa):
    tokens = fa.split('_')
    isForward = tokens[-1].startswith('True')
    isReverse = tokens[-1].startswith('False')
    assert isForward or isReverse 
    return tokens[-3], int(tokens[-2]), isForward 

  def unique_regions(self, primer_list):
    unique_primers = dict()
    for p in primer_list:
      if p.region_idx in unique_primers.viewkeys():
        p0 = unique_primers[region_idx]
        assert p.orient==p0.orient
        if p.orient and p.pos<p0.pos: continue
        if (not p.orient) and p.pos>p0.pos: continue

      unique_primers[p.region_idx] = p 
    return unique_primers

  def set_primers(self):
    with open(self.fa_fpath, 'rb') as f:
      fas = [line.strip() for line in f]
    for name, seq in izip(fas[::2], fas[1::2]):  
      sln_idx, region_idx = self.get_id(name)
      contig, pos, orient = self.get_pos_orient(name)
      self.slns[sln_idx].add(region_idx)
      p = Primer(region_idx=region_idx, contig=contig, pos=pos, orient=orient, name=name, seq=seq)
      self.slns_full[sln_idx].append(p)

  def _verify_region(self, p):
    return self.regions[p.region_idx].primer_in(p)

  def get_complete_sln(self, soi=None):
    if not soi is None:
      sln = self.slns_full[soi]
      uniq_sln = self.unique_regions(sln)
      p_verify =  all([self._verify_region(p) for p in uniq_sln.itervalues()])
      if p_verify and len(uniq_sln)==len(self.regions):
        return soi, uniq_sln
      else:
        return None 

    for soi, sln in sorted(self.slns_full.iteritems()):
      uniq_sln = self.unique_regions(sln)
      #print uniq_sln, soi
      p_verify =  all([self._verify_region(p) for p in uniq_sln.itervalues()])
      #print p_verify
      if p_verify and len(uniq_sln)==len(self.regions):
        return soi, uniq_sln
    return None
  
  def product_lengths(self, uniq_sln):
    link_lengths = []
    for i,j,n in self.links:
      pi =self.regions[i].synthesize_length(uniq_sln[i], seq_ext=self.seq_ext)
      pj =self.regions[j].synthesize_length(uniq_sln[j], seq_ext=self.seq_ext)
      link_lengths.append(pi+pj)
    return link_lengths

  def analyze(self):
    self.set_primers()
    soi = self.get_complete_sln()
    if soi is None:
      raise Exception('No complete solution for %s'%self.regions_fpath)
    soi, uniq_sln = soi
    
    product_lengths = self.product_lengths(uniq_sln)
    i,j,n = zip(*self.links)
    
    print "{0:d} - {1}".format(soi, ' '.join(["{0}:{1}".format(name, length) for name, length in zip(n,product_lengths)]))


class PhasingLinks(NucleicAcidsLinks):
  def __init__(self, ref_fasta_fpath, regions_fpath, region_links_fpath, fa_fpath, seq_extension):
    NucleicAcidsLinks.__init__(self, ref_fasta_fpath, regions_fpath, region_links_fpath, fa_fpath, seq_extension)
    self.snp_distance = None
    self.snp_idx = None
    self.inner_idx = None

  def set_distance_to_snp(self):
    t = []
    inner_link_idx = None
    for i,j,n in self.links:
      if not 'snp' in n: continue
      if 'inner' in n:
        inner_link_idx = len(t)
      t.append(set([i,j]))
    snp_regions = list(t[0].intersection(t[1]))
    assert len(snp_regions)==1
    self.snp_idx = snp_regions[0]
    self.inner_idx = [i for i in t[inner_link_idx] if not i==self.snp_idx][0]
    self.snp_distance = self.regions[self.inner_idx].min_dist_from_region(self.regions[self.snp_idx])

  def get_dist_to_snp(self, uniq_sln):
    # region closest to snp
    return self.regions[self.snp_idx].synthesize_length(uniq_sln[self.snp_idx], seq_ext=self.seq_ext)

  def product_lengths(self, uniq_sln):
    new_lengths = []
    for (i,j,n),l in zip(self.links, NucleicAcidsLinks.product_lengths(self, uniq_sln)):
      if not 'snp' in n:
        new_lengths.append(l)
      else:
        new_lengths.append(l+self.snp_distance)
    return new_lengths

  def get_primer_seq(self):
    regions = sorted(self.regions, key=attrgetter('pos'))

  def output_primers_in_sln(self, uniq_sln, row):
    region_labels = ['']*len(self.regions)
    
    if self.regions[self.snp_idx].pos < self.regions[self.inner_idx].pos:
      order = 'F'
    else:
      order = 'R'

    region_labels[self.snp_idx] = 'outer{0}'.format(order)
    region_labels[self.inner_idx] = 'inner'
    region_labels[2] = 'outerDelF'
    region_labels[3] = 'outerDelR'
    #print region_labels, self.inner_idx, self.snp_idx
    names = ['outerF', 'outerDelF', 'inner', 'outerDelR', 'outerR']
    primers = {}
   
     
    for i,n in enumerate(region_labels):
      primers[n] = uniq_sln[i]
      #primers[n] = self.regions[i]
   
    #uniq_sln[self.inner_idx].print_pos() 
    #uniq_sln[self.inner_idx].info()
    #uniq_sln[self.snp_idx].info()
    
    #out_str = []
    #out_str2 = []
    #for n in names:
    #  p = primers[n]
    #  if p is None:
    #    out_str.append('')
    #  else:
    #    out_str.append(p.seq)
        #out_str.append(str(p.orient))
    #    out_str2.append(p.name[:-7])
    #print '\t'.join(out_str)
    if row=='A':
      print "WellPosition\tName\tSequence\tNotes"
    for i,n in enumerate(names):
      if not n in primers: continue
      p = primers[n]  
      print '{0}{1}\t{2}\t{3}\t{4}'.format(row,i+1, p.name[5:-2],p.seq,os.path.basename(self.fa_fpath[:-3]))

  def analyze(self, soi=None):
    self.set_primers()
    self.set_distance_to_snp()
    soi = self.get_complete_sln(soi)
    if soi is None:
      raise Exception('No complete solution for %s'%self.regions_fpath)
    soi, uniq_sln = soi
    
    product_lengths = self.product_lengths(uniq_sln)
    i,j,n = zip(*self.links)
   
    product_length_str =  ' '.join(["{0}:{1:06d}".format(name, length) for name, length in zip(n,product_lengths)])
    print "#{0:02d} - PRIMER_SNP_DIST:{1:04d} SNP_DIST:{2:06d} {3}".format(soi, self.get_dist_to_snp(uniq_sln), self.snp_distance, product_length_str)
    return uniq_sln    

if __name__ == '__main__':
  import sys
  regions_fpath = sys.argv[1]
  region_links_fpath = sys.argv[2]
  fa_fpath = sys.argv[3]
  row = sys.argv[4]
  try:
    soi = int(sys.argv[5])
  except IndexError:
    soi = None
 
  w = PhasingLinks('/media/T02/data/hg/hg18.fa', regions_fpath, region_links_fpath, fa_fpath, 32)
  uniq_sln = w.analyze(soi)
  w.output_primers_in_sln(uniq_sln, row)
   
