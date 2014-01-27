'''
#  ambre.design.workflow.py
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
from ambre.config import CONFIG

from subprocess import Popen, PIPE
from collections import namedtuple, defaultdict
import os
import sys
import argparse
import numpy as na
import shlex
from operator import attrgetter
import tempfile
from itertools import izip
import ambre.design.parse_primer3 as parse_primer3
import ambre.design.sa_cost as sa_cost
from ambre.utils import reference

import time

# Might be useless if all the programs
# are linux only.

#match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes      qStarts  tStarts
#        match   match           count   bases   count   bases           name            size    start   end     name            size    start   end     count
BLAT_ALIGNMENT_H="match mismatch repmatch n_s qgap qgap_bases sgap tgap_bases strand qseqid qlen qstart qend sseqid slen sstart send block_count block_sizes q_starts t_starts"

BlatAlignment=namedtuple('BlatAlignment', BLAT_ALIGNMENT_H)
PSolution=namedtuple('PSolution', "param cost time primers")

class Primer3FileManager(object):
  def __init__(self, temptag, regions_on_seq, seq, primers_per_kbp):
    self.temptag = temptag
    self.regions = regions_on_seq
    self.seq = seq
    self.primers_per_kbp = primers_per_kbp
    self.name = None
  def exists(self):
    for fpath in self.iter_input_fpaths():
      if not os.path.isfile(fpath):
        return False
    return True
  def remove_in(self):
    for fpath in iter_input_fpaths():
      os.remove(fpath)
  def remove_out(self):
    for fpath in iter_output_fpaths():
      os.remove(fpath)

  def _prepare_input(self, fpath, seq, forward_flag):
    if forward_flag:
      ok_reg = ",,%d,5"%(len(seq)-5)
    else:
      ok_reg = "0,5,,"

    with open(fpath, 'wb') as primer3_in:
      print >>primer3_in, "PRIMER_NUM_RETURN=%d"%((len(seq)*self.primers_per_kbp)/1000)

      print >>primer3_in, "SEQUENCE_ID=pamp_workflow"
      print >>primer3_in, "SEQUENCE_TEMPLATE=%s"%seq
      print >>primer3_in, "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=%s"%(ok_reg)
      print >>primer3_in, "="

  def prepare_inputs(self):
    for (p,l,fflag), fpath in izip(self.regions,self.iter_input_fpaths()):
      self._prepare_input(fpath, self.seq[p:p+l], fflag)
    self.name = ",".join(self.iter_input_fpaths())
  def iter_input_fpaths(self):
    for i in xrange(len(self.regions)):
      yield '%s.%04d.primer3'%(self.temptag, i)
  def iter_output_fpaths(self):
    for i in xrange(len(self.regions)):
      yield '%s.%04d.primer3.out'%(self.temptag, i)
  def iter_output(self):
    for region, fpath in izip( self.regions,self.iter_output_fpaths()):
      yield region, fpath 

class PrimerDesignWorkflow(object):
  '''
  classdocs
  '''

  def __init__(self, primer3_path, primer3_param, aligner, multiplx, d=6500, rho=0.2):
    '''
    Constructor
    '''
    self.primer3_path = primer3_path
    self.primer3_param = primer3_param
  
    assert aligner.endswith('blat')
    self.aligner = aligner
    self.multiplx = multiplx
    
    self.alignments_dict = None
    self.cross_amp = None
    self.alignments_in_regions = None
    self.regions = None
    self.regions_on_seq = None
    self.solutions = []
    self.forward_dict = {}
    self.reverse_dict = {}
    self.map_pos_to_primer_idx = {}
    self.is_forward_dict = {True:'forward', False:'reverse'}
    self.d = d
    self.rho = rho
  def get_original_region_position(self, position):
    a,b,c = zip(*self.regions_on_seq)
    a = na.array(a)
    b = a+na.array(b)
    # slow way to get the original position with original region.
    region_idx = na.nonzero(na.logical_and(position>=a,position<b))[0][0]
    contig, original_region_a, original_region_b, original_region_orient = self.regions[region_idx]
    norm_pos = position-self.regions_on_seq[region_idx][0]
    return contig, original_region_a+norm_pos, original_region_orient

  def get_original_region(self, position):
    a,b,c = zip(*self.regions_on_seq)
    a = na.array(a)
    b = na.array(b)+a
    region_idx = na.nonzero(na.logical_and(position>=a,position<b))[0][0]
    return region_idx, self.regions[region_idx]
  
  def set_alignments_dict(self, *args, **kwargs):
    if self.aligner.endswith('blat'):
      self.set_alignments_dict_from_blat(*args, **kwargs)
  
  def set_alignments_dict_from_blat(self, blat_out_file, max_align=5, min_align_len=15, init_tolerance=1):
    alignments_dict = defaultdict(list)
    for i in blat_out_file:
      try:
        match = BlatAlignment._make(i.split('\t'))
      except:
        print "Blast alignment Failure %s"%i
        continue
      
      if int(match.send)-int(match.sstart)<min_align_len:
        continue
      
      is_forward = match.strand=='+'
      # Ensures left is a suffix match
      if match.qseqid.startswith("LEFT") and (not is_forward or int(match.qlen) - int(match.qend)>init_tolerance):
        continue
      
      # Ensures right is a prefix match
      elif match.qseqid.startswith("RIGHT") and (is_forward or int(match.qstart)>init_tolerance):
        continue
      primer_pos = int(match.qseqid.split(':')[-1])
      
      # Each alignment is stored as (contig, alignment start, is_forward)
      alignments_dict[primer_pos].append((match.sseqid, int(match.sstart), is_forward))
  
    # Remove primers with too many alignments    
    for k,v in alignments_dict.items():
      if len(v)>max_align:
        alignments_dict.pop(k)
    self.alignments_dict = alignments_dict
    print "#Number of Primers after Blat: %d"%len(self.alignments_dict)

    
  def check_cross_amplification(self, max_dist=20000):
    
    alignments = []
    for primer_pos, v in self.alignments_dict.iteritems():
      primer_poss = [primer_pos]*len(v)
      c,a,s = zip(*v)
      alignments.extend(zip(c,a,s,primer_poss))
    
    alignments.sort()
    
    chrs, pos, orient, primer_pos = zip(*alignments)
    
    primer_pos = na.array(primer_pos)
    chrs_x, chrs_y = na.meshgrid(chrs, chrs)
    chrs_eq = chrs_x==chrs_y
    chrs_x, chrs_y = None, None
    
    orient_x,orient_y = na.meshgrid(orient, orient)
    orient_fr = orient_x==na.logical_not(orient_y)
    orient_x, orient_y = None, None
    
    pos = na.array(pos, dtype=na.int32)
    pos_x,pos_y = na.meshgrid(pos, pos)
    
    # Filter by minimum distance.
    # And if orientation is forward/reverse
    # Amps is [Alignments i X Alignments j ]
    amps = na.logical_and(na.logical_and(na.abs(pos_x-pos_y)<max_dist, 
                      na.logical_and(orient_fr, pos_x<pos_y)), # Forward reverse if one position is less than the other and is forward, while the other is reverse              
              chrs_eq)
    
    # Convert each alignment back to their original primer_pos
    # [ Primers i X Primers j]
    
    self.cross_amp = set()
    primers_i, primers_j = na.nonzero(amps)
    if len(primers_i)==0:
      print "#Number of Cross Amplifications 0"
      return
    for i,j in zip(primer_pos[primers_i], primer_pos[primers_j]):
      # Note that if multiple alignments of a single primer generate a unique amplicon,
      # then i==j and is still included as "cross amplification"
      self.cross_amp.add((i,j))
    
    print "#Number of Cross Amplifications %d"%len(self.cross_amp)
    
  def set_regions(self, regions_fpath):
    self.regions = []
    with open(regions_fpath, 'rb') as f:
      for l in f:
        if l.startswith('#'): continue
        tokens = l.strip().split('\t')
         
        c,a,b, s = tokens[0], int(tokens[1]), int(tokens[2]), tokens[3]
        try:
          assert s=='forward' or s=='reverse'
        except AssertionError:
          print >>sys.stderr, "Regions file orientation is not forward or reverse"
        if s=='forward':
          self.regions.append((c, a, b, True))
        else:
          self.regions.append((c, a, b, False))
          
  def set_sequence(self):
    ref = reference.Reference(CONFIG.param['reference_fpath'])
    
    self.seq = ""
    self.regions_on_seq = []
    for c,a,b,s in self.regions:
      
      start = len(self.seq)
      self.seq += ref.get_sequence(c,a,a+b)
      
      region_length = len(self.seq)-start
      assert region_length == b
      self.regions_on_seq.append((start, region_length, s))
    
    ref.close()


  def _get_primers(self, primer3_out, 
                  seq_offset=0,
                  forward_offset=0, 
                  reverse_offset=0, 
                  max_penalty=1.5):
    
    forward_dict, reverse_dict = parse_primer3.get_primer_dict(primer3_out)
      
    for primer_idx, info_dict in forward_dict.iteritems():
      if float(info_dict['PENALTY'])>max_penalty: continue

      info_dict[''] = (info_dict[''][0]+seq_offset, info_dict[''][1])
      self.forward_dict[primer_idx+forward_offset] = info_dict
      self.map_pos_to_primer_idx[info_dict[''][0]] = (self.forward_dict, primer_idx+forward_offset)

    for primer_idx, info_dict in reverse_dict.iteritems():
      if float(info_dict['PENALTY'])>max_penalty: continue
      
      info_dict[''] = (info_dict[''][0]+seq_offset, info_dict[''][1])
      self.reverse_dict[primer_idx+reverse_offset] = info_dict
      self.map_pos_to_primer_idx[info_dict[''][0]] = (self.reverse_dict, primer_idx+reverse_offset)

  def get_primers(self, fm, max_penalty=1.5):
    self.forward_dict = {}
    self.reverse_dict = {}
    self.map_pos_to_primer_idx = {}
 
    for (start_on_seq,length_on_seq,isforward), primer3_out_fpath in fm.iter_output():
      with open(primer3_out_fpath, 'rb') as f:
        self._get_primers(f, seq_offset=start_on_seq, forward_offset=len(self.forward_dict), reverse_offset=len(self.reverse_dict), max_penalty=max_penalty)     

  def add_cross_amp_to_multiplx(self, edges_fpath):
    BAD_SCORES_STR = "-100\t-100\t-100"
    with open(edges_fpath, 'a') as f:
      for (s_i,s_j) in self.cross_amp:
        print >>f, "%d\t%d\t%s"%(s_i, s_j, BAD_SCORES_STR)
  
  def parse_pamp_output(self, output_str):
    slns = []
    for line in output_str.split('\n')[:-1]:
      sln = PSolution._make(line.split('\t'))
      
      f_cost = int(sln.cost)
      
      f_primers = [map(int, [l for l in reg.split(',') if len(l)>0]) for reg in sln.primers.split(':')]
      f_sln = PSolution(sln.param, cost=f_cost, time=sln.time, primers=f_primers)
      
      slns.append(f_sln)

    self.solutions = sorted(slns, key=attrgetter('cost'))
  
  def run(self, regions_fpath, 
        temp_tag=None, delete_flag=False, 
        primers_per_kbp=75, max_cross_amp_dist=20000,
        pamp_off=False,
        max_primer_penalty=1.5,
        max_alignment_count=10,
        min_alignment_len=18,
        pamp_max_iterations= 100000000,pamp_repeats= 3,
        pamp_t_ms= (-1,-(10**-1),-(10**-2),-(10**-3)), pamp_t_bs=(10**4,10**5,10**6)):
    cur_dir = os.path.abspath(os.curdir)
    self.set_regions(regions_fpath)
    self.set_sequence()
    
    if temp_tag is None:
      out_dir=CONFIG.dir['default_temp']
      temp_tag = tempfile.mktemp(prefix='', dir=out_dir)
    elif os.path.isdir(temp_tag):
      out_dir=temp_tag
      temp_tag = tempfile.mktemp(prefix='', dir=out_dir)
    
    # Get Primer3 Output
    primer3_fm, primer3_out = None, None
    aligner_in, aligner_out_fpath = None, None    
    multiplx_in, multiplx_out_fpath = None, None
    pamp_out = None
    
    try:
      pre_primer3_time = time.time()
      primer3_fm = Primer3FileManager(temp_tag, self.regions_on_seq, self.seq, primers_per_kbp)
      if not primer3_fm.exists():

        primer3_fm.prepare_inputs()

        os.chdir(self.primer3_path)

        for primer3_infpath, primer3_outfpath in izip(primer3_fm.iter_input_fpaths(),primer3_fm.iter_output_fpaths()):
  
          primer3_cmd = "%s -p3_settings_file=%s %s"%(os.path.join(self.primer3_path, 'src', 'primer3_core'), self.primer3_param, primer3_infpath)
          print primer3_cmd
          with open(primer3_outfpath, 'wb') as primer3_out:
            p_primer3 = Popen(shlex.split(primer3_cmd), stdout=primer3_out, stderr=PIPE)
        
            (primer3_out_str, primer3_err) = p_primer3.communicate()
        
        os.chdir(cur_dir)
      
      self.get_primers(primer3_fm, max_penalty=max_primer_penalty)


      print "#Stage1 Primer3Time: %.4f"%(time.time()-pre_primer3_time)
       
      pre_align_time = time.time()
      aligner_out_fpath = '%s.align.out'%temp_tag
      if not (os.path.isfile('%s.align.fa'%temp_tag) and os.path.isfile(aligner_out_fpath)):
        
        aligner_in = open('%s.align.fa'%temp_tag, 'wb')
        primers_fa = [">LEFT:%d\n%s"%(info_dict[''][0],info_dict['SEQUENCE']) for info_dict in self.forward_dict.itervalues()]
        primers_fa += [">RIGHT:%d\n%s"%(info_dict[''][0],info_dict['SEQUENCE']) for info_dict in self.reverse_dict.itervalues()]
        
        print >>aligner_in, "\n".join(primers_fa)
        aligner_in.close()

        if self.aligner.endswith('blat'):
          align_cmd = '%s -t=dna -q=dna -stepSize=2 -tileSize=8 -repMatch=1048576 -minScore=15 -noHead %s %s %s'%(self.aligner, CONFIG.param['reference_fpath'], aligner_in.name, aligner_out_fpath)
        print len(primers_fa), align_cmd
        p_align = Popen(shlex.split(align_cmd), stdout=PIPE, stderr=PIPE)
        
        (aligner_out, aligner_err) = p_align.communicate()
        
      with open(aligner_out_fpath, 'rb') as aligner_out_f:
        self.set_alignments_dict(aligner_out_f, max_align=max_alignment_count, min_align_len=min_alignment_len)
      print "#Stage2 AlignTime: %.4f"%(time.time()-pre_align_time)
      
      pre_multiplx_time = time.time()
      multiplx_out_fpath = '%s.multiplx.out'%temp_tag
      if not (os.path.isfile('%s.multiplx'%temp_tag) and os.path.isfile(multiplx_out_fpath)):
        multiplx_in = open('%s.multiplx'%temp_tag, 'wb') 
        
        for primer_pos in self.alignments_dict.iterkeys():
          primer_dict, primer_idx = self.map_pos_to_primer_idx[primer_pos]
          print >>multiplx_in, "%d\t%s\tN"%(primer_pos, primer_dict[primer_idx]['SEQUENCE'])
        multiplx_in.close()
        
        os.chdir(os.path.dirname(self.multiplx))
        
        multiplx_cmd = '%s -primers %s -calcscores 123 -savescores %s'%(self.multiplx, multiplx_in.name, multiplx_out_fpath)
        p_multiplx = Popen(shlex.split(multiplx_cmd), stdout=PIPE, stderr=PIPE)
        
        (multiplx_out_str, multiplx_err) = p_multiplx.communicate();
        os.chdir(cur_dir)
      print "#Stage2 MultiPlxTime: %.4f"%(time.time()-pre_multiplx_time)
      
      pre_crossamp_time = time.time()
      self.check_cross_amplification(max_dist=max_cross_amp_dist)
      self.add_cross_amp_to_multiplx(multiplx_out_fpath)
      print "#Stage2 CrossAmpTime: %.4f"%(time.time()-pre_crossamp_time)
      
      regions_pamp = ''.join(["%d\t%d\t%s\n"%(a,b,self.is_forward_dict[s]) for (a,b,s) in self.regions_on_seq])
      
      pre_pamp_time = time.time()
      primer_graph = sa_cost.get_primer_graph(multiplx_out_fpath, regions_pamp, d=self.d, rho=self.rho)
      
      if not pamp_off:
        pamp_out = open('%s.sa'%temp_tag, 'wb')
      
        sa_cost.multi_vary_temperature_schedule(graph=primer_graph, 
                        max_iterations=pamp_max_iterations, repeats=pamp_repeats,
                        t_ms=pamp_t_ms, t_bs=pamp_t_bs,
                        output=pamp_out)
        pamp_out.close()
        print "#Stage3 PAMPTime: %.4f"%(time.time()-pre_pamp_time)
      
        with open(pamp_out.name, 'rb') as pamp_out:  
          self.parse_pamp_output(pamp_out.read())
        
    finally:
      if not primer3_out is None:
        primer3_out.close()
        print primer3_out.name
      if not aligner_in is None:
        aligner_in.close()
        print aligner_in.name
      if not aligner_out_fpath is None:
        print aligner_out_fpath
      if not multiplx_in is None:
        multiplx_in.close()
        print multiplx_in.name
      if not multiplx_out_fpath is None:
        print multiplx_out_fpath
      
      if not pamp_out is None:
        pamp_out.close()
        print pamp_out.name
        
      if(delete_flag):
        primer3_fm.remove_in()
        primer3_fm.remove_out()
        os.remove(aligner_in.name)
        os.remove(aligner_out_fpath)
        os.remove(multiplx_in.name)
        os.remove(multiplx_out_fpath)
        os.remove(pamp_out.name)
  def set_compatible_primers_to_sln(self, regions_fpath, temp_tag):
    assert not temp_tag is None
    regions_pamp = ''.join(["%d\t%d\t%s\n"%(a,b,self.is_forward_dict[s]) for (a,b,s) in self.regions_on_seq])
    multiplx_out_fpath = '%s.multiplx.out'%temp_tag
    primer_graph = sa_cost.get_primer_graph(multiplx_out_fpath, regions_pamp, d=self.d, rho=self.rho)
    slns = []
    print >>sys.stderr, na.array(primer_graph.primer_graph_mat, dtype=int)
    for i in xrange(len(primer_graph.combined_primers)):
      primer_i_dict, primer_i_idx = self.map_pos_to_primer_idx[primer_graph.combined_primers[i]] 
      contig_i, norm_pos_i, orientation_i = self.get_original_region_position(primer_i_dict[primer_i_idx][''][0]) 
      print >>sys.stderr, contig_i, norm_pos_i, orientation_i
    for (i,j), v in na.ndenumerate(primer_graph.primer_graph_mat):
      if (not v) and i<j:
        primer_i_dict, primer_i_idx = self.map_pos_to_primer_idx[primer_graph.combined_primers[i]]
        contig_i, norm_pos_i, orientation_i = self.get_original_region_position(primer_i_dict[primer_i_idx][''][0])
        primer_j_dict, primer_j_idx = self.map_pos_to_primer_idx[primer_graph.combined_primers[j]]
        contig_j, norm_pos_j, orientation_j = self.get_original_region_position(primer_j_dict[primer_j_idx][''][0])
        if orientation_i==orientation_j:
          continue
        slns.append(PSolution("%d,%d"%(i,j),
                    #cost=primer_i_dict[primer_i_idx]['PENALTY']+primer_j_dict[primer_j_idx]['PENALTY'], 
                    #cost=abs(norm_pos_i-norm_pos_j),
                    cost=abs(primer_graph.combined_primers[i]-primer_graph.combined_primers[j]),
                    time="0", 
                    primers=[[primer_graph.combined_primers[i]], [primer_graph.combined_primers[j]]]))
        
    self.solutions = sorted(slns, key=attrgetter('cost'))
    
  def print_solutions(self, out_fpath=None, top=5):
    j = min(top, len(self.solutions))
    if out_fpath is None:
      out = sys.stdout
    else:
      out = open(out_fpath, 'wb')
    
    for i in range(j):
      sln = self.solutions[i]
      print >>out, "%d\t%s\t%s"%(sln.cost, sln.time, sln.param)
      for reg in sln.primers:
        for p in reg:
          primer_dict, primer_idx = self.map_pos_to_primer_idx[p]
          contig, norm_pos, orientation = self.get_original_region_position(primer_dict[primer_idx][''][0])
          print >>out, "%s\t%d\t%s\t%s\t%.4f"%(contig, norm_pos, orientation, primer_dict[primer_idx]['SEQUENCE'], primer_dict[primer_idx]['PENALTY'])
        
      print >>out, "="
    
    if not out_fpath is None:
      out.close()
    else:
      out.flush()
 
  def print_fasta_solutions(self, out_fpath=None, top=5):
    n = min(top, len(self.solutions))
    
    if out_fpath is None:
      out = sys.stdout
    else:
      out = open(out_fpath, 'wb')

    for i in xrange(n):
      sln = self.solutions[i]
      for reg in sln.primers:
        for p in reg:
          primer_dict, primer_idx = self.map_pos_to_primer_idx[p]
          contig, norm_pos, orientation = self.get_original_region_position(primer_dict[primer_idx][''][0])
          region_idx, regions = self.get_original_region(primer_dict[primer_idx][''][0])
           
          print >>out, ">sln-%02d-%02d_%s_%d_%s\t%.4f\n%s"%(i, region_idx, contig, norm_pos, orientation, primer_dict[primer_idx]['PENALTY'], primer_dict[primer_idx]['SEQUENCE'])
    if not out_fpath is None:
      out.close()
    else:
      out.flush()  

  def check(self, regions_fpath, temp_tag, out_fpath, fasta_flag=False):
    pamp_out_fpath = '%s.sa'%temp_tag

    self.set_regions(regions_fpath)
    self.set_sequence()

    primer3_fm = Primer3FileManager(temp_tag, self.regions_on_seq, self.seq, 0)
    self.get_primers(primer3_fm)

    with open(pamp_out_fpath, 'rb') as pamp_out:  
      self.parse_pamp_output(pamp_out.read())
    if fasta_flag:
      self.print_fasta_solutions(out_fpath=out_fpath, top=15)
    else:
      self.print_solutions(out_fpath=out_fpath, top=15)
    try: 
      import matplotlib
      
      self.validate()
    except ImportError:
      pass

  def validate(self):
    from matplotlib import font_manager,pyplot as plt
    
    # Two figures
    # Coverage of primer3 primers if Not deleted
    # 
    # Top 2 Solution Coverages
    fig = plt.figure()
    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,
                    wspace=None, hspace=0.27)
    num_solutions = min(3, len(self.solutions))
    num_plot = num_solutions+1
    ax_primer3 = fig.add_subplot(num_plot,1,1)
    x_range = (-2,len(self.seq)+2)
    prop = font_manager.FontProperties(size=10) 
    
    
    f_pos = na.array([i[''][0] for i in self.forward_dict.itervalues()])
    r_pos = na.array([i[''][0] for i in self.reverse_dict.itervalues()])
    
    colors = ['b','r']*(len(self.regions)/2+1)
    
    bin_width = len(self.seq)/150
    for idx, ((c,a,b,s), (s_a,s_b,s_s)) in enumerate(zip(self.regions, self.regions_on_seq)):
      
      tag = "%s %s:%d-%d"%(self.is_forward_dict[s], c,a,a+b)
      if s:
        pos = f_pos
        n_pos = r_pos
      else:
        pos = r_pos
        n_pos = f_pos
        
      ax_primer3.hist(pos[na.logical_and(pos>s_a, pos<(s_b+s_a))], bins=(b/bin_width), color=colors[idx], label=tag)
      try:
        print pos.shape, "NPos", n_pos.shape, "SA", s_a, "SB", s_b
        wrong_primers = pos[na.logical_and(n_pos>s_a, n_pos<s_b)] 
        if len(wrong_primers)>0:
          print "Region %s has opposing oriented primers, %s"%(tag, ",".join(['%d'%i for i in wrong_primers]))
      except:
        pass
    
    
    ax_primer3.set_xlim(*x_range)
    ax_primer3.legend(loc='upper right', prop=prop)
    ax_primer3.set_title('histogram of filtered primers')
    ax_primer3.set_ylabel('number of primers')
    ax_primer3.xaxis.set_major_locator(plt.NullLocator())
    colors = {True:'b', False:'r'}
    
    for i in range(num_solutions):
      sln = self.solutions[i]
      ax_sln = fig.add_subplot(num_plot, 1, i+2)
      ax_sln.set_title('Solution %d %d'%(i+1, sln.cost))
      
      for (a,b,s), reg in zip(self.regions_on_seq, sln.primers):
        
        left = na.array(reg)
        
        direction = 1 if s else -1
        
        ax_sln.bar(left, [100*direction]*len(left), bottom=[0]*len(left), width=direction*self.d, color=colors[s], alpha=0.5)
    
        ax_sln.vlines([a,a+b],[0,0],[100,100], lw=3)
      ax_sln.set_xlim(*x_range)
      ax_sln.yaxis.set_major_locator(plt.NullLocator())
      if i<num_solutions-1:
        ax_sln.xaxis.set_major_locator(plt.NullLocator())
        
    ax_sln.set_xlabel('DNA positions covered')  
    plt.show()

  
def test_cross_amp(regions_fpath, temp_tag):
  w = PrimerDesignWorkflow(primer3_path=CONFIG.dir['primer3'], 
          primer3_param=CONFIG.param['primer3_long'], 
          aligner=CONFIG.bin['aligner'], 
          multiplx=CONFIG.bin['multiplx'])
  
  w.set_regions(regions_fpath)
  w.set_sequence()
  
  primer3_fm = Primer3FileManager(temp_tag, w.regions_on_seq, w.seq, 0) 
  w.get_primers(primer3_fm)
 
  with open('%s.align.out'%temp_tag, 'rb') as blast_out_f:
    
    w.set_alignments_dict(blast_out_f,
                          max_align=int(CONFIG.param['design_max_alignments']),
                          min_align_len=int(CONFIG.param['design_3end_len_alignment']))
  
    
  w.check_cross_amplification(max_dist=int(CONFIG.param['design_max_cross_amp_dist']))
  
  multiplx_out_fpath = '%s.multiplx.out'%temp_tag
  if (os.path.isfile('%s.multiplx'%temp_tag) and os.path.isfile(multiplx_out_fpath)):
    regions_pamp = ''.join(["%d\t%d\t%s\n"%(a,b,w.is_forward_dict[s]) for (a,b,s) in w.regions_on_seq])
    
    primer_graph = sa_cost.get_primer_graph(multiplx_out_fpath, regions_pamp, d=w.d, rho=w.rho)
    print "#Incompatible primer density: ", primer_graph.get_dimerization_density()
    print "#Primers without dimers: ", primer_graph.get_number_of_primers_without_dimers() 
  if len(w.cross_amp)>1:
    cross_amp_index =  na.array(list(w.cross_amp))
    if len(w.cross_amp)>30:
      cross_amp_index = cross_amp_index[na.random.random_integers(0,len(w.cross_amp)-1,30)]
    for amp in cross_amp_index:
      print "CrossAmp:", amp  
      for p in amp:
        primer_dict, primer_idx = w.map_pos_to_primer_idx[p]
        contig, norm_pos, orientation = w.get_original_region_position(primer_dict[primer_idx][''][0])
        print "%s\t%d\t%s\t%s"%(contig, norm_pos, orientation, primer_dict[primer_idx]['SEQUENCE'])
  
       
  
  
def check(regions_fpath, temp_tag, out_fpath, fasta_flag):
  w = PrimerDesignWorkflow(primer3_path=CONFIG.dir['primer3'], 
          primer3_param=CONFIG.param['primer3_long'], 
          aligner=CONFIG.bin['aligner'], 
          multiplx=CONFIG.bin['multiplx'],
          d=int(CONFIG.param['design_primer_spacing']),
          rho=float(CONFIG.param['design_primer_density']))

  w.check(regions_fpath, temp_tag, out_fpath, fasta_flag=fasta_flag)

def ambre_run(regions_fpath, temp_tag=None, out_fpath=None,  fasta_flag=False):
  w = PrimerDesignWorkflow(primer3_path=CONFIG.dir['primer3'], 
          primer3_param=CONFIG.param['primer3_long'], 
          aligner=CONFIG.bin['aligner'], 
          multiplx=CONFIG.bin['multiplx'],
          d=int(CONFIG.param['design_primer_spacing']),
          rho=float(CONFIG.param['design_primer_density']))

  w.run(regions_fpath,
        delete_flag=(CONFIG.param['cleanup_flag']=="True"),
        temp_tag=temp_tag,
        primers_per_kbp=int(CONFIG.param['design_primer3_primers_per_kbp']),
        max_primer_penalty=float(CONFIG.param['design_max_primer3_penalty']),
        max_cross_amp_dist=int(CONFIG.param['design_max_cross_amp_dist']),
        max_alignment_count=int(CONFIG.param['design_max_alignments']),
        min_alignment_len=int(CONFIG.param['design_3end_len_alignment']),
        pamp_max_iterations= int(CONFIG.param['design_sa_max_iterations']),
        pamp_repeats= int(CONFIG.param['design_sa_repetitions']),
        pamp_t_ms= map(float,CONFIG.param['design_sa_ms'].split(',')),
        pamp_t_bs=map(float,CONFIG.param['design_sa_bs'].split(',')))
  if fasta_flag:
    w.print_fasta_solutions(out_fpath=out_fpath)
  else:
    w.print_solutions(out_fpath=out_fpath)
  try:
    import matplotlib
    w.validate()
  except ImportError:
    pass

def main(): 
  parser = argparse.ArgumentParser(prog='ambre-design.py', description='Select compatible primers covering reference region.')
  parser.add_argument('-c, --check', dest='check', action='store_const', const=check, default=None, help='checks if the temptag solution is valid')
  parser.add_argument('-f, --fasta-out', dest='ofasta', action='store_const', const=True, default=False, help='output in fasta')
  parser.add_argument('-a, --check-align', dest='check_align', action='store_const', const=test_cross_amp, default=None, help='Checks for cross amplifications')
  
  parser.add_argument('reference', type=str, nargs=1, help='Fasta with multiple sequence entries.')
  parser.add_argument('regions', type=str, nargs=1, help='TAB-delimited regions file. A row is of the form:<fasta_seqid>  <start_position>  <region_length>  <primer_orientation>')
  
  parser.add_argument('temptag', type=str, nargs='?', default=None, help='Prefix for the run id / Directory to store Temps of a run id')
  parser.add_argument('-o, --primers-out', type=str, nargs=1, default=[None], dest="primer_fpath", help='Output primers to tab-delimited file. Default is to output to STDOUT')
  parser.add_argument('--config', type=str, nargs=1, default=None, dest="config_fpath", help='Update parameters in default config file with new config file.')
   
  args = parser.parse_args()
  
  CONFIG.param['reference_fpath'] = args.reference[0]
  if not args.config_fpath is None:
    CONFIG.update(args.config_fpath[0])
  
  if not args.check is None:
    try: 
      assert not args.temptag is None
      args.check(args.regions[0], args.temptag, args.primer_fpath[0], fasta_flag=args.ofasta)
    except AssertionError:
      print "No temp id prefix specified"
      sys.exit()
  elif args.check_align is not None:
    try: 
      assert not args.temptag is None
      args.check_align(args.regions[0], args.temptag)
    except AssertionError:
      print "No temp id prefix specified"
      sys.exit()
    
  else:
    ambre_run(args.regions[0], args.temptag, args.primer_fpath[0], fasta_flag=args.ofasta)

if __name__=='__main__':
  main()

