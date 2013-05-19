'''
Created on May 13, 2013

@author: anand
'''
import argparse
from subprocess import Popen,PIPE
import os, sys
import time
import shlex
import tempfile
import numpy as na
from ambre.config import CONFIG
from ambre.utils import reference
from collections import defaultdict
from ambre.design.workflow import PrimerDesignWorkflow

class PrimerBinarySearchWorkflow(PrimerDesignWorkflow):
  '''
  classdocs
  '''
  def __init__(self, primer3_path, primer3_param, aligner):
    PrimerDesignWorkflow.__init__(self,primer3_path, primer3_param, aligner, '')

    self.ref = reference.Reference(CONFIG.param['reference_fpath']) 
    self.pairs_per_iter = None
    self.length = None
  def primer_pairs_even_spaced(self, length, number_of_pairs):
    from bisect import bisect
    
    positions, sorted_idx = zip(*sorted([(0,None), (length, None)]+[((self.reverse_dict[idx][''][0]+self.forward_dict[idx][''][0])/2,idx) for idx in self.pairs_idx]))
    
    spacing = length/(number_of_pairs+1)
    uniform_pos =  na.arange(spacing,length,step=spacing)[:number_of_pairs]
    bisect_idx = na.array([bisect(positions,p) for p in uniform_pos])
    positions = na.array(positions)
    sorted_idx = na.array(sorted_idx)
    
    assert bisect_idx.size == number_of_pairs
    
    # dist from left side
    left_dist = uniform_pos - positions[bisect_idx-1]
    
    # dist from right side
    right_dist = positions[bisect_idx] - uniform_pos
    
    primer_chose = na.copy(bisect_idx)
    primer_chose[left_dist<right_dist] -= 1
    
    primer_idx_chose = sorted_idx[primer_chose]
    # primer_idx_chose may have NANs
    return primer_idx_chose
 
  def recurse_2annotate(self, n, iteration=1, prefix=''):
    if n==1:
      return [(iteration, prefix),]
    m = (n-1)/2
    left = self.recurse_2annotate(m,iteration+1, prefix+'L')
    right = self.recurse_2annotate(m,iteration+1, prefix+'R')
    return left + [(iteration, prefix),] + right
 
  def recurse_3annotate(self, n, iteration=1, prefix=''):
    if n==1:
      return [(iteration, prefix),]
    m = (n-2)/3
    left = self.recurse_2annotate(m,iteration+1, prefix+'L')
    middle = self.recurse_2annotate(m,iteration+1, prefix+'M')
    right = self.recurse_2annotate(m,iteration+1, prefix+'R')
    return left + [(iteration, prefix),] + middle + [(iteration, prefix),] + right

  def check_cross_amplification(self, max_dist):
    pop_set = set()
    for pos in self.pairs_pos:
      p1,p2 = pos

      c1,a1,s1 = zip(*self.alignments_dict[p1])
      c2,a2,s2 = zip(*self.alignments_dict[p2])
     
      chrs_x, chrs_y = na.meshgrid(c1, c2)
      chrs_eq = chrs_x==chrs_y
    
      orient_x,orient_y = na.meshgrid(s1, s2)
      orient_fr = orient_x==na.logical_not(orient_y)
    
      pos_x,pos_y = na.meshgrid(a1, a2)
    
      # Filter by minimum distance.
      # And if orientation is forward/reverse
      # Amps is [Alignments i X Alignments j ]
      amps = na.logical_and(na.logical_and(na.abs(pos_x-pos_y)<max_dist, 
                      na.logical_and(orient_fr, pos_x<pos_y)), # Forward reverse if one position is less than the other and is forward, while the other is reverse              
              chrs_eq)
      number_of_amps = na.sum(amps)
      #assert number_of_amps!=0
      if number_of_amps>1:
        pop_set.add(pos)
        
    self.pairs_pos.difference_update(pop_set)

  def run(self, contig, start, end, number_of_pairs, 
          temp_tag=None, delete_flag=False, out_fpath=None,
          max_cross_amp_dist=20000,
          max_primer_penalty=1.2,
          max_alignment_count=10,
          min_alignment_len=18,
          primers_per_kbp=50):
    
    valid_2counts = set(list(na.power(2, na.arange(1,10))))
    valid_3counts = set(list(na.power(3, na.arange(1,10))))
    assert number_of_pairs+1 in valid_2counts or number_of_pairs+1 in valid_3counts
    
    if number_of_pairs+1 in valid_2counts:
      self.pairs_per_iter = (2,int(na.log2(number_of_pairs+1)))
    else:
      self.pairs_per_iter = (3,int(na.log(number_of_pairs+1)/na.log(3) ))
    
    length = end-start
    cur_dir = os.path.abspath(os.curdir)
    
    seq = self.ref.get_sequence(contig, start, end)
    print "#Stage1 Region Length: %d"%(length)
    if temp_tag is None:
      out_dir=CONFIG.dir['default_temp']
      temp_tag = tempfile.mktemp(prefix='', dir=out_dir)
    elif os.path.isdir(temp_tag):
      out_dir=temp_tag
      temp_tag = tempfile.mktemp(prefix='', dir=out_dir)
    
    primer3_in, primer3_out = None, None
    try:
      pre_primer3_time = time.time()
      if not (os.path.isfile('%s.primer3'%temp_tag) and os.path.isfile('%s.primer3.out'%temp_tag)):
        os.chdir(self.primer3_path)  
        primer3_in = open('%s.primer3'%temp_tag, 'wb')
        print >>primer3_in, "PRIMER_NUM_RETURN=%d"%((length/1000)*primers_per_kbp)

        print >>primer3_in, "SEQUENCE_ID=%s:%d-%d"%(contig, start, end)
        print >>primer3_in, "SEQUENCE_TEMPLATE=%s"%seq
        
        print >>primer3_in, "="
        primer3_in.close()
        print os.path.abspath(os.path.curdir)
        primer3_cmd = "%s -p3_settings_file=%s %s"%(os.path.join(self.primer3_path, 'src', 'primer3_core'), self.primer3_param, primer3_in.name)
        print primer3_cmd
        primer3_out = open("%s.primer3.out"%temp_tag, 'wb')
        p_primer3 = Popen(shlex.split(primer3_cmd), stdout=primer3_out, stderr=PIPE)
        
        (primer3_out_str, primer3_err) = p_primer3.communicate()
        
        primer3_out.close()
        os.chdir(cur_dir)
        
      with open("%s.primer3.out"%temp_tag, 'rb') as primer3_out_f:
        self.get_primers(primer3_out_f, max_penalty=max_primer_penalty)
      print "#Stage1 Primer3Time: %.4f"%(time.time()-pre_primer3_time)

      # Compress
      map_repeat_seqs_f = dict([(info_dict['SEQUENCE'],info_dict[''][0]) for info_dict in self.forward_dict.itervalues()])
      map_repeat_seqs_r = dict([(info_dict['SEQUENCE'],info_dict[''][0]) for info_dict in self.reverse_dict.itervalues()])

      pre_align_time = time.time()
      aligner_out_fpath = '%s.align.out'%temp_tag
      if not (os.path.isfile('%s.align.fa'%temp_tag) and os.path.isfile(aligner_out_fpath)):
        
        aligner_in = open('%s.align.fa'%temp_tag, 'wb')

        primers_fa = [">LEFT:%d\n%s"%(pos,seq) for seq,pos in map_repeat_seqs_f.iteritems()]
        primers_fa += [">RIGHT:%d\n%s"%(pos,seq) for seq,pos in map_repeat_seqs_r.iteritems()]

        #primers_fa = [">LEFT:%d\n%s"%(info_dict[''][0],info_dict['SEQUENCE']) for info_dict in self.forward_dict.itervalues()]
        #primers_fa += [">RIGHT:%d\n%s"%(info_dict[''][0],info_dict['SEQUENCE']) for info_dict in self.reverse_dict.itervalues()]
        
        print >>aligner_in, "\n".join(primers_fa)
        aligner_in.close()

        if self.aligner.endswith('blat'):
          align_cmd = '%s -t=dna -q=dna -stepSize=2 -tileSize=11 -repMatch=1048576 -minScore=15 -noHead %s %s %s'%(self.aligner, CONFIG.param['reference_fpath'], aligner_in.name, aligner_out_fpath)
        print len(primers_fa), align_cmd
        p_align = Popen(shlex.split(align_cmd), stdout=PIPE, stderr=PIPE)
        
        (aligner_out, aligner_err) = p_align.communicate()
      
      with open(aligner_out_fpath, 'rb') as aligner_out_f:
        self.set_alignments_dict(aligner_out_f, max_align=max_alignment_count, min_align_len=min_alignment_len)
      
      print "#Stage2 AlignTime: %.4f"%(time.time()-pre_align_time)
      
      for_filtered_idx = [k for k, info_dict in self.forward_dict.iteritems() if info_dict[''][0] in self.alignments_dict.viewkeys()]
      rev_filtered_idx = [k for k, info_dict in self.reverse_dict.iteritems() if info_dict[''][0] in self.alignments_dict.viewkeys()]
      self.pairs_idx = set(for_filtered_idx).intersection(set(rev_filtered_idx))
      self.pairs_pos = set([(self.forward_dict[idx][''][0],self.reverse_dict[idx][''][0]) for idx in self.pairs_idx])
      print "#Number of Primer Pairs after BLAT: %d"%len(self.pairs_pos)
      self.check_cross_amplification(max_dist=max_cross_amp_dist)      
 
      pre_log_select_time = time.time()
      self.selected_pairs = self.primer_pairs_even_spaced(length, number_of_pairs)
      print "#Stage3 LogSelectTime: %.4f"%(time.time()-pre_log_select_time)
      
      self.print_solutions(out_fpath=out_fpath)
  
      
    finally:
      if not primer3_in is None:
        primer3_in.close()
        print primer3_in.name
      if not primer3_out is None:
        primer3_out.close()
        print primer3_out.name
      
      if(delete_flag):
        os.remove(primer3_in.name)
        os.remove(primer3_out.name)
  def print_solutions(self, out_fpath):
    if out_fpath is None:
      out = sys.stdout
    else:
      out = open(out_fpath, 'wb')
    
    if self.pairs_per_iter[0] == 2:
      annot = self.recurse_2annotate(self.selected_pairs.size)
    elif self.pairs_per_iter[0] == 3:
      annot = self.recurse_3annotate(self.selected_pairs.size)
     
    for idx, (iteration, label) in zip(self.selected_pairs, annot):
      if na.isnan(idx): continue 
      for_info_dict = self.forward_dict[idx]
      print >>out, ">%d_%d%s_T\n%s"%(for_info_dict[''][0], iteration, label, for_info_dict['SEQUENCE'])
      rev_info_dict = self.reverse_dict[idx]
      print >>out, ">%d_%d%s_F\n%s"%(rev_info_dict[''][0], iteration, label, rev_info_dict['SEQUENCE'])
    
    if not out_fpath is None:
      out.close()
    else:
      out.flush()  

def ambre_run(region_str, number_of_pairs, temp_tag=None, out_fpath=None):
  w = PrimerBinarySearchWorkflow(primer3_path=CONFIG.dir['primer3'], 
          primer3_param=CONFIG.param['primer3_short_param'],
          aligner=CONFIG.bin['aligner'])

  contig, pos_str = region_str.split(':')
  start, end = map(int,pos_str.split('-'))
  
  w.run(contig, start, end, number_of_pairs,
        delete_flag=(CONFIG.param['cleanup_flag']=="True"),
        temp_tag=temp_tag,
        out_fpath=out_fpath,
        max_primer_penalty=float(CONFIG.param['design_max_primer3_penalty']),
        primers_per_kbp=int(CONFIG.param['design_primer3_primers_per_kbp']),
        max_cross_amp_dist=int(CONFIG.param['design_max_cross_amp_dist']),
        max_alignment_count=int(CONFIG.param['design_max_alignments']),
        min_alignment_len=int(CONFIG.param['design_3end_len_alignment']))
  
def main(): 
  parser = argparse.ArgumentParser(prog='ambre-design.py', description='Select compatible primers covering reference region.')
  
  parser.add_argument('reference', type=str, nargs=1, help='Fasta with multiple sequence entries.')
  parser.add_argument('region', type=str, nargs=1, help='contig:start-end format string')
  parser.add_argument('primer_count', type=int, nargs=1, help='Number of primer pairs to select (must be one less than a power of 2 or 3)')
  parser.add_argument('temptag', type=str, nargs='?', default=None, help='Prefix for the run id / Directory to store Temps of a run id')  
  parser.add_argument('-o, --primers-out', type=str, nargs=1, default=[None], dest="primer_fpath", help='Output primers to tab-delimited file. Default is to output to STDOUT')
  parser.add_argument('--config', type=str, nargs=1, default=None, dest="config_fpath", help='Update parameters in default config file with new config file.')
   
  args = parser.parse_args()
  
  CONFIG.param['reference_fpath'] = args.reference[0]
  if not args.config_fpath is None:
    CONFIG.update(args.config_fpath[0])
  
  primer_count = args.primer_count[0]
  
  ambre_run(args.region[0], primer_count, args.temptag, args.primer_fpath[0])

if __name__=='__main__':
  main()
