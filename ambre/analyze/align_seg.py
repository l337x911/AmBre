'''
#  ambre.analyze.align_seg.py
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

import ambre.utils.SAM as SAM
import numpy as na
import sys
import bisect
from collections import defaultdict
from utils import cigar_listparse, get_frag_interval, CIGARAlignmentScoring, get_mdi_cigar
import itertools
debug_log = sys.stderr


class Alignments(object):
  def __init__(self):
    self.frags = defaultdict(lambda : defaultdict(list))
    self.frag_lengths = {}
    self.align_scorer = None
  def get_frag_length(self, qname):
    return self.frag_lengths[qname]
  def get_alignments(self, qname, rname):
    return self.frags[qname][rname]
  def get_aligned_frag(self, fp, fe, fl):
    frag_matched_bp = na.zeros((fl,), dtype=na.int16)
    if fp>fe:
      frag_matched_bp[fe:fp+1] = True
    else:
      frag_matched_bp[fp:fe+1] = True
    return frag_matched_bp
  def get_ascore(self, astr_parts): 
    return self.align_scorer.get_score(astr_parts)
  def get_segment_rp(self, md_array, i_array, fp_c, ref_start, ref_end, rc):
    
    rp_c = []
    for i in xrange(1,len(fp_c)):
      rp_c.append(na.sum(md_array[fp_c[i-1]:fp_c[i]]!='-')+na.sum(i_array[fp_c[i-1]:fp_c[i]]))
    cumsum_rp_c = na.cumsum(rp_c)
    if rc:
      nrp = na.concatenate(([ref_start],(ref_start-cumsum_rp_c)))
      #assert nrp[-1] == ref_end
    else:
      nrp = na.concatenate(([ref_start],(ref_start+cumsum_rp_c)))
      #assert nrp[-1] == ref_end
    return nrp
    
  def str_concatenate(self, a_strs_list):
    pass
class SAMAlignments(Alignments):
  def __init__(self, sam_fpath):
    Alignments.__init__(self)
    parser = SAM.SAMParser()
    
    self.sam_fpath = sam_fpath
    rev_flag_idx = SAM.SAM_FLAGS_H['rev_strand_of_query']

    for frag in parser.parse(sam_fpath):
      cig_parts = cigar_listparse(frag)
      c = SAM.get_flag(rev_flag_idx, int(frag.flags))
      start, end, frag_length = get_frag_interval(cig_parts, c)
      #if abs(int(frag.tlen))<150:
      #  continue
      self.frags[frag.qname][frag.rname].append((int(frag.pos), int(frag.pos) + int(frag.tlen), start, end, cig_parts, c))
      
      self.frag_lengths[frag.qname] = frag_length
      
    self.align_scorer = CIGARAlignmentScoring()
    
  def get_mdi(self, a_str):
    return get_mdi_cigar(a_str)
  
  def get_segment_a_str(self, fp, a_str, ref_start, ref_end, rc):
    if len(a_str)>0:
      if a_str[0][1]=='H' or a_str[0][1]=='S':
        a_str = a_str[1:]
      if len(a_str)>0:
        if a_str[-1][1]=='H' or a_str[-1][1]=='S':
          a_str = a_str[:-1]
          
    if rc:
      a_str = a_str[::-1]
    md_array, i_array, idx_array = get_mdi_cigar(a_str)
    
    assert (fp[-1]-fp[0])==md_array.size

    fp_c = fp-fp[0]
    # In the context of the fragment
    nrp = self.get_segment_rp(md_array, i_array, fp_c, ref_start, ref_end, rc)
    
    # Requires cigar informative splitting of the a_str at segment parts
    # Create a new a_str with proper breaks and idx array adjustments
    # with 
    segment_a_str = []
    for i in xrange(1,len(fp_c)-1):
      a_str_idx = idx_array[fp_c[i]] 
      if a_str_idx==idx_array[fp_c[i]-1]:
        count_to_split = na.sum(idx_array[:fp_c[i]]==a_str_idx)
        
        c, t = a_str.pop(a_str_idx)
        a_str.insert(a_str_idx, (c-count_to_split, t))
        a_str.insert(a_str_idx, (count_to_split, t))
        idx_array[fp_c[i]:] += 1
        
      a_str_idx = idx_array[fp_c[i]] 
      assert a_str_idx!=idx_array[fp_c[i]-1]
      segment_a_str.append(a_str[idx_array[fp_c[i-1]]:idx_array[fp_c[i]]]) 
    
    segment_a_str.append(a_str[idx_array[fp_c[-2]]:])
    return segment_a_str, nrp
  def get_a_str(self, cig_parts):
    return ''.join(["%d%s"%a for a in cig_parts])
  def str_concatenate(self, a_strs_list):
    
    for i in xrange(1, len(a_strs_list)):
      c,t = a_strs_list[i][0]
      
      if a_strs_list[i-1][-1][1]==t:
        a_strs_list[i][0] = (a_strs_list[i-1][-1][0]+c, t)
        a_strs_list[i-1].pop(-1)
        
        
    return ''.join([self.get_a_str(a_str) for a_str in a_strs_list])
      
class MaxAlignmentScoreFragFiltering(object):
  def __init__(self, fpath):
    self.fpath = fpath
    self.a = SAMAlignments(fpath)

  def check_overlaps(self, frag_alignments, top_a_scores=None, remove_encompassed=False):
    
    rnames, fps, fes, rps, res, a_strs, rcs = zip(*frag_alignments)
    a_scores = [self.a.get_ascore(a_str) for a_str in a_strs]
    
    frag_alignments = zip(a_scores, rnames, fps, fes, a_strs, rcs)
    frag_alignments.sort()
    if not top_a_scores is None:
      frag_alignments = frag_alignments[:top_a_scores]
    n = len(frag_alignments)
    a_scores, rnames, fps, fes, a_strs, rcs = zip(*frag_alignments)
    
    breakpoints = list(set(zip(rnames, fps)+zip(rnames, fes)))
    breakpoints.sort()
    
    m = len(breakpoints)
    nodes = na.zeros((m, n), dtype=na.int)
    
    for idx, rname, fp, fe in zip(range(n), rnames, fps, fes):
      i,j = bisect.bisect_left(breakpoints, (rname, fp)), bisect.bisect_right(breakpoints, (rname, fe))
      nodes[i:j,idx] = 1
    if remove_encompassed:
      mark_for_deletion = []
      for j in xrange(n):
        a = nodes[:,j]  
        if na.sum(na.dot(a.reshape((1,m)), nodes)==na.sum(a))>1:
          mark_for_deletion.append(j)
      mask = na.ones(n, dtype=na.bool)
      mask[mark_for_deletion] = False
      nodes = nodes[:, mask]
    n = nodes.shape[1]
    # count overlaps
    # count encompasses
    c_overlaps, c_encompasses, encompassment = 0,0,0

    for j in xrange(n):
      a = nodes[:,j]  
      similar_nodes = na.dot(a.reshape((1,m)), nodes)
      c = na.sum(similar_nodes==na.sum(a))
      c_encompasses += c-1
      c_overlaps += na.sum(similar_nodes>0)-1
      if c>1:
        encompassment += 1

    print >>debug_log, n, c_overlaps, c_encompasses, encompassment
    
  def max_scoring_path(self, frag_alignments, breakpoint_weight=-50):
    n = len(frag_alignments)
    frag_alignments.sort()
    rnames, fps, fes, rps, res, a_strs, rcs = zip(*frag_alignments)
    
    breakpoints = list(set(fps).union(set(map(lambda x:x+1, fes))))
    breakpoints.sort()
    m = len(breakpoints)
    
    nodes = na.zeros((m, n), dtype=na.int)
      
    for idx, fp, fe in zip(range(n), fps, fes):
      i,j = bisect.bisect_left(breakpoints, fp), bisect.bisect_right(breakpoints, fe+1)
      nodes[i:j,idx] = 1
      
    # Removes alignments that are encompassed by another alignment.
    mark_for_deletion = []
    for j in xrange(n):
      a = nodes[:,j]  
      if na.sum(na.dot(a.reshape((1,m)), nodes)==na.sum(a))>1:
        mark_for_deletion.append(j)
        
    mask = na.ones(n, dtype=na.bool)
    mask[mark_for_deletion] = False
    
    if not mask.any():
      return [], -1
    rnames, fps, fes, rps, res, a_strs, rcs = zip(*itertools.compress(frag_alignments, mask))
    
    n = len(rnames)
    
    if n==1:
      a_score = self.a.get_ascore(a_strs[0])
      if rcs[0]:
        return [(rnames[0], res[0], rps[0], fps[0], fes[0], rcs[0], self.a.get_a_str(a_strs[0]), a_score)], a_score
      else:
        return [(rnames[0], rps[0], res[0], fps[0], fes[0], rcs[0], self.a.get_a_str(a_strs[0]), a_score)], a_score
    
    breakpoints = list(set(fps).union(set(map(lambda x:x+1, fes))))
    breakpoints.sort()
    m = len(breakpoints)
    
    nodes = na.zeros((m, n), dtype=na.int)
    break_rev_edges = defaultdict(set)
    frag_rev_edges = dict()
    frag_rev_edge_a_scores = dict()
    prev_alignment_termination = None

    max_node_dict = dict()
    min_node_dict = dict()
    
    all_segment_strs, all_n_rps, all_segment_scores = [],[],[]
    for idx, a_str, fp, fe, ref_start, ref_end, rc in zip(range(n), a_strs, fps, fes, rps, res, rcs):
      i,j = bisect.bisect_left(breakpoints, fp), bisect.bisect_right(breakpoints, fe+1) 
      nodes[i:j,idx] = 1
      max_node_dict[idx] = j-1
      min_node_dict[idx] = i
      
      # Note that segment_strs length is j-i-1 and length of n_rps is j-1
      bps = na.array(breakpoints[i:j])
      #bps[-1]-=1
      
      if rc:
        segment_strs, n_rps = self.a.get_segment_a_str(bps, a_str, ref_end, ref_start, rc)
      else:
        segment_strs, n_rps = self.a.get_segment_a_str(bps, a_str, ref_start, ref_end, rc)
      
      segment_scores = na.array(map(self.a.get_ascore, segment_strs))
      all_segment_strs.append(segment_strs)
      all_n_rps.append(n_rps)
      all_segment_scores.append(segment_scores)
      # segment edges
      for k in range(i+1,j):
        frag_rev_edges[(k,idx)] = (k-1,idx)
        frag_rev_edge_a_scores[(k,idx)] = segment_scores[k-(i+1)]
      
      # breakpoint edges
      if prev_alignment_termination is None:
        prev_alignment_termination = (j-1,idx)
      else:
        # get more than two overlaps
        a = nodes[:,idx]
        assert a.any()
        first_bp_idx = na.nonzero(a)[0][0]
        
        similar_nodes = na.dot(a.reshape((1,m)), nodes[:,prev_alignment_termination[1]:idx])
        idx_overlaps = na.nonzero(similar_nodes>1)[1]+prev_alignment_termination[1]
        
        # Add termination connection edge
        if idx_overlaps.size == 0:
          prev_alignment_termination = (max_node_dict[idx-1], idx-1)
          break_rev_edges[(i,idx)].add(prev_alignment_termination)
        else:
          # Handle all the overlapping alignments
          for o_idx in idx_overlaps:
            break_rev_edges[(max_node_dict[o_idx],idx)].add((first_bp_idx, o_idx))
          
          prev_term_idx = idx_overlaps[0]-1
          if prev_term_idx>=0:
            prev_alignment_termination = (max_node_dict[prev_term_idx], prev_term_idx)
            break_rev_edges[(i,idx)].add(prev_alignment_termination)

    dp_score = na.zeros_like(nodes)
    dp_pointer = na.ones((m,n,2), dtype=na.int)*-1
    for k,idx in zip(*na.nonzero(nodes)):
      mx = []
      if (k,idx) in frag_rev_edges.viewkeys():
        i,j = frag_rev_edges[(k,idx)]
        mx.append((dp_score[i,j]+frag_rev_edge_a_scores[(k,idx)], (i,j)))
      for (i,j) in break_rev_edges[(k,idx)]:
        overlap_break_score = 0
        if (k,j) in frag_rev_edges.viewkeys() and (k,idx) in frag_rev_edges.viewkeys():
          overlap_break_score = (frag_rev_edge_a_scores[(k,idx)]+frag_rev_edge_a_scores[(k,j)])/2
        mx.append((dp_score[i,j]+overlap_break_score+breakpoint_weight, (i,j)))
      if len(mx)>0:
        mx.sort()
        dp_score[k,idx] = mx[-1][0]
        dp_pointer[k,idx,:] = mx[-1][1]
    mx = [(dp_score[k,idx], (k,idx)) for idx,k in max_node_dict.iteritems()]
    mx.sort()
    
    path = [[mx[-1][1][1], mx[-1][1][0], mx[-1][1][0]]]
    while True:
      node = path[-1]
      i,j = dp_pointer[node[1], node[0], :]
      if j == path[-1][0]:
        path[-1][1] = i
        continue
      
      if i>=0 and j>=0:
        path.append([j,i,i])
      else:
        break
    
    return [(rnames[n_idx], 
             all_n_rps[n_idx][bp_i_idx-min_node_dict[n_idx]], 
             all_n_rps[n_idx][bp_j_idx-min_node_dict[n_idx]], 
             breakpoints[bp_i_idx], 
             breakpoints[bp_j_idx]-1,
             rcs[n_idx],
             self.a.str_concatenate(all_segment_strs[n_idx][(bp_i_idx-min_node_dict[n_idx]):(bp_j_idx-min_node_dict[n_idx])]),
             na.sum(all_segment_scores[n_idx][(bp_i_idx-min_node_dict[n_idx]):(bp_j_idx-min_node_dict[n_idx])])) for (n_idx,bp_i_idx, bp_j_idx) in path[::-1]], mx[-1][0]
 
  def frag_isoform_filtering(self, frag_fpath=None):
    if frag_fpath is None:
      frag_file = sys.stdout
    else:
      frag_file = open(frag_fpath, 'wb')
      
    print >>frag_file, "#qname\trname\tref_start\tref_end\tfrag_start\tfrag_end\tfrag_len\talign_str\tascore"
    
    for frag_name in self.a.frags.iterkeys():
      frag_alignments = []
      for rname, alignments in self.a.frags[frag_name].iteritems():
        n = len(alignments)
        ref_start, ref_end, fp, fe, a_strs, rc = zip(*alignments)
        
        fe, fp, rc = na.array(fe, dtype=na.int), na.array(fp, dtype=na.int), na.array(rc, dtype=na.bool)
        t = na.array(fe)
        fe[rc] = fp[rc]
        fp[rc] = t[rc]
        # both fragment and reference are ascending intervals.
      
        frag_alignments.extend(zip([rname]*n, fp, fe, ref_start, ref_end, a_strs, rc))
      frag_alignments.sort()
      
      #self.check_overlaps(frag_alignments, top_a_scores=None, remove_encompassed=True)
      f_alignments, tot_a_score = self.max_scoring_path(frag_alignments)
      
      for rname, ref_start, ref_end, fp, fe, rc, a_str, a_score in f_alignments:
        # Note that talign is in the context of the fragment and not the reference.
        
        # arises based on trivial overlaps. 
        if fp==fe:
          continue
        print >>frag_file, "\t".join((frag_name, rname, "%d"%ref_start, "%d"%ref_end, "%d"%fp, "%d"%fe, "%d"%self.a.get_frag_length(frag_name), a_str, "%0.2f"%a_score))
    if not frag_fpath is None:
      frag_file.close()
  
