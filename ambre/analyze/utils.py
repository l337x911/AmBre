'''
'''
#  ambre.analyze.utils.py
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
'''
import re
import numpy as na
#PRECISE_BP=50
#PRECISE_EMPTY_BP= 20
#PRECISE_ALIGN_BP= 100
from ambre.config import CONFIG

cigar_re = re.compile('(\d+)([MIDNSHPX])')
# -A is gap in query, A- is gap in subject
btop_re = re.compile('(\d+)|((?:[ACGTN][ACGTN])+)|((?:-[ACGTN])+)|((?:[ACGTN]-)+)')

def cigar_listparse(record):
  if record.cigar == '*':
    return []
  elif record.cigar == '=':
    return [(len(int(record.tlen)),'M') ]
  else:
    return [(int(m.group(1)), m.group(2)) for m in re.finditer(cigar_re, record.cigar)]

def cigar_str_listparse(record):
  if record == '*':
    return []
  return [(int(m.group(1)), m.group(2)) for m in re.finditer(cigar_re, record)]

def get_mdi_cigar(cig_list):
  # Assumes insertion is inserted after a fragment position
  # The first position must be a Match or a Deletion
  # in the context of the fragment
  # i_str length is n-1
  md_str = ''
  i_str = []
  idx_str = []
  for idx, (i,t) in enumerate(cig_list):
    if t=='H' or t=='S': continue
    elif t=='M':
      md_str += 'M'*i
      i_str.extend([0]*i)
      idx_str.extend([idx]*i)
    elif t=='I':
      md_str += '-'*i
      i_str.extend([0]*i)
      idx_str.extend([idx]*i)
    elif t=='D':
      i_str[-1]= i
  return na.fromstring(md_str, dtype='S1'), na.array(i_str[:-1]), na.array(idx_str)
  
def get_frag_interval(cig_list, rc_flag):
  frag_length = sum([a for a,b in cig_list if b!='D'])
  if rc_flag:
    start = frag_length if cig_list[0][1] != 'H' and cig_list[0][1] != 'S' else frag_length-cig_list[0][0]
    end = 1 if cig_list[-1][1] != 'H' and cig_list[-1][1] != 'S' else cig_list[-1][0]+1
  else:
    start = 1 if cig_list[0][1] != 'H' and cig_list[0][1] != 'S' else cig_list[0][0]+1
    end = frag_length if cig_list[-1][1] != 'H' and cig_list[-1][1] !='S' else frag_length-cig_list[-1][0]
  
  return start, end, frag_length

def get_subalign_interval(align_range, cig_list, rc_flag):
  frag_length = sum([a for a,b in cig_list if b!='D'])
  
  sidx, eidx = align_range
  
  start_length = sum([a for a,b in cig_list[:sidx] if b!='D'])
  end_length = sum([a for a,b in cig_list[:eidx] if b!='D'])
  
  if rc_flag:
    return frag_length-start_length, frag_length-end_length
  else:
    return start_length, end_length

def get_ref_length(cig_list):
  return sum([a for a,b in cig_list if b!='I' and b!='H' and b!='S'])



MISMATCH_FRAC=float(CONFIG.param['analyze_mismatch_fraction'])
MATCH_FRAC=1-MISMATCH_FRAC
MATCH_SCORE=float(CONFIG.param['analyze_match_score'])
MISMATCH_PENALTY=float(CONFIG.param['analyze_mismatch_penalty'])
GAPOPEN_PENALTY=float(CONFIG.param['analyze_gapopen_penalty'])
GAPEXT_PENALTY=float(CONFIG.param['analyze_gapext_penalty'])

class AlignmentScoring(object):
  def __init__(self, m=MATCH_SCORE, mis=MISMATCH_PENALTY,go=GAPOPEN_PENALTY, ge=GAPEXT_PENALTY):
    self.match = m
    self.mismatch = mis
    self.gap_opening = go
    self.gap_extension = ge
  def no_score(self, length): return 0
  def match_score(self, length): return self.match*length  
  def mismatch_score(self, length): return self.mismatch*length
  def gap_score(self, length): return self.gap_opening+self.gap_extension*length

class CIGARAlignmentScoring(AlignmentScoring):
  def __init__(self, m=MATCH_SCORE, mis=MISMATCH_PENALTY,go=GAPOPEN_PENALTY, ge=GAPEXT_PENALTY):
    AlignmentScoring.__init__(self, m, mis, go, ge)
    # Unfortunately, BLASR only reports an alignment match (not whether the basepair is a mismatch or a match)
    # hence the empirical 2% mismatch error rate.
    self.score_dict = {'M': self.match_score, 'D': self.gap_score, 'I': self.gap_score, 'H': self.no_score, 'S':self.no_score}
  def match_score(self, length): return (MATCH_FRAC*self.match+MISMATCH_FRAC*self.mismatch)*length  
  def get_score(self, cigar_list):
    return sum([self.score_dict[b](a) for a,b in cigar_list])
  
