'''
#  ambre.test.reference.py
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
import os
import string

BASE_COMPLEMENT = string.maketrans("ACGTNacgtn", "TGCANtgcan")

class Reference(object):
  '''
  classdocs
  '''

  def __init__(self, fpath):
    '''
    Constructor.
    Relies on .fai index (consistent with samtools faidx format)
    '''
    self.ref_fpath = fpath
    self.fasta = None
    self.fasta_byte_offset_dict = {}
    self.index(fpath + '.fai')
    
  def index(self, idx_fpath):
    self.fasta_byte_offset_dict = {}
    if os.path.isfile(idx_fpath):
      self.get_index_from_file(idx_fpath)
    else:
      self.set_index(idx_fpath)
  
  def get_index_from_file(self, idx_fpath):
    
    with open(idx_fpath, 'rb') as idx_file:
      for l in idx_file:
        tokens = l.strip().split('\t')
        self.fasta_byte_offset_dict[tokens[0]] = (int(tokens[2]), int(tokens[1]), int(tokens[3]), int(tokens[4])) 
    
  def set_index(self, idx_fpath):
    keys, starts, segments, bytes_in_lines, bases_in_lines = [], [], [], [], []
    
    self.fasta = open(self.ref_fpath, 'rb')

    bases_in_line = 0
    bytes_in_line = 0
    
    position = 0
    prev = 0
    for l in self.fasta:
      prev = position
      position += len(l)
      
      if l[0] == '>': 
        keys.append(l.strip()[1:])
        starts.append(position)
        # top off previous entry
        segments.append(prev) 
        bytes_in_lines.append(bytes_in_line)
        bases_in_lines.append(bases_in_line)
        
        bases_in_line = 0
        bytes_in_line = 0
      else:
        bytes_in_line = max(bytes_in_line, len(l))
        bases_in_line = max(bases_in_line, len(l.strip()))
    
    bytes_in_lines.append(bytes_in_line)
    bases_in_lines.append(bases_in_line)
    segments.append(position)

    self.fasta.close()
    
    out_idx_file = open(idx_fpath, 'wb')

    for k, start_offset, end_offset, bytes_in_line, bases_in_line in zip(keys, starts, segments[1:], bytes_in_lines[1:], bases_in_lines[1:]):
      contig_length = end_offset - start_offset
      if contig_length % bytes_in_line ==0:
        contig_length -= (contig_length / bytes_in_line)
      else:
        contig_length -= 1 + ( contig_length / bytes_in_line)
      print >> out_idx_file, "\t".join([k, "%d" % contig_length, "%d" % start_offset, "%d" % bases_in_line, "%d" % bytes_in_line])
      self.fasta_byte_offset_dict[k] = (start_offset, contig_length, bases_in_line, bytes_in_line)

    out_idx_file.close()
  def get_contig_names(self):
    return self.fasta_byte_offset_dict.keys()
  def get_contig(self, key):
    try:
      k_s, k_e, bases_in_line, bytes_in_line = self.fasta_byte_offset_dict[key]
    except KeyError:
      print >>sys.stderr, "Contig not in reference %s found %s."%(self.ref_fpath, key)
      sys.exit(1)
    return self.get_sequence(key, 0, k_e)
  def get_contig_length(self, key):
    k_s, k_e, bases_in_line, bytes_in_line = self.fasta_byte_offset_dict[key]
    return k_e
  def get_sequence(self, key, start, end, rev_complement=False):
    k_s, k_e, bases_in_line, bytes_in_line = self.fasta_byte_offset_dict[key]
    
    new_lines_s = start / bases_in_line
    new_lines_e = end / bases_in_line
    
    seq_length = end - start + new_lines_e - new_lines_s
  
    assert end <= k_e
    with open(self.ref_fpath, 'rb') as self.fasta: 
      self.fasta.seek(k_s + start + new_lines_s)
      if rev_complement:
        seq = self.fasta.read(seq_length).replace('\n', '')[::-1].translate(BASE_COMPLEMENT)
      seq = self.fasta.read(seq_length).replace('\n', '')
    return seq
  def close(self):
    if not self.fasta is None:
      self.fasta.close()

def test():
  r = Reference('/home/anand/data/hg/hg19/full.fa')
  y = r.get_sequence('chr1', 9980, 10080)
  x = r.get_sequence('chr1', 23000000, 23000100)

  print len(y), y
  print len(x), x
  
if __name__ == '__main__':

  import sys
  fasta_fpath, filter_fpath = sys.argv[1], sys.argv[2]
  if filter_fpath == '-':
    f = sys.stdin
  else:
    f = open(filter_fpath, 'rb')
  
  contigs = Reference(fasta_fpath)
  for line in f:
    try:
      seq = contigs.get_contig(line.strip())
    except KeyError:
      print >> sys.stderr, "Not Present in FASTA, %s" % (line.strip())
  f.close()
  contigs.close()
  
