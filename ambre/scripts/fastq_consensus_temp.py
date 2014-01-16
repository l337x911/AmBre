import numpy as na

with open('consensus.fastq', 'rb') as f:
  lines = [line.strip() for line in f]

templates = {}

for i in xrange(0,len(lines),4):
  templates[lines[i]] = (lines[i+1],lines[i+3])

def longest_string_of_nonchar(nonchar, seq_arr):
  prefix_count_arr = na.zeros(seq_arr.size, dtype=na.int)
  nonmatch = na.logical_not(seq_arr==nonchar)
  for i in xrange(seq_arr.size):
    if nonmatch[i]:
      prefix_count_arr[i] = int(nonmatch[i])+prefix_count_arr[i-1]
  argmax = na.argmax(prefix_count_arr)
  return prefix_count_arr[argmax], argmax-prefix_count_arr[argmax]+1, argmax+1

for tname, (seq,qual) in templates.iteritems():
  qual_arr = na.fromstring(qual, dtype='S1')
 
  max_lc,i,j = longest_string_of_nonchar('!',qual_arr)
  if max_lc==0: continue
  print ">%s;%d-%d(%d)"%(tname, i,j,max_lc)
  print "%s"%seq[i:j]
  print "+"
  print "%s"%qual[i:j]





