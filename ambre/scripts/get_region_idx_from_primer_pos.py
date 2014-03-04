from ambre.design.workflow import ambre_run, CONFIG
from collections import defaultdict
import sys

CONFIG.param['reference_fpath'] = sys.argv[2]
CONFIG.update(sys.argv[1])
regions_fpath, temp_tag = sys.argv[3], sys.argv[4]
multiplex_fpath = "%s.multiplx.out"%temp_tag

w = ambre_run(regions_fpath, temp_tag)

primer_set = set()
with open(multiplex_fpath, 'rb') as f:
  for line in f:
    if line.startswith('Name'):
      continue
    tokens = line.strip().split('\t')
    primer_set.add(int(tokens[0]))
    primer_set.add(int(tokens[1]))

max_primer_lengths = defaultdict(int)
   
primers = []
for v in w.forward_dict.itervalues():
  pos,seq = v[''][0], v['SEQUENCE']
  if not pos in primer_set: continue
  idx, regions = w.get_original_region(pos)
  if v[''][1]>max_primer_lengths[pos]:
    max_primer_lengths[pos] = v[''][1]
  primers.append((pos, ">%d-%d-forward"%(pos,idx), seq))

for v in w.reverse_dict.itervalues():
  pos,seq = v[''][0], v['SEQUENCE']
  if not pos in primer_set: continue
  idx, regions = w.get_original_region(pos)
  if v[''][1]>max_primer_lengths[pos]:
    max_primer_lengths[pos] = v[''][1]
  primers.append((pos, ">%d-%d-reverse"%(pos,idx), seq))

idx =0
for p,n,s in sorted(primers):
  if max_primer_lengths[p]>len(s): continue
  idx += 1
  print n
  print s

print len(primer_set), idx
