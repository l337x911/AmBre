import sys
from collections import defaultdict
from itertools import izip

# assume mrna

fa_fpath = sys.argv[1]
region_count = int(sys.argv[2])

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

with open(fa_fpath, 'rb') as f:
  fas = [line.strip() for line in f]

t = zip(fas[::2], fas[1::2])

slns = defaultdict(set)

for name, seq in izip(fas[::2], fas[1::2]):
  sln_idx, primer_idx = get_id(name)
  slns[sln_idx].add(primer_idx)

print "# Slns with a primer in each region"
for sln_idx in sorted(slns.iterkeys()):
  region_with_primers = slns[sln_idx]
  if len(region_with_primers) != region_count: continue
  print "%02d: "%sln_idx, sorted(list(region_with_primers))





  
