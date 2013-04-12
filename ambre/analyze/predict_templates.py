import ambre.utils.reference as reference
from ambre.config import CONFIG
import sys
from operator import itemgetter
from itertools import count
from collections import defaultdict
from frag_clustering import BreakPoint

# breakpoint positions left and right breakpoints
# negative left is forward (positive is forward)
# positive right is forward (negative is reverse)

class PredictedTemplates(object):
  def __init__(self, fpath, ref_fpath, contig):
    self.fpath = fpath
    # split fasta file by 80 characters
    self.line_width = 80
    
    self.min_sv_dist = int(CONFIG.param['analyze_min_sv_dist'])
    self.max_d = int(CONFIG.param['analyze_max_d'])
    
    self.k = int(CONFIG.param['analyze_default_segment'])
    # Contains subset of bp that will generate into
    # templates. 
    self.bps_set = {}
    
    ref = reference.Reference(ref_fpath)
    self.contig = ref.get_contig(contig)
    
    # Recalls ordering of bps.
    self.bps_ordering = {}
    self._parse_bps()
  def _filtering(self, tokens, bp):
    return tokens[2]>self.max_d or abs(abs(tokens[0])-abs(tokens[1]))<self.min_sv_dist
  def _parse_bps(self):
    
    with open(self.fpath, 'rb') as f:

      for idx in count():
        bp = BreakPoint()
        try:
          bp.from_file(f)
        except StopIteration:
          break
        
        tokens = list(bp.get_pair_mode())
        a,b = tokens[:2]
    
        # assume opposite was there
        try:
          m = self.bps_set[(b,a)]
          if m<tokens[2]:
            tokens[2]=m
          if abs(b)<abs(a):
            tokens[0], tokens[1] = b,a
          del self.bps_set[(b,a)]
        except KeyError:
          pass

        # noise filtering
        if self._filtering(tokens, bp): continue
        
        self.bps_set[tokens[0],tokens[1]] = tokens[2]
        self.bps_ordering[(tokens[0],tokens[1], tokens[2])] = (idx, bp)
        
    x,d = zip(*self.bps_set.iteritems())
    a,b = zip(*x)

    self.bps_set = set(zip(a,b,d))
  def _forward_extract(self, i, offset=None):
    if offset is None: offset = self.k
    return self.contig[i-1:i+offset]
  
  def _reverse_extract(self, i, offset=None):
    if offset is None: offset = self.k
    return self.contig[i-1:i+offset][::-1].translate(reference.BASE_COMPLEMENT)
  def _forward_seg_extract(self, i, j):
    # extract i,j are 1-based positions 
    return self.contig[i-1:j]
  def _reverse_seg_extract(self, i, j):
    return self.contig[i-1:j][::-1].translate(reference.BASE_COMPLEMENT)
  
  def _format_fasta_str(self, template_str):
    l = len(template_str)+self.line_width
    return '\n'.join([template_str[i:j] for i,j in zip(range(0,l,self.line_width), range(self.line_width,l,self.line_width))])
  def _modify_junction(self, d,offset_d):
      if offset_d is None:
        return ''
      else:
        spacing = max(d+offset_d,0)
        return 'N'*spacing

  def pred_to_fasta(self, fpath=None, offset_d=0):
    # split fasta file by 80 characters

    if fpath is None:
      out_f = sys.stdout
    else:
      out_f = open(fpath, 'wb')
      
    for bp_coord, (idx,bp) in sorted(self.bps_ordering.items(), key=itemgetter(1)):
      if not bp_coord in self.bps_set: continue
      a,b,d = bp_coord
      
      if a<0:
        pref = self._forward_extract(abs(a)-self.k)
      else:
        pref = self._reverse_extract(a)

      if b>0:
        suff = self._forward_extract(b)
      else:
        suff = self._reverse_extract(abs(b)-self.k)
      
      template_str = pref+self._modify_junction(d, offset_d) +suff
      print >>out_f, ">%s,%d,%d,%d\n%s"%(bp.name, a,b,d,self._format_fasta_str(template_str))
    if not fpath is None:
      out_f.close()
  
class PredictedMultiTemplates(PredictedTemplates):
  # Includes functionality based on "naming" of breakpoints
  def __init__(self, fpath, ref_fpath, contig):
    PredictedTemplates.__init__(self, fpath, ref_fpath, contig)
    self._create_multi_sets()
  
  def _filtering(self, tokens, bp):
    return tokens[2]>self.max_d or abs(abs(tokens[0])-abs(tokens[1]))<self.min_sv_dist or bp.name==''
  def _create_multi_sets(self):
    self.bp_multi_sets = defaultdict(list)
    # assumes bp names are sortable by template position.
    for bp_coord, (idx, bp) in sorted(self.bps_ordering.items(), key=lambda x: x[1][1].name):
      if not bp_coord in self.bps_set: continue
      
      self.bp_multi_sets[bp.name.split('_')[0]].append((bp_coord, bp))
      
  def pred_to_fasta(self, fpath=None, offset_d=0):
    if fpath is None:
      out_f = sys.stdout
    else:
      out_f = open(fpath, 'wb')
      
    for bp_name, l_bp_template in sorted(self.bp_multi_sets.items()):
      t = l_bp_template[:]
      
      bp_coord,bp = t[0]
      a,b,d = bp_coord
      if a<0:
        false_start = ((None, abs(a)-self.k, 0), None)
      else:
        false_start = ((None, -a-self.k, 0), None)
      
      bp_coord,bp = t[-1]
      a,b,d = bp_coord
      if b<0:
        false_end = ((abs(b)-self.k, None, 0), None)
      else:
        false_end = ((-(b+self.k), None, 0), None)
        
      t.append(false_end)
      t.insert(0, false_start)
      
      bp_coords, bps = zip(*t)
      starts, ends, ds = zip(*bp_coords)
      template_str = []
      for i,j,d in zip(ends[:-1],starts[1:],ds[1:]):
        assert i*j<0
        
        if i>j:
          seg = self._forward_seg_extract(i, abs(j))
        else:
          seg = self._reverse_seg_extract(j, abs(i))
          
        template_str.append(seg)
        template_str.append(self._modify_junction(d, offset_d))
      
      bp_rename = '+'.join(['%s,%d,%d,%d'%(bp.name, a,b,d) for (a,b,d),bp in t[1:-1]])
      print >>out_f, ">%s\n%s"%(bp_rename, self._format_fasta_str(''.join(template_str)))
      
    if not fpath is None:
      out_f.close()  
