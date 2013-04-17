'''
#  ambre.analyze.workflow.py
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

import os
import tempfile
import time

from ambre.analyze.align_seg import MaxAlignmentScoreFragFiltering
from ambre.analyze.frag_clustering import AssemblyTAlign
from ambre.analyze.predict_templates import PredictedTemplates

class BreakpointAnalyzeWorkflow(object):
  '''
  classdocs
  '''
  def __init__(self, reference_fpath, contig):
    '''
    Constructor
    '''
    self.reference_fpath = reference_fpath
    self.contig = contig
  def run(self, sam_fpath, temp_tag=None, delete_flag=False):
    if temp_tag is None:
      out_dir=CONFIG.dir['default_temp']
      temp_tag = tempfile.mktemp(prefix='', dir=out_dir)
    elif os.path.isdir(temp_tag):
      out_dir=temp_tag
      temp_tag = tempfile.mktemp(prefix='', dir=out_dir)
    
    talign_fpath = "%s.talign"%temp_tag
    bp_fpath = "%s.bp"%temp_tag
    template_fpath = "%s.fasta"%temp_tag
    
    try:
      
      pre_filtering_time=time.time()
      if not os.path.isfile(talign_fpath):
        filtering = MaxAlignmentScoreFragFiltering(sam_fpath)
        filtering.frag_isoform_filtering(talign_fpath)
        assert os.path.isfile(talign_fpath)
      print "#Stage1 Filtering: %.4f"%(time.time()-pre_filtering_time)
        
      pre_clustering_time=time.time()
      if not os.path.isfile(bp_fpath):
        clustering = AssemblyTAlign(talign_fpath)
        try:
          max_d=int(CONFIG.param['analyze_cluster_max_d'])
        except KeyError:
          max_d=None
          
        bps = clustering.condense_fragments_clustering(max_d)
        clustering.output_breakpoints(bps,
                                      breakpoint_fpath=bp_fpath, 
                                      min_cluster_size=int(CONFIG.param['analyze_min_cluster_size']))
      print "#Stage2 Clustering: %.4f"%(time.time()-pre_clustering_time)
      
      pre_templating_time = time.time()
      templating = PredictedTemplates(bp_fpath, self.reference_fpath, self.contig)
      try:
        offset_d=int(CONFIG.param['analyze_template_offset_d'])
      except KeyError:
        offset_d=None
      templating.pred_to_fasta(template_fpath, offset_d)
      print "#Stage3 Templating: %.4f"%(time.time()-pre_templating_time)
        
    finally:
      if(delete_flag):
        if os.path.isfile(talign_fpath):
          os.remove(talign_fpath)
            
def ambre_run(reference_fpath, contig, sam_fpath, temp_tag=None):  
  w = BreakpointAnalyzeWorkflow(reference_fpath, contig)
  w.run(sam_fpath, temp_tag=temp_tag, delete_flag=(CONFIG.param['cleanup_flag']=="True"))

def main():
  import argparse
  
  parser = argparse.ArgumentParser(description='Analyzes local alignments on sequencing fragments to call breakpoints and form template sequences.')
  parser.add_argument('reference', type=str, nargs=1, help='Fasta file with multiple sequence entries.')
  parser.add_argument('contig', type=str, nargs=1, help='Fasta file with multiple sequence entries.')
  parser.add_argument('alignments', type=str, nargs=1, help='BLASR output in SAM format')
  
  parser.add_argument('temptag', type=str, nargs='?', default=None, help='Prefix for the run id / Directory to store Temps of a run id')
   
  args = parser.parse_args()
  ambre_run(args.reference[0], args.contig[0], args.alignments[0], args.temptag)
if __name__=='__main__':
  main()

