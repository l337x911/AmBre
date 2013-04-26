'''
#  ambre.test.analyze_unit.py
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
import unittest
from pkg_resources import resource_filename
import os
from ambre.analyze.align_seg import MaxAlignmentScoreFragFiltering
from ambre.analyze.frag_clustering import AssemblyTAlign
from ambre.analyze.predict_templates import PredictedTemplates
from ambre.analyze.workflow import ambre_run


class TestAnalyze(unittest.TestCase):
  def setUp(self):
    self.alignment_fpath = resource_filename('ambre',  os.path.join('examples', 'aligned_blasr_h5.sam'))
    self.ref_fpath = resource_filename('ambre',  os.path.join('examples', 'reference.fasta'))
    
    os.link(self.alignment_fpath, os.path.join(CONFIG.dir['examples'], os.path.basename(self.alignment_fpath)))
    os.link(self.ref_fpath, os.path.join(CONFIG.dir['examples'], os.path.basename(self.ref_fpath)))
    
    if CONFIG.dir['examples'] is None:
      self.ex_temp = resource_filename('ambre', os.path.join('examples','regions_ex', 'analyze_ex'))
    else:
      self.ex_temp = os.path.join(CONFIG.dir['examples'], 'analyze_ex')
    
    self.talignment_fpath = "%s.talign"%self.ex_temp
    self.bp_fpath = '%s.bp'%self.ex_temp
    self.template_fpath = '%s.fasta'%self.ex_temp

  def test_filtering(self):
    worker = MaxAlignmentScoreFragFiltering(self.alignment_fpath)
    worker.frag_isoform_filtering(self.talignment_fpath)
    self.assertTrue(os.path.isfile(self.talignment_fpath), msg='Filtering file creation failed')
  def test_clustering(self):
    assembler = AssemblyTAlign(self.talignment_fpath)
    bps = assembler.condense_fragments_clustering()
    assembler.output_breakpoints(bps, breakpoint_fpath=self.bp_fpath)
    self.assertTrue(os.path.isfile(self.bp_fpath), msg='Clustering bp file creation failed')
    
  def test_template_generating(self): 
    # currently only supports SVs on one contig
    worker = PredictedTemplates(self.bp_fpath, self.ref_fpath, 'chr9:21815000-22135000')
    worker.pred_to_fasta(self.template_fpath, None)
    assert os.path.isfile(self.template_fpath)
  def test_workflow(self):
    CONFIG.param['analyze_min_cluster_size'] = 25
    ambre_run(self.ref_fpath,'chr9:21815000-22135000',self.alignment_fpath,
              temp_tag=self.ex_temp)
    self.assertTrue(os.path.isfile(self.talignment_fpath),msg='Filtered alignments file %s not created.'%self.talignment_fpath)
    self.assertTrue(os.path.isfile(self.bp_fpath),msg='BP file %s not created.'%self.bp_fpath)
    self.assertTrue(os.path.isfile(self.template_fpath),msg='Template fasta file %s not created.'%self.template_fpath)

def three_step_pipeline_suite():
  suite = unittest.TestSuite()
  suite.addTest(TestAnalyze('test_filtering'))
  suite.addTest(TestAnalyze('test_clustering'))
  suite.addTest(TestAnalyze('test_template_generating'))
  return suite

def workflow_suite():
  suite = unittest.TestSuite()
  suite.addTest(TestAnalyze('test_workflow'))
  return suite
if __name__ == "__main__":
  #unittest.TextTestRunner(verbosity=2).run(three_step_pipeline_suite())
  unittest.TextTestRunner(verbosity=2).run(workflow_suite())