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
    self.talignment_fpath = resource_filename('ambre',  os.path.join('examples', 'aligned_blasr_h5.sam.talign'))
    self.bp_fpath = resource_filename('ambre',  os.path.join('examples', 'aligned_blasr_h5.sam.bp'))
    self.template_fpath = resource_filename('ambre',  os.path.join('examples', 'aligned_blasr_h5.sam.templates.fasta'))

  def test_filtering(self):
    worker = MaxAlignmentScoreFragFiltering(self.alignment_fpath)
    worker.frag_isoform_filtering(self.talignment_fpath)
    assert os.path.isfile(self.talignment_fpath)
  def test_clustering(self):
    assembler = AssemblyTAlign(self.talignment_fpath)
    bps = assembler.condense_fragments_clustering()
    assembler.output_breakpoints(bps, breakpoint_fpath=self.bp_fpath)
    assert os.path.isfile(self.bp_fpath)
    
  def test_template_generating(self):
    ref_fpath = '/home/anand/data/hg/hg19/full.fa' 
    # currently only supports SVs on one contig
    worker = PredictedTemplates(self.bp_fpath, ref_fpath, 'chr9')
    worker.pred_to_fasta(self.template_fpath, None)
    assert os.path.isfile(self.template_fpath)
  def test_workflow(self):
    ambre_run('/home/anand/data/hg/hg19/full.fa','chr9',self.alignment_fpath)

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