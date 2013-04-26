'''
#  ambre.test.design_unit.py
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
from ambre.design.workflow import PrimerDesignWorkflow
from pkg_resources import resource_filename
import unittest
import os

class TestDesignWorkflow(unittest.TestCase):
  def setUp(self):
    pass
  
  def test_workflow(self):
    CONFIG.param['reference_fpath'] = resource_filename('ambre', os.path.join('examples', 'reference.fasta'))

    w = PrimerDesignWorkflow(primer3_path=CONFIG.dir['primer3'],
                     primer3_param=CONFIG.param['primer3_long'], 
                     aligner=CONFIG.bin['aligner'], 
                     multiplx=CONFIG.bin['multiplx'],
                     d=1000,rho=0.1)
    ex_regions = resource_filename('ambre', os.path.join('examples','regions.test'))
    if CONFIG.dir['examples'] is None:
      ex_temp = resource_filename('ambre', os.path.join('examples','regions_ex', 'design_ex'))
    else:
      ex_temp = os.path.join(CONFIG.dir['examples'], 'design_ex')
    
    os.link(ex_regions, os.path.join(CONFIG.dir['examples'], os.path.basename(ex_regions)))
    os.link(CONFIG.param['reference_fpath'], os.path.join(CONFIG.dir['examples'], os.path.basename(CONFIG.param['reference_fpath'])))
    
    w.run(ex_regions,temp_tag=ex_temp,
        primers_per_bp=75, max_cross_amp_dist=20000,
        max_primer_penalty=1.5,
        max_alignment_count=10,
        min_alignment_len=18,
        pamp_max_iterations= 1000000,pamp_repeats=2,
        pamp_t_ms= (-(10**-1),), pamp_t_bs=(10**4,))
    w.print_solutions(out_fpath = "%s.out"%ex_temp)
    #w.validate()
  
  def test_workflow_valid(self):
    CONFIG.param['reference_fpath'] = resource_filename('ambre', os.path.join('examples', 'reference.fasta'))

    w = PrimerDesignWorkflow(primer3_path=CONFIG.dir['primer3'],
                     primer3_param=CONFIG.param['primer3_long'], 
                     aligner=CONFIG.bin['aligner'], 
                     multiplx=CONFIG.bin['multiplx'],
                     d=1000,rho=0.1)
    ex_regions = resource_filename('ambre', os.path.join('examples','regions.test'))
    ex_temp = resource_filename('ambre', os.path.join('examples','regions_ex', 'ex'))
    w.check(ex_regions,temp_tag=ex_temp)

if __name__ == '__main__':
  suite = unittest.TestLoader().loadTestsFromTestCase(TestDesignWorkflow)
  unittest.TextTestRunner(verbosity=2).run(suite)
