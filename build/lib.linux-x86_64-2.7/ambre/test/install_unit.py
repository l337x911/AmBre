'''
#  ambre.test.valid.py
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
import os
from ambre.config import CONFIG
from ambre.test.analyze_unit import TestAnalyze
from ambre.test.design_unit import TestDesignWorkflow


class TestInstallationValid(unittest.TestCase):
  def setUp(self):
    self.config = CONFIG
  def testTmpDir(self):
    try:
      path = os.path.join(CONFIG.dir['default_temp'], 'writable.txt')
      with open(path, 'wb') as f:
        print >>f, "written"
      os.remove(path)
    except:
      self.fail("Writable temp directory is not configured. %s"%CONFIG.dir['default_temp'])
    
  def testPrimer3Installation(self):
    self.assertTrue(os.path.isdir(self.config.dir['primer3']),
                         msg="Primer3 directory not found.")
    self.assertTrue(os.path.isfile(os.path.join(self.config.dir['primer3'],'src','primer3_core')),
                         msg="Primer3 binary not found. Check compatible primer3 version.")
  def testMultiPlxInstallation(self):
    self.assertTrue(os.path.isfile(self.config.bin['multiplx']),
                         msg="Multiplx binary not found.\n")
  def testBlatInstallation(self):
    self.assertTrue(os.path.isfile(self.config.bin['aligner']),
                         msg="Blat binary not found.\n")
      
def install_check():
  import argparse
  
  parser = argparse.ArgumentParser(prog='install_unit.py', description='Unit tests to check installation. Runs small example for AmBre-design and AmBre-analyze')
  parser.add_argument('--examples-dir', type=str, nargs='?', default=None, dest="examples_dir", help='Specify directory to output example data')
  parser.add_argument('--config', type=str, nargs='?', default=None, dest="config_fpath", help='Update parameters in default config file with new config file')
  
  args = parser.parse_args()
  if not args.config_fpath is None:
    CONFIG.update(args.config_fpath[0])
  if not args.examples_dir is None:
    CONFIG.dir['examples'] = args.examples_dir[0]
    
  ic_suite = unittest.TestLoader().loadTestsFromTestCase(TestInstallationValid)
  ic_suite.addTest(TestDesignWorkflow('test_workflow'))
  ic_suite.addTest(TestAnalyze('test_workflow'))
  
  unittest.TextTestRunner(verbosity=2).run(ic_suite)
  
if __name__ == "__main__":
  install_check()