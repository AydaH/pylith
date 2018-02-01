#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

from pylith.tests.FullTestApp import FullTestApp

import unittest

class TestApp(FullTestApp):
  """
  Test application.
  """

  def __init__(self):
    """
    Constructor.
    """
    FullTestApp.__init__(self)
    return


  def _suite(self):
    """
    Create test suite.
    """
    suite = unittest.TestSuite()

    from TestPetscApp import TestPetscApp
    suite.addTest(unittest.makeSuite(TestPetscApp))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.main()

  
# End of file 
