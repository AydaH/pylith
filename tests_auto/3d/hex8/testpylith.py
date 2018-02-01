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

    from TestAxialDisp import TestAxialDisp
    suite.addTest(unittest.makeSuite(TestAxialDisp))

    from TestShearDisp import TestShearDisp
    suite.addTest(unittest.makeSuite(TestShearDisp))

    from TestShearDispNoSlip import TestShearDispNoSlip
    suite.addTest(unittest.makeSuite(TestShearDispNoSlip))

    from TestShearDispNoSlipRefine import TestShearDispNoSlipRefine
    suite.addTest(unittest.makeSuite(TestShearDispNoSlipRefine))

    from TestSlipOneFault import TestSlipOneFault
    suite.addTest(unittest.makeSuite(TestSlipOneFault))

    from TestSlipTwoFaults import TestSlipTwoFaults
    suite.addTest(unittest.makeSuite(TestSlipTwoFaults))

    from TestFaultsIntersect import TestFaultsIntersect
    suite.addTest(unittest.makeSuite(TestFaultsIntersect))

    from TestFrictionNoSlip import TestFrictionNoSlip
    suite.addTest(unittest.makeSuite(TestFrictionNoSlip))

    from TestFrictionNoSlipHalo import TestFrictionNoSlipHalo
    suite.addTest(unittest.makeSuite(TestFrictionNoSlipHalo))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.main()

  
# End of file 
