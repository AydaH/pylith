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

## @file unittests/feassemble/testfeassemble.py

## @brief Python application for testing feassemble code.

from pylith.tests.UnitTestApp import UnitTestApp

import unittest

class TestApp(UnitTestApp):
  """
  Test application.
  """

  def __init__(self):
    """
    Constructor.
    """
    UnitTestApp.__init__(self)
    return


  def _suite(self):
    """
    Setup the test suite.
    """

    suite = unittest.TestSuite()

    from TestCellGeometry import TestCellGeometry
    suite.addTest(unittest.makeSuite(TestCellGeometry))

    from TestFIATSimplex import TestFIATSimplex
    suite.addTest(unittest.makeSuite(TestFIATSimplex))

    from TestFIATLagrange import TestFIATLagrange
    suite.addTest(unittest.makeSuite(TestFIATLagrange))

    from TestMeshQuadrature import TestMeshQuadrature
    suite.addTest(unittest.makeSuite(TestMeshQuadrature))

    from TestElasticityImplicit import TestElasticityImplicit
    suite.addTest(unittest.makeSuite(TestElasticityImplicit))

    from TestElasticityExplicit import TestElasticityExplicit
    suite.addTest(unittest.makeSuite(TestElasticityExplicit))

    from TestElasticityImplicitLgDeform import TestElasticityImplicitLgDeform
    suite.addTest(unittest.makeSuite(TestElasticityImplicitLgDeform))

    from TestElasticityExplicitLgDeform import TestElasticityExplicitLgDeform
    suite.addTest(unittest.makeSuite(TestElasticityExplicitLgDeform))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
