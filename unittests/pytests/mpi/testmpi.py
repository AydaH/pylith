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

## @file unittests/bc/testmpi.py

## @brief Python application for testing mpi code.

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

    from TestCommunicator import TestCommunicator
    suite.addTest(unittest.makeSuite(TestCommunicator))

    from TestReduce import TestReduce
    suite.addTest(unittest.makeSuite(TestReduce))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
