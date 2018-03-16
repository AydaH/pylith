#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

# @file pylith/bc/Dirichlet.py
##
# @brief Python object for managing a Dirichlet (prescribed
# values) boundary condition.
##
# Factory: boundary_condition

from .BoundaryCondition import BoundaryCondition
from pylith.feassemble.ConstraintPointwise import ConstraintPointwise

# Dirichlet class


class Dirichlet(BoundaryCondition,
                ConstraintPointwise):
    """
    Python object for managing a Dirichlet (prescribed values)
    boundary condition.

    Factory: boundary_condition
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="dirichletnew"):
        """
        Constructor.
        """
        BoundaryCondition.__init__(self, name)
        ConstraintPointwise.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Do pre-initialization setup.
        """
        BoundaryCondition.preinitialize(self, mesh)
        ConstraintPointwise.preinitialize(self, mesh)
        return

    def finalize(self):
        """
        Cleanup after running problem.
        """
        print(":TODO: @brad Implement Dirichlet.finalize() once output manager is added.")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        try:
            BoundaryCondition._configure(self)
            ConstraintPointwise._configure(self)
        except ValueError as err:
            aliases = ", ".join(self.aliases)
            raise ValueError("Error while configuring Dirichlet boundary condition "
                             "(%s):\n%s" % (aliases, err.message))
        return

# End of file
