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

# @file pyre/meshio/OutputSoln.py
##
# @brief Python object for managing output of finite-element
# solution information.
##
# Factory: output_manager

from .OutputManager import OutputManager
from .meshio import OutputSoln as ModuleOutputSoln

# OutputSoln class


class OutputSoln(OutputManager, ModuleOutputSoln):
    """
    Python object for managing output of finite-element solution
    information.

    @class Inventory
    Python object for managing OutputSoln facilities and properties.

    \b Properties
    @li \b vertex_data_fields Names of vertex data fields to output.

    \b Facilities
    @li None

    Factory: output_manager
    """

    # INVENTORY //////////////////////////////////////////////////////////

    import pyre.inventory

    vertexDataFields = pyre.inventory.list("vertex_data_fields", default=["all"])
    vertexDataFields.meta['tip'] = "Names of vertex data fields to output."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="OutputSoln"):
        """
        Constructor.
        """
        OutputManager.__init__(self, name)
        return

    def preinitialize(self):
        """
        Do
        """
        OutputManager.preinitialize(self)
        ModuleOutputSoln.vertexDataFields(self, self.vertexDataFields)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        OutputManager._configure(self)
        self.vertexDataFields = self.inventory.vertexDataFields
        return

    def _createModuleObj(self):
        """
        Create handle to C++ object.
        """
        ModuleOutputSoln.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def output_manager():
    """
    Factory associated with OutputManager.
    """
    return OutputSoln()


# End of file
