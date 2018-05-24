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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

# @file pylith/problems/Problem.py
##
# @brief Python abstract base class for crustal dynamics problems.
##
# Factory: problem.

from pylith.utils.PetscComponent import PetscComponent
from .problems import Problem as ModuleProblem

from pylith.utils.NullComponent import NullComponent
from pylith.meshio.OutputSoln import OutputSoln
from pylith.feassemble.IntegratorPointwise import IntegratorPointwise
from pylith.feassemble.ConstraintPointwise import ConstraintPointwise

# ITEM FACTORIES ///////////////////////////////////////////////////////


def materialFactory(name):
    """
    Factory for material items.
    """
    from pyre.inventory import facility
    from pylith.materials.IsotropicLinearElasticityPlaneStrain import IsotropicLinearElasticityPlaneStrain
    return facility(name, family="material", factory=IsotropicLinearElasticityPlaneStrain)


def bcFactory(name):
    """
    Factory for boundary condition items.
    """
    from pyre.inventory import facility
    from pylith.bc.DirichletTimeDependent import DirichletTimeDependent
    return facility(name, family="boundary_condition", factory=DirichletTimeDependent)


def faultFactory(name):
    """
    Factory for fault items.
    """
    from pyre.inventory import facility
    from pylith.faults.FaultCohesiveKin import FaultCohesiveKin
    return facility(name, family="fault", factory=FaultCohesiveKin)


def outputFactory(name):
    """
    Factory for output items.
    """
    from pyre.inventory import facility
    return facility(name, family="output_manager", factory=OutputSoln)


class Problem(PetscComponent, ModuleProblem):
    """
    Python abstract base class for crustal dynamics problems.

    Factory: problem.
    """

    # INVENTORY //////////////////////////////////////////////////////////
    #
    # @class Inventory
    # Python object for managing Problem facilities and properties.
    #
    # \b Properties
    # @li \b dimension Spatial dimension of problem space.
    # @li \b solver Type of solver to use.
    #
    # \b Facilities
    # @li \b solution Solution field.
    # @li \b normalizer Nondimensionalizer for problem.
    # @li \b materials Array of materials (governing equations) in the problem.
    # @li \b bc Array of boundary conditions.
    # @li \b interfaces Array of interior surfaces with relative displacement constraints or constitutive models.
    # @li \b solution_outputs Array of output managers for solution.
    # @li \b gravityField Gravity field for problem (SpatialDB).

    import pyre.inventory
    from pylith.utils.EmptyBin import EmptyBin

    solverChoice = pyre.inventory.str("solver", default="linear", validator=pyre.inventory.choice(["linear", "nonlinear"]))
    solverChoice.meta['tip'] = "Type of solver to use ['linear', 'nonlinear']."

    from Solution import Solution
    solution = pyre.inventory.facility("solution", family="solution", factory=Solution)
    solution.meta['tip'] = "Solution field for problem."

    from spatialdata.units.NondimElasticQuasistatic import NondimElasticQuasistatic
    normalizer = pyre.inventory.facility("normalizer", family="nondimensional", factory=NondimElasticQuasistatic)
    normalizer.meta['tip'] = "Nondimensionalizer for problem."

    from pylith.materials.Homogeneous import Homogeneous
    materials = pyre.inventory.facilityArray("materials", itemFactory=materialFactory, factory=Homogeneous)
    materials.meta['tip'] = "Materials in problem."

    bc = pyre.inventory.facilityArray("bc", itemFactory=bcFactory, factory=EmptyBin)
    bc.meta['tip'] = "Boundary conditions."

    interfaces = pyre.inventory.facilityArray("interfaces", itemFactory=faultFactory, factory=EmptyBin)
    interfaces.meta['tip'] = "Interior surfaces with constraints or constitutive models."

    from pylith.meshio.SingleOutput import SingleOutput
    outputs = pyre.inventory.facilityArray("solution_outputs", itemFactory=outputFactory, factory=SingleOutput)
    outputs.meta['tip'] = "Output managers for solution."

    gravityField = pyre.inventory.facility("gravity_field", family="spatial_database", factory=NullComponent)
    gravityField.meta['tip'] = "Database used for gravity field."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="problem"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="problem")
        self.mesh = None
        return

    def preinitialize(self, mesh):
        """
        Do minimal initialization.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Performing minimal initialization before verifying configuration.")

        ModuleProblem.identifier(self, self.aliases[-1])
        ModuleProblem.solverType(self, self.solverType)
        ModuleProblem.normalizer(self, self.normalizer)
        if not isinstance(self.gravityField, NullComponent):
            ModuleProblem.gravityField(self, self.gravityField)

        # Do minimal setup of solution.
        self.solution.preinitialize(mesh, self.normalizer)
        ModuleProblem.solution(self, self.solution.field)

        # Preinitialize materials
        for material in self.materials.components():
            material.preinitialize(mesh)

        # Preinitialize BC
        for bc in self.bc.components():
            bc.preinitialize(mesh)

        # Preinitialize interfaces
        for interface in self.interfaces.components():
            interface.preinitialize(mesh)

        # Preinitialize solution output.
        for output in self.outputs.components():
            output.preinitialize()

        # Set integrators and constraints.
        self._setIntegratorsConstraints()

        # Set solution outputs.
        self._setSolutionOutputs()

        ModuleProblem.preinitialize(self, mesh)
        return

    def verifyConfiguration(self):
        """
        Verify compatibility of configuration.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Verifying compatibility of problem configuration.")

        # Check to make sure ids of materials and interfaces are unique
        materialIds = {}
        for material in self.materials.components():
            if material.id() in materialIds.keys():
                raise ValueError("ID values for materials '%s' and '%s' are both '%d'. Material id values must be unique." %
                                 (material.label, materialIds[material.materialId], material.materialId))
            materialIds[material.materialId] = material.label

        for interface in self.interfaces.components():
            if interface.id() in materialIds.keys():
                raise ValueError("ID values for material '%s' and interface '%s' are both '%d'. Material and interface id values must be unique." %
                                 (materialIds[interface.id()], interface.label(), interface.id()))
            materialIds[interface.id()] = interface.label()

        # Check to make sure material-id for each cell matches the id of a material
        import numpy
        idValues = numpy.array(materialIds.keys(), dtype=numpy.int32)
        ModuleProblem.verifyConfiguration(self, idValues)

        self._printInfo()
        return

    def initialize(self):
        """
        Initialize integrators and constraints.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Initializing problem.")

        ModuleProblem.initialize(self)
        return

    def run(self, app):
        """
        Solve the problem.
        """
        raise NotImplementedError("run() not implemented.")
        return

    def finalize(self):
        """
        Cleanup after running problem.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Finalizing problem.")

        # Close solution output.
        for output in self.outputs.components():
            output.close()
        return

    def checkpoint(self):
        """
        Save problem state for restart.
        """
        raise NotImplementedError("checkpoint() not implemented.")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        PetscComponent._configure(self)

        if self.solverChoice == "linear":
            self.solverType = ModuleProblem.LINEAR
        elif self.solverChoice == "nonlinear":
            self.solverType = ModuleProblem.NONLINEAR
        else:
            raise ValueError("Unknown solver choice '%s'." % self.solverChoice)

        return

    def _setIntegratorsConstraints(self):
        integrators = []
        constraints = []

        for material in self.materials.components():
            if not isinstance(material, IntegratorPointwise):
                raise TypeError("Material '%s' fails integrator implementation test." % material.name)
            integrators.append(material)

        for interface in self.interfaces.components():
            if not isinstance(interface, IntegratorPointwise):
                raise TypeError("Interface '%s' fails integrator implementation test." % interface.name)
            integrators.append(interface)

        for bc in self.bc.components():
            if isinstance(bc, IntegratorPointwise):
                integrators.append(bc)
            elif isinstance(bc, ConstraintPointwise):
                constraints.append(bc)
            else:
                raise TypeError("Unable to classify bc '%s' into an in integrator or constraint." % bc)

        ModuleProblem.integrators(self, integrators)
        ModuleProblem.constraints(self, constraints)
        return

    def _setSolutionOutputs(self):
        ModuleProblem.outputs(self, self.outputs.components())
        return

    def _printInfo(self):
        """Write overview of problem to info journal.
        """
        msg = (
            "Scales for nondimensionalization:",
            "    Length scale: {}".format(self.normalizer.lengthScale()),
            "    Time scale: {}".format(self.normalizer.timeScale()),
            "    Pressure scale: {}".format(self.normalizer.pressureScale()),
            "    Density scale: {}".format(self.normalizer.densityScale()),
            "    Temperature scale: {}".format(self.normalizer.temperatureScale()),
        )
        self._info.log("\n".join(msg))
        return

    def _setupLogging(self):
        """
        Setup event logging.
        """
        if not "_loggingPrefix" in dir(self):
            self._loggingPrefix = ""

        from pylith.utils.EventLogger import EventLogger
        logger = EventLogger()
        logger.className("Problem")
        logger.initialize()

        self._eventLogger = logger
        return


# FACTORIES ////////////////////////////////////////////////////////////

def problem():
    """
    Factory associated with Problem.
    """
    return Problem()


# End of file
