// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "Problem.hh" // implementation of class methods

#include "pylith/problems/IntegrationData.hh" // HOLDSA IntegrationData
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // HASA Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/materials/Material.hh" // USES Material
#include "pylith/faults/FaultCohesive.hh" // USES FaultCohesive
#include "pylith/bc/BoundaryCondition.hh" // USES BoundaryCondition
#include "pylith/feassemble/Integrator.hh" // USES Integrator
#include "pylith/feassemble/Constraint.hh" // USES Constraint
#include "pylith/problems/ObserversSoln.hh" // USES ObserversSoln
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor::optimizeClosure()

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Problem::Problem() :
    _integrationData(new pylith::problems::IntegrationData),
    _normalizer(NULL),
    _gravityField(NULL),
    _observers(new pylith::problems::ObserversSoln),
    _formulation(pylith::problems::Physics::QUASISTATIC),
    _solverType(LINEAR) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::Problem::~Problem(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::Problem::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        delete _integrators[i];_integrators[i] = NULL;
    } // for

    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        delete _constraints[i];_constraints[i] = NULL;
    } // for

    delete _integrationData;_integrationData = NULL;
    delete _normalizer;_normalizer = NULL;
    _gravityField = NULL; // Held by Python. :KLUDGE: :TODO: Use shared pointer.
    delete _observers;_observers = NULL;

    pylith::topology::FieldOps::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set formulation for solving equation.
void
pylith::problems::Problem::setFormulation(const pylith::problems::Physics::FormulationEnum value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setFormulation(value="<<value<<")");

    _formulation = value;

    PYLITH_METHOD_END;
} // setFormulation


// ---------------------------------------------------------------------------------------------------------------------
// Get formulation for solving equation.
pylith::problems::Physics::FormulationEnum
pylith::problems::Problem::getFormulation(void) const {
    return _formulation;
} // getFormulation


// ---------------------------------------------------------------------------------------------------------------------
// Set problem type.
void
pylith::problems::Problem::setSolverType(const SolverTypeEnum value) {
    PYLITH_COMPONENT_DEBUG("Problem::setSolverType(value="<<value<<")");

    _solverType = value;
} // setSolverType


// ---------------------------------------------------------------------------------------------------------------------
// Get problem type.
pylith::problems::Problem::SolverTypeEnum
pylith::problems::Problem::getSolverType(void) const {
    return _solverType;
} // getSolverType


// ---------------------------------------------------------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::problems::Problem::setNormalizer(const spatialdata::units::Nondimensional& dim) {
    PYLITH_COMPONENT_DEBUG("Problem::setNormalizer(dim="<<typeid(dim).name()<<")");

    if (!_normalizer) {
        _normalizer = new spatialdata::units::Nondimensional(dim);
    } else {
        *_normalizer = dim;
    } // if/else
} // setNormalizer


// ---------------------------------------------------------------------------------------------------------------------
// Set gravity field.
void
pylith::problems::Problem::setGravityField(spatialdata::spatialdb::GravityField* const g) {
    PYLITH_COMPONENT_DEBUG("Problem::setGravityField(g="<<typeid(g).name()<<")");

    _gravityField = g;
} // setGravityField


// ----------------------------------------------------------------------
// Register observer to receive notifications.
void
pylith::problems::Problem::registerObserver(pylith::problems::ObserverSoln* observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("registerObserver(observer="<<typeid(observer).name()<<")");

    assert(_observers);
    assert(_normalizer);
    _observers->registerObserver(observer);
    _observers->setTimeScale(_normalizer->getTimeScale());

    PYLITH_METHOD_END;
} // registerObserver


// ----------------------------------------------------------------------
// Remove observer from receiving notifications.
void
pylith::problems::Problem::removeObserver(pylith::problems::ObserverSoln* observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("removeObserver(observer="<<typeid(observer).name()<<")");

    assert(_observers);
    _observers->removeObserver(observer);

    PYLITH_METHOD_END;
} // removeObserver


// ---------------------------------------------------------------------------------------------------------------------
// Set solution field.
void
pylith::problems::Problem::setSolution(pylith::topology::Field* field) {
    PYLITH_COMPONENT_DEBUG("Problem::setSolution(field="<<typeid(*field).name()<<")");

    assert(_integrationData);
    _integrationData->setField("solution", field);
} // setSolution


// ---------------------------------------------------------------------------------------------------------------------
// Get solution field.
const pylith::topology::Field*
pylith::problems::Problem::getSolution(void) const {
    PYLITH_METHOD_BEGIN;

    assert(_integrationData);
    pylith::topology::Field* solution = NULL;
    if (_integrationData->hasField(IntegrationData::solution)) {
        solution = _integrationData->getField(IntegrationData::solution);
    } // if

    PYLITH_METHOD_RETURN(solution);
}


// ---------------------------------------------------------------------------------------------------------------------
// Set materials.
void
pylith::problems::Problem::setMaterials(pylith::materials::Material* materials[],
                                        const int numMaterials) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::setMaterials("<<materials<<", numMaterials="<<numMaterials<<")");

    assert( (!materials && 0 == numMaterials) || (materials && 0 < numMaterials) );

    _materials.resize(numMaterials);
    for (int i = 0; i < numMaterials; ++i) {
        _materials[i] = materials[i];
    } // for

    PYLITH_METHOD_END;
} // setMaterials


// ---------------------------------------------------------------------------------------------------------------------
// Set boundary conditions.
void
pylith::problems::Problem::setBoundaryConditions(pylith::bc::BoundaryCondition* bc[],
                                                 const int numBC) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::setBoundaryConditions("<<bc<<", numBC="<<numBC<<")");

    assert( (!bc && 0 == numBC) || (bc && 0 < numBC) );

    _bc.resize(numBC);
    for (int i = 0; i < numBC; ++i) {
        _bc[i] = bc[i];
    } // for

    PYLITH_METHOD_END;
} // setBoundaryConditions


// ---------------------------------------------------------------------------------------------------------------------
// Set materials.
void
pylith::problems::Problem::setInterfaces(pylith::faults::FaultCohesive* interfaces[],
                                         const int numInterfaces) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::setInterfaces("<<interfaces<<", numInterfaces="<<numInterfaces<<")");

    assert( (!interfaces && 0 == numInterfaces) || (interfaces && 0 < numInterfaces) );

    _interfaces.resize(numInterfaces);
    for (int i = 0; i < numInterfaces; ++i) {
        _interfaces[i] = interfaces[i];
    } // for

    PYLITH_METHOD_END;
} // setInterfaces


// ----------------------------------------------------------------------
// Do minimal initialization.
void
pylith::problems::Problem::preinitialize(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::preinitialzie(mesh="<<typeid(mesh).name()<<")");

    assert(_normalizer);

    const size_t numMaterials = _materials.size();
    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        _materials[i]->setNormalizer(*_normalizer);
        _materials[i]->setGravityField(_gravityField);
        _materials[i]->setFormulation(_formulation);
    } // for

    const size_t numInterfaces = _interfaces.size();
    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        _interfaces[i]->setNormalizer(*_normalizer);
        _interfaces[i]->setFormulation(_formulation);
    } // for

    const size_t numBC = _bc.size();
    for (size_t i = 0; i < numBC; ++i) {
        assert(_bc[i]);
        _bc[i]->setNormalizer(*_normalizer);
        _bc[i]->setFormulation(_formulation);
    } // for

    PYLITH_METHOD_END;
} // preinitialize


// ----------------------------------------------------------------------
// Verify configuration.
void
pylith::problems::Problem::verifyConfiguration(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::verifyConfiguration(void)");

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);

    // Check to make sure materials are compatible with the solution.
    const size_t numMaterials = _materials.size();
    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        _materials[i]->verifyConfiguration(*solution);
    } // for

    // Check to make sure interfaces are compatible with the solution.
    const size_t numInterfaces = _interfaces.size();
    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        _interfaces[i]->verifyConfiguration(*solution);
    } // for

    // Check to make sure boundary conditions are compatible with the solution.
    const size_t numBC = _bc.size();
    for (size_t i = 0; i < numBC; ++i) {
        assert(_bc[i]);
        _bc[i]->verifyConfiguration(*solution);
    } // for

    _checkMaterialIds();

    assert(_observers);
    _observers->verifyObservers(*solution);

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Initialize.
void
pylith::problems::Problem::initialize(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::initialize()");

    assert(_integrationData);
    pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);

    // Initialize solution field.
    PetscErrorCode err = DMSetFromOptions(solution->getDM());PYLITH_CHECK_ERROR(err);
    _setupSolution();

    const pylith::topology::Mesh& mesh = solution->getMesh();
    pylith::topology::CoordsVisitor::optimizeClosure(mesh.getDM());

    // Initialize integrators.
    _createIntegrators();
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        assert(_integrators[i]);
        _integrators[i]->initialize(*solution);
    } // for

    // Initialize constraints.
    _createConstraints();
    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        assert(_constraints[i]);
        _constraints[i]->initialize(*solution);
    } // for

    solution->allocate();
    solution->createGlobalVector();
    solution->createOutputVector();

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG("Displaying solution field layout");
        solution->view("Solution field", pylith::topology::Field::VIEW_LAYOUT);
    } // if

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Check material and interface ids.
void
pylith::problems::Problem::_checkMaterialIds(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::_checkMaterialIds()");

    const size_t numMaterials = _materials.size();
    const size_t numInterfaces = _interfaces.size();

    pylith::int_array materialIds(numMaterials + numInterfaces);
    size_t count = 0;
    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        materialIds[count++] = _materials[i]->getMaterialId();
    } // for
    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        materialIds[count++] = _interfaces[i]->getInterfaceId();
    } // for

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);
    pylith::topology::MeshOps::checkMaterialIds(solution->getMesh(), materialIds);

    PYLITH_METHOD_END;
} // _checkMaterialIds


// ---------------------------------------------------------------------------------------------------------------------
// Create array of integrators from materials, interfaces, and boundary conditions.
void
pylith::problems::Problem::_createIntegrators(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::_createIntegrators()");

    const size_t numMaterials = _materials.size();
    const size_t numInterfaces = _interfaces.size();
    const size_t numBC = _bc.size();

    const size_t maxSize = numMaterials + numInterfaces + numBC;
    _integrators.resize(maxSize);
    size_t count = 0;

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);

    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        pylith::feassemble::Integrator* integrator = _materials[i]->createIntegrator(*solution);
        assert(count < maxSize);
        if (integrator) { _integrators[count++] = integrator;}
    } // for

    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        pylith::feassemble::Integrator* integrator = _interfaces[i]->createIntegrator(*solution);
        assert(count < maxSize);
        if (integrator) { _integrators[count++] = integrator;}
    } // for

    // Check to make sure boundary conditions are compatible with the solution.
    for (size_t i = 0; i < numBC; ++i) {
        assert(_bc[i]);
        pylith::feassemble::Integrator* integrator = _bc[i]->createIntegrator(*solution);
        assert(count < maxSize);
        if (integrator) { _integrators[count++] = integrator;}
    } // for

    _integrators.resize(count);

    PYLITH_METHOD_END;
} // _createIntegrators


// ---------------------------------------------------------------------------------------------------------------------
// Create array of constraints from materials, interfaces, and boundary conditions.
void
pylith::problems::Problem::_createConstraints(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::_createConstraints()");

    const size_t numMaterials = _materials.size();
    const size_t numInterfaces = _interfaces.size();
    const size_t numBC = _bc.size();

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);

    _constraints.resize(0); // insure we start with an empty array.

    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        std::vector<pylith::feassemble::Constraint*> constraints = _materials[i]->createConstraints(*solution);
        _constraints.insert(_constraints.end(), constraints.begin(), constraints.end());

    } // for

    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        std::vector<pylith::feassemble::Constraint*> constraints = _interfaces[i]->createConstraints(*solution);
        _constraints.insert(_constraints.end(), constraints.begin(), constraints.end());

    } // for

    for (size_t i = 0; i < numBC; ++i) {
        assert(_bc[i]);
        std::vector<pylith::feassemble::Constraint*> constraints = _bc[i]->createConstraints(*solution);
        _constraints.insert(_constraints.end(), constraints.begin(), constraints.end());

    } // for

    PYLITH_METHOD_END;
} // _createConstraints


// ---------------------------------------------------------------------------------------------------------------------
// Setup solution subfields and discretization.
void
pylith::problems::Problem::_setupSolution(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::_setupSolution()");

    assert(_integrationData);
    pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);
    solution->subfieldsSetup();
    solution->createDiscretization();

    // Mark fault fields as implicit.
    const pylith::string_vector& subfieldNames = solution->getSubfieldNames();
    for (size_t i = 0; i < subfieldNames.size(); ++i) {
        const pylith::topology::Field::SubfieldInfo& subfieldInfo = solution->getSubfieldInfo(subfieldNames[i].c_str());
        if (subfieldInfo.fe.isFaultOnly) {
            PetscErrorCode err;
            PetscDS ds = NULL;
            PetscInt cStart = 0, cEnd = 0;
            PetscDM dmSoln = solution->getDM();assert(dmSoln);
            err = DMPlexGetHeightStratum(dmSoln, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
            PetscInt cell = cStart;
            for (; cell < cEnd; ++cell) {
                if (pylith::topology::MeshOps::isCohesiveCell(dmSoln, cell)) { break; }
            } // for
            err = DMGetCellDS(dmSoln, cell, &ds);PYLITH_CHECK_ERROR(err);
            assert(ds);
            err = PetscDSSetImplicit(ds, subfieldInfo.index, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
        } // if
    } // for

    PYLITH_METHOD_END;
} // _setupSolution


// End of file
