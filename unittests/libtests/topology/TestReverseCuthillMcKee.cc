// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestReverseCuthillMcKee.hh" // Implementation of class methods

#include "pylith/topology/ReverseCuthillMcKee.hh" // USES ReverseCuthillMcKee

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveStub.hh" // USES FaultCohesiveStub
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::topology::TestReverseCuthillMcKee::setUp(void)
{ // setUp
    PYLITH_METHOD_BEGIN;

    _data = new TestReverseCuthillMcKee_Data; CPPUNIT_ASSERT(_data);
    _mesh = NULL;

    PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::topology::TestReverseCuthillMcKee::tearDown(void)
{ // tearDown
    PYLITH_METHOD_BEGIN;

    delete _data; _data = NULL;
    delete _mesh; _mesh = NULL;

    PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test reorder().
void
pylith::topology::TestReverseCuthillMcKee::testReorder(void)
{ // testReorder
    PYLITH_METHOD_BEGIN;

    _initialize();
    CPPUNIT_ASSERT(_mesh);

    // Get original DM and create Mesh for it
    const PetscDM dmOrig = _mesh->dmMesh();
    PetscObjectReference((PetscObject) dmOrig);
    Mesh meshOrig;
    meshOrig.dmMesh(dmOrig);

    ReverseCuthillMcKee::reorder(_mesh);

    const PetscDM& dmMesh = _mesh->dmMesh(); CPPUNIT_ASSERT(dmMesh);

    // Check vertices (size only)
    topology::Stratum verticesStratumE(dmOrig, topology::Stratum::DEPTH, 0);
    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    CPPUNIT_ASSERT_EQUAL(verticesStratumE.size(), verticesStratum.size());

    // Check cells (size only)
    topology::Stratum cellsStratumE(dmOrig, topology::Stratum::HEIGHT, 0);
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    CPPUNIT_ASSERT_EQUAL(cellsStratumE.size(), cellsStratum.size());

    // Check groups
    PetscInt numGroupsE, numGroups;
    PetscErrorCode err;
    err = DMGetNumLabels(dmOrig, &numGroupsE); CPPUNIT_ASSERT(!err);
    err = DMGetNumLabels(dmMesh, &numGroups); CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL(numGroupsE, numGroups);

    for (PetscInt iGroup = 0; iGroup < numGroups; ++iGroup) {
        const char *name = NULL;
        err = DMGetLabelName(dmMesh, iGroup, &name); CPPUNIT_ASSERT(!err);

        PetscInt numPointsE, numPoints;
        err = DMGetStratumSize(dmOrig, name, 1, &numPointsE); CPPUNIT_ASSERT(!err);
        err = DMGetStratumSize(dmMesh, name, 1, &numPoints); CPPUNIT_ASSERT(!err);
        CPPUNIT_ASSERT_EQUAL(numPointsE, numPoints);
    } // for

    // Check element centroids
    PylithScalar coordsCheckOrig = 0.0;
    { // original
        Stratum cellsStratum(dmOrig, Stratum::HEIGHT, 0);
        const PetscInt cStart = cellsStratum.begin();
        const PetscInt cEnd = cellsStratum.end();
        topology::CoordsVisitor coordsVisitor(dmOrig);
        for (PetscInt cell = cStart; cell < cEnd; ++cell) {
            PetscScalar* coordsCell = NULL;
            PetscInt coordsSize = 0;
            PylithScalar value = 0.0;
            coordsVisitor.getClosure(&coordsCell, &coordsSize, cell);
            for (int i = 0; i < coordsSize; ++i) {
                value += coordsCell[i];
            } // for
            coordsCheckOrig += value*value;
            coordsVisitor.restoreClosure(&coordsCell, &coordsSize, cell);
        } // for
    } // original
    PylithScalar coordsCheck = 0.0;
    { // reordered
        Stratum cellsStratum(dmMesh, Stratum::HEIGHT, 0);
        const PetscInt cStart = cellsStratum.begin();
        const PetscInt cEnd = cellsStratum.end();
        topology::CoordsVisitor coordsVisitor(dmMesh);
        for (PetscInt cell = cStart; cell < cEnd; ++cell) {
            PetscScalar* coordsCell = NULL;
            PetscInt coordsSize = 0;
            PylithScalar value = 0.0;
            coordsVisitor.getClosure(&coordsCell, &coordsSize, cell);
            for (int i = 0; i < coordsSize; ++i) {
                value += coordsCell[i];
            } // for
            coordsCheck += value*value;
            coordsVisitor.restoreClosure(&coordsCell, &coordsSize, cell);
        } // for
    } // reordered
    const PylithScalar tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(coordsCheckOrig, coordsCheck, tolerance*coordsCheckOrig);

    // Verify reduction in Jacobian bandwidth
    Field fieldOrig(meshOrig);
    Field::Description description;
    description.label = "solution";
    description.vectorFieldType = FieldBase::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = "field";
    description.scale = 1.0;
    description.validator = NULL;

    Field::Discretization discretization;
    discretization.basisOrder = 1;
    discretization.quadOrder = 1;
    discretization.isBasisContinuous = true;
    discretization.feSpace = FieldBase::POLYNOMIAL_SPACE;
    fieldOrig.subfieldAdd(description, discretization);
    fieldOrig.subfieldsSetup();
    fieldOrig.allocate();
    PetscMat matrix = NULL;
    PetscInt bandwidthOrig = 0;
    err = DMCreateMatrix(fieldOrig.dmMesh(), &matrix); CPPUNIT_ASSERT(!err);
    err = MatComputeBandwidth(matrix, 0.0, &bandwidthOrig); CPPUNIT_ASSERT(!err);
    err = MatDestroy(&matrix); CPPUNIT_ASSERT(!err);

    Field field(*_mesh);
    field.subfieldAdd(description, discretization);
    field.subfieldsSetup();
    field.allocate();
    PetscInt bandwidth = 0;
    err = DMCreateMatrix(field.dmMesh(), &matrix); CPPUNIT_ASSERT(!err);
    err = MatComputeBandwidth(matrix, 0.0, &bandwidth); CPPUNIT_ASSERT(!err);
    err = MatDestroy(&matrix); CPPUNIT_ASSERT(!err);

    CPPUNIT_ASSERT(bandwidthOrig > 0);
    CPPUNIT_ASSERT(bandwidth > 0);
    CPPUNIT_ASSERT(bandwidth <= bandwidthOrig);

    PYLITH_METHOD_END;
} // testReorder


// ----------------------------------------------------------------------
void
pylith::topology::TestReverseCuthillMcKee::_initialize()
{ // _initialize
    PYLITH_METHOD_BEGIN;
    CPPUNIT_ASSERT(_data);

    delete _mesh; _mesh = new Mesh; CPPUNIT_ASSERT(_mesh);

    meshio::MeshIOAscii iohandler;
    iohandler.filename(_data->filename);
    iohandler.read(_mesh);
    CPPUNIT_ASSERT(_mesh->numCells() > 0);
    CPPUNIT_ASSERT(_mesh->numVertices() > 0);

    // Adjust topology if necessary.
    if (_data->faultLabel) {
        pylith::faults::FaultCohesiveStub fault;
        fault.id(100);
        fault.label(_data->faultLabel);
        fault.adjustTopology(_mesh);
    } // if

    PYLITH_METHOD_END;
} // _initialize

// ----------------------------------------------------------------------
// Constructor
pylith::topology::TestReverseCuthillMcKee_Data::TestReverseCuthillMcKee_Data(void) :
    filename(NULL),
    faultLabel(NULL)
{   // constructor
}   // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::topology::TestReverseCuthillMcKee_Data::~TestReverseCuthillMcKee_Data(void)
{   // destructor
}   // destructor


// End of file
