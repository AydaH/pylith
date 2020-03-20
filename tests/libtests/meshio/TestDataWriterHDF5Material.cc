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

#include <portinfo>

#include "TestDataWriterHDF5Material.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/DataWriterHDF5.hh" // USES DataWriterHDF5

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5Material::setUp(void) {
    PYLITH_METHOD_BEGIN;

    TestDataWriterMaterial::setUp();
    _data = NULL;

    PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterHDF5Material::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    TestDataWriterMaterial::tearDown();
    delete _data;_data = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ----------------------------------------------------------------------
// Test open() and close()
void
pylith::meshio::TestDataWriterHDF5Material::testOpenClose(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_materialMesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5 writer;

    writer.filename(_data->opencloseFilename);

    const bool isInfo = false;
    writer.open(*_materialMesh, isInfo);
    writer.close();

    checkFile(_data->opencloseFilename);

    PYLITH_METHOD_END;
} // testOpenClose


// ----------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterHDF5Material::testWriteVertexField(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_domainMesh);
    CPPUNIT_ASSERT(_materialMesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5 writer;

    pylith::topology::Fields vertexFields(*_domainMesh);
    _createVertexFields(&vertexFields);

    writer.filename(_data->vertexFilename);

    const PylithScalar timeScale = 4.0;
    writer.setTimeScale(timeScale);
    const PylithScalar t = _data->time / timeScale;

    const bool isInfo = false;
    writer.open(*_materialMesh, isInfo);
    writer.openTimeStep(t, *_materialMesh);

    const int numFields = 4;
    const char* fieldNames[4] = {"scalar", "vector", "tensor", "other"};
    for (int i = 0; i < numFields; ++i) {
        pylith::topology::Field& field = vertexFields.get(fieldNames[i]);
        writer.writeVertexField(t, field, *_materialMesh);
    } // for

    writer.closeTimeStep();
    writer.close();

    checkFile(_data->vertexFilename);

    PYLITH_METHOD_END;
} // testWriteVertexField


// ----------------------------------------------------------------------
// Test writeCellField.
void
pylith::meshio::TestDataWriterHDF5Material::testWriteCellField(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_materialMesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5 writer;

    pylith::topology::Fields cellFields(*_materialMesh);
    _createCellFields(&cellFields);

    writer.filename(_data->cellFilename);

    const PylithScalar timeScale = 4.0;
    writer.setTimeScale(timeScale);
    const PylithScalar t = _data->time / timeScale;

    const bool isInfo = false;
    writer.open(*_materialMesh, isInfo);
    writer.openTimeStep(t, *_materialMesh);

    const int numFields = 4;
    const char* fieldNames[4] = {"scalar", "vector", "tensor", "other"};
    for (int i = 0; i < numFields; ++i) {
        pylith::topology::Field& field = cellFields.get(fieldNames[i]);
        writer.writeCellField(t, field);
    } // for

    writer.closeTimeStep();
    writer.close();

    checkFile(_data->cellFilename);

    PYLITH_METHOD_END;
} // testWriteCellField


// ----------------------------------------------------------------------
// Get test data.
pylith::meshio::TestDataWriterMaterial_Data*
pylith::meshio::TestDataWriterHDF5Material::_getData(void) { // _getData
    return _data;
} // _getData


// End of file
