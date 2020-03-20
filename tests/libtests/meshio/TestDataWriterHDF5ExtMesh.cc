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

#include "TestDataWriterHDF5ExtMesh.hh" // Implementation of class methods

#include "pylith/utils/types.hh" // HASA PylithScalar

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/DataWriterHDF5Ext.hh" // USES DataWriterHDF5Ext

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMesh::setUp(void)
{ // setUp
    PYLITH_METHOD_BEGIN;

    TestDataWriterMesh::setUp();
    _data = NULL;

    PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterHDF5ExtMesh::tearDown(void)
{ // tearDown
    PYLITH_METHOD_BEGIN;

    TestDataWriterMesh::tearDown();
    delete _data; _data = NULL;

    PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testConstructor(void)
{ // testConstructor
    PYLITH_METHOD_BEGIN;

    DataWriterHDF5Ext writer;

    CPPUNIT_ASSERT(writer._h5);

    PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testFilename(void)
{ // testDebug
    PYLITH_METHOD_BEGIN;

    DataWriterHDF5Ext writer;

    const char* filename = "data.h5";
    writer.filename(filename);
    CPPUNIT_ASSERT_EQUAL(std::string(filename), writer._filename);

    PYLITH_METHOD_END;
} // testFilename

// ----------------------------------------------------------------------
// Test open() and close()
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testOpenClose(void)
{ // testOpenClose
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5Ext writer;

    writer.filename(_data->opencloseFilename);

    const bool isInfo = false;
    writer.open(*_mesh, isInfo);
    writer.close();

    checkFile(_data->opencloseFilename);

    PYLITH_METHOD_END;
} // testOpenClose

// ----------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testWriteVertexField(void)
{ // testWriteVertexField
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5Ext writer;

    topology::Fields vertexFields(*_mesh);
    _createVertexFields(&vertexFields);

    writer.filename(_data->vertexFilename);

    const PylithScalar timeScale = 4.0;
    writer.setTimeScale(timeScale);
    const PylithScalar t = _data->time / timeScale;

    const bool isInfo = false;
    writer.open(*_mesh, isInfo);
    writer.openTimeStep(t, *_mesh);

    const int numFields = 4;
    const char* fieldNames[4] = {"scalar", "vector", "tensor", "other"};
    for (int i = 0; i < numFields; ++i) {
        pylith::topology::Field& field = vertexFields.get(fieldNames[i]);
        writer.writeVertexField(t, field, *_mesh);
    } // for

    writer.closeTimeStep();
    writer.close();

    checkFile(_data->vertexFilename);

    PYLITH_METHOD_END;
} // testWriteVertexField

// ----------------------------------------------------------------------
// Test writeCellField.
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testWriteCellField(void)
{ // testWriteCellField
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    DataWriterHDF5Ext writer;

    topology::Fields cellFields(*_mesh);
    _createCellFields(&cellFields);

    writer.filename(_data->cellFilename);

    const PylithScalar timeScale = 4.0;
    writer.setTimeScale(timeScale);
    const PylithScalar t = _data->time / timeScale;

    const bool isInfo = false;
    writer.open(*_mesh, isInfo);
    writer.openTimeStep(t, *_mesh);

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
// Test hdf5Filename().
void pylith::meshio::TestDataWriterHDF5ExtMesh::testHdf5Filename(void)
{ // testHdf5Filename
  PYLITH_METHOD_BEGIN;

  DataWriterHDF5Ext writer;

  // Append info to filename if number of time steps is 0.
  writer._isInfo = true;
  writer._filename = "output.h5";
  CPPUNIT_ASSERT_EQUAL(std::string("output_info.h5"), writer.hdf5Filename());
		       
  writer._isInfo = false;
  writer._filename = "output_abc.h5";
  CPPUNIT_ASSERT_EQUAL(std::string("output_abc.h5"), writer.hdf5Filename());
  
  writer._isInfo = false;
  writer._filename = "output_abcd.h5";
  CPPUNIT_ASSERT_EQUAL(std::string("output_abcd.h5"), writer.hdf5Filename());

  PYLITH_METHOD_END;
} // testHdf5Filename

// ----------------------------------------------------------------------
// Test _datasetFilename().
void
pylith::meshio::TestDataWriterHDF5ExtMesh::testDatasetFilename(void)
{ // testDatasetFilename
    PYLITH_METHOD_BEGIN;

    DataWriterHDF5Ext writer;

    // Append info to filename if info.
    writer._isInfo = true;
    writer._filename = "output.h5";
    CPPUNIT_ASSERT_EQUAL(std::string("output_info_ABCD.dat"), writer._datasetFilename("ABCD"));

    writer._isInfo = false;
    writer._filename = "output_abc.h5";
    CPPUNIT_ASSERT_EQUAL(std::string("output_abc_field1.dat"), writer._datasetFilename("field1"));

    writer._isInfo = false;
    writer._filename = "output_abcd.h5";
    CPPUNIT_ASSERT_EQUAL(std::string("output_abcd_field2.dat"), writer._datasetFilename("field2"));

    PYLITH_METHOD_END;
} // testDatasetFilename


// ----------------------------------------------------------------------
// Get test data.
pylith::meshio::TestDataWriter_Data*
pylith::meshio::TestDataWriterHDF5ExtMesh::_getData(void)
{ // _getData
    return _data;
} // _getData


// End of file
