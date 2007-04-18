// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFaultCohesive.hh" // Implementation of class methods

#include "data/CohesiveData1D.hh" // USES CohesiveData1D

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultsCohesiveKin

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh
#include "pylith/utils/array.hh" // USES int_array
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include <stdexcept> // TEMPORARY

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesive );

// ----------------------------------------------------------------------
// Test adjustTopology() with 1-D line element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyLine(void)
{ // testAdjustTopologyLine
  CohesiveData1D data;
  _testAdjustTopology(data);
} // testAdjustTopologyLine

// ----------------------------------------------------------------------
// Test adjustTopology().
void
pylith::faults::TestFaultCohesive::_testAdjustTopology(const CohesiveData& data)
{ // _testAdjustTopology
  ALE::Obj<ALE::Mesh> mesh;
  meshio::MeshIOAscii iohandler;
  iohandler.filename(data.filename);
  iohandler.debug(true);
  iohandler.interpolate(false);
  iohandler.read(&mesh);

  FaultCohesiveKin fault;
  fault.id(0);
  fault.label("fault");
  fault.adjustTopology(&mesh);

  mesh->view("Mesh");

  CPPUNIT_ASSERT_EQUAL(data.cellDim, mesh->getDimension());

  // Check vertices
  const ALE::Obj<Mesh::label_sequence>& vertices = mesh->depthStratum(0);
  const ALE::Obj<Mesh::real_section_type>& coordsField =
    mesh->getRealSection("coordinates");
  const int numVertices = vertices->size();
  CPPUNIT_ASSERT_EQUAL(data.numVertices, numVertices);
  CPPUNIT_ASSERT_EQUAL(data.spaceDim, 
		       coordsField->getFiberDimension(*vertices->begin()));
  int i = 0;
  const int spaceDim = data.spaceDim;
  for(Mesh::label_sequence::iterator v_iter = 
	vertices->begin();
      v_iter != vertices->end();
      ++v_iter) {
    const Mesh::real_section_type::value_type *vertexCoords = 
      coordsField->restrictPoint(*v_iter);
    const double tolerance = 1.0e-06;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      if (data.vertices[i] < 1.0) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(data.vertices[i++], vertexCoords[iDim],
				   tolerance);
      } else {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, vertexCoords[iDim]/data.vertices[i++],
				   tolerance);
      }
  } // for

  // check cells
  const ALE::Obj<sieve_type>& sieve = mesh->getSieve();
  const ALE::Obj<Mesh::label_sequence>& cells = mesh->heightStratum(0);

  const int numCells = cells->size();
  CPPUNIT_ASSERT_EQUAL(data.numCells, numCells);

  i = 0;
  const int offset = numCells;
  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter) {
    const int numCorners = sieve->nCone(*c_iter, mesh->depth())->size();
    CPPUNIT_ASSERT_EQUAL(data.numCorners[i], numCorners);
    const ALE::Obj<sieve_type::traits::coneSequence>& cone = 
      sieve->cone(*c_iter);
    for(sieve_type::traits::coneSequence::iterator v_iter = cone->begin();
	v_iter != cone->end();
	++v_iter)
      CPPUNIT_ASSERT_EQUAL(data.cells[i++], *v_iter - offset);
  } // for

  // check materials
  const ALE::Obj<Mesh::label_type>& labelMaterials = 
    mesh->getLabel("material-id");
  const int idDefault = -999;
  const int size = numCells;
  int_array materialIds(size);
  i = 0;
  for(Mesh::label_sequence::iterator c_iter = cells->begin();
      c_iter != cells->end();
      ++c_iter)
    materialIds[i++] = mesh->getValue(labelMaterials, *c_iter, idDefault);
  
  for (int iCell=0; iCell < numCells; ++iCell)
    CPPUNIT_ASSERT_EQUAL(data.materialIds[iCell], materialIds[iCell]);

  // Check groups
  const ALE::Obj<std::set<std::string> >& groupNames = 
    mesh->getIntSections();
  int iGroup = 0;
  int index = 0;
  for (std::set<std::string>::const_iterator name=groupNames->begin();
       name != groupNames->end();
       ++name, ++iGroup) {
    const ALE::Obj<int_section_type>& groupField = mesh->getIntSection(*name);
    CPPUNIT_ASSERT(!groupField.isNull());
    const int_section_type::chart_type& chart = groupField->getChart();
    const Mesh::point_type firstPoint = *chart.begin();
    std::string groupType = 
      (mesh->height(firstPoint) == 0) ? "cell" : "vertex";
    const int numPoints = chart.size();
    int_array points(numPoints);
    int i = 0;
    const int offset = ("vertex" == groupType) ? numCells : 0;
    for(int_section_type::chart_type::iterator c_iter = chart.begin();
	c_iter != chart.end();
	++c_iter)
      points[i++] = *c_iter - offset;

    CPPUNIT_ASSERT_EQUAL(std::string(data.groupNames[iGroup]), *name);
    CPPUNIT_ASSERT_EQUAL(std::string(data.groupTypes[iGroup]), groupType);
    CPPUNIT_ASSERT_EQUAL(data.groupSizes[iGroup], numPoints);
    for (int i=0; i < numPoints; ++i)
      CPPUNIT_ASSERT_EQUAL(data.groups[index++], points[i]);
  } // for
  
} // _testAdjustTopology

// End of file 
