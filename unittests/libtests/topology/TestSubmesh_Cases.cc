// -*- C++ -*-
//
// -----------------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// -----------------------------------------------------------------------------
//

#include <portinfo>

#include "TestSubmesh.hh"   // Implementation of class methods

// -----------------------------------------------------------------------------
namespace pylith {
    namespace topology {

        // ---------------------------------------------------------------------
        class TestSubmesh_Tri : public TestSubmesh {

            CPPUNIT_TEST_SUB_SUITE( TestSubmesh_Tri, TestSubmesh );
            CPPUNIT_TEST_SUITE_END();

            void setUp(void) {
                TestSubmesh::setUp();

                _data->cellDim = 2;
                _data->numVertices = 4;
                _data->numCells = 2;
                _data->numCorners = 3;
                static const int _cells[2*3] = {
                    0, 1, 3,
                    0, 3, 2,
                };
                _data->cells = const_cast<int*>(_cells);
                static const PylithScalar _coordinates[4*2] = {
                    0.0, 0.0,
                    1.0, 0.0,
                    0.0, 1.0,
                    1.0, 1.0,
                };
                _data->coordinates = const_cast<PylithScalar*>(_coordinates);
                _data->label = "bc";
                _data->groupSize = 3;
                static const int _groupVertices[3] = {
                    1, 2, 3,
                };
                _data->groupVertices = const_cast<int*>(_groupVertices);
                _data->submeshNumCorners = 2;
                _data->submeshNumVertices = 3;
                static const int _submeshVertices[3] = {
                    2, 3, 4,
                };
                _data->submeshVertices = const_cast<int*>(_submeshVertices);
                _data->submeshNumCells = 2;
                static const int _submeshCells[2] = {
                    0, 1,
                };
                _data->submeshCells = const_cast<int*>(_submeshCells);
            }   // setUp

        };  // _TestSubmesh_Tri
        CPPUNIT_TEST_SUITE_REGISTRATION( TestSubmesh_Tri );

    }   // topology
}   // pylith


// End of file
