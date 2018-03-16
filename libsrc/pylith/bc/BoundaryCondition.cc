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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "BoundaryCondition.hh" // implementation of object methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cstring> // USES strlen()
#include <stdexcept> \
    // USES std::runtime_error()

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::BoundaryCondition::BoundaryCondition(void) :
    _label("")
{ // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::bc::BoundaryCondition::~BoundaryCondition(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::BoundaryCondition::deallocate(void) {}


// ----------------------------------------------------------------------
// Set mesh label associated with boundary condition surface.
void
pylith::bc::BoundaryCondition::label(const char* value) {
    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for boundary condition label.");
    } // if

    _label = value;
} // label


// ----------------------------------------------------------------------
// Get mesh label associated with boundary condition surface.
const char*
pylith::bc::BoundaryCondition::label(void) const {
    return _label.c_str();
} // label


// ----------------------------------------------------------------------
// Set name of field in solution to constrain.
void
pylith::bc::BoundaryCondition::field(const char* value) {
    PYLITH_METHOD_BEGIN;

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for name of solution field for boundary condition.");
    } // if
    _field = value;

    PYLITH_METHOD_END;
}  // field


// ----------------------------------------------------------------------
// Get name of field in solution to constrain.
const char*
pylith::bc::BoundaryCondition::field(void) const {
    journal::debug_t debug("boundarycondition");
    debug << journal::at(__HERE__)
          << "BoundaryCondition::field()" << journal::endl;

    return _field.c_str();
} // field


// End of file
