// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "OutputSoln.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <typeinfo> // USES typeid()

// ----------------------------------------------------------------------
const char* pylith::meshio::OutputSoln::_pyreComponent = "outputsoln";

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSoln::OutputSoln(void)
{ // constructor
    PyreComponent::name(_pyreComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSoln::~OutputSoln(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputSoln::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    OutputManager::deallocate();

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set names of solution fields to output.
void
pylith::meshio::OutputSoln::vertexDataFields(const char* names[],
                                             const int numNames)
{ // vertexDataFields
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputSoln::vertexDataFields(names="<<names<<", numNames="<<numNames<<")");

    assert((names && numNames) || (!names && !numNames));

    _vertexDataFields.resize(numNames);
    for (int i = 0; i < numNames; ++i) {
        assert(names[i]);
        _vertexDataFields[i] = names[i];
    } // for

    PYLITH_METHOD_END;
} // vertexDataFields

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::meshio::OutputSoln::verifyConfiguration(const pylith::topology::Field& solution,
                                                const pylith::topology::Field& auxField) const
{ // verifyConfiguration
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputSoln::verifyConfiguration(solution="<<solution.label()<<")");

    const size_t numFields = _vertexDataFields.size();
    if ((numFields > 0) && (std::string("all") != _vertexDataFields[0])) {
        for (size_t iField = 0; iField < numFields; iField++) {
            if (!solution.hasSubfield(_vertexDataFields[iField].c_str())) {
                std::ostringstream msg;
                msg << "Could not find field '" << _vertexDataFields[iField] << "' in solution '" << solution.label() << "' for output.";
                throw std::runtime_error(msg.str());
            } // if
        } // for
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration

// ----------------------------------------------------------------------
// Write solution at time step.
void
pylith::meshio::OutputSoln::writeTimeStep(const PylithReal t,
                                          const PylithInt tindex,
                                          const pylith::topology::Field& solution,
                                          const pylith::topology::Field& auxField)
{ // writeTimeStep
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputSoln::writeTimeStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.label()<<")");

    if (!this->shouldWrite(t, tindex)) {
        PYLITH_METHOD_END;
    } // if

    const pylith::string_vector& subfieldNames = (1 == _vertexDataFields.size() && std::string("all") == _vertexDataFields[0]) ? solution.subfieldNames() : _vertexDataFields;

    this->openTimeStep(t, solution.mesh());
    const size_t numFields = subfieldNames.size();
    for (size_t iField = 0; iField < numFields; iField++) {
        if (!solution.hasSubfield(subfieldNames[iField].c_str())) {
            std::ostringstream msg;
            msg << "Could not find field '" << subfieldNames[iField] << "' in solution for output.";
            throw std::runtime_error(msg.str());
        } // if

        pylith::topology::Field& fieldBuffer = this->getBuffer(solution, subfieldNames[iField].c_str());
        this->appendVertexField(t, fieldBuffer, fieldBuffer.mesh());
    } // for
    this->closeTimeStep();

    PYLITH_METHOD_END;
} // writeTimeStep


// End of file
