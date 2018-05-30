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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "FieldQuery.hh" // implementation of class methods

#include "Field.hh" // USES Field
#include "Mesh.hh" // USES Mesh

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include "pylith/utils/error.hh" \
    // USES PYLITH_CHECK_ERROR

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::FieldQuery::FieldQuery(const Field& field) :
    _field(field),
    _functions(NULL),
    _contexts(NULL),
    _contextPtrs(NULL)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::FieldQuery::~FieldQuery(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::FieldQuery::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    delete[] _functions; _functions = NULL;
    delete[] _contexts; _contexts = NULL;
    delete[] _contextPtrs; _contextPtrs = NULL;

    _queryFns.clear();
    _queryDBs.clear();

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Set query function information for subfield.
void
pylith::topology::FieldQuery::queryFn(const char* subfield,
                                      const queryfn_type fn,
                                      spatialdata::spatialdb::SpatialDB* db)
{ // queryFn
    PYLITH_METHOD_BEGIN;

    assert(subfield);

    _queryFns[subfield] = fn;
    _queryDBs[subfield] = db;

    PYLITH_METHOD_END;
} // queryFn


// ----------------------------------------------------------------------
// Get query function information for subfield.
const pylith::topology::FieldQuery::queryfn_type
pylith::topology::FieldQuery::queryFn(const char* subfield) const
{ // queryFn
    PYLITH_METHOD_BEGIN;

    assert(subfield);
    const queryfn_map_type::const_iterator& iter = _queryFns.find(subfield);
    if (iter == _queryFns.end()) {
        std::ostringstream msg;
        msg << "Could not find query function for subfield '" << subfield << "'." << std::endl;
        throw std::logic_error(msg.str());
    } // if

    PYLITH_METHOD_RETURN(*iter->second);
} // queryFn

// ----------------------------------------------------------------------
// Initialize query with default query functions.
void
pylith::topology::FieldQuery::initializeWithDefaultQueryFns(void) {
    PYLITH_METHOD_BEGIN;

    _queryFns.clear();
    _queryDBs.clear();

    const pylith::string_vector& subfields = _field.subfieldNames();
    const unsigned int numSubfields = subfields.size();
    for (unsigned int i = 0; i < numSubfields; ++i) {
        _queryFns[subfields[i]] = dbQueryGeneric;
        _queryDBs[subfields[i]] = NULL;
    } // for

    PYLITH_METHOD_END;
} // initializeWithDefaultQueryFns


// ----------------------------------------------------------------------
// Get spatial database used to get values for subfield.
const spatialdata::spatialdb::SpatialDB*
pylith::topology::FieldQuery::queryDB(const char* subfield) const
{ // queryDB
    PYLITH_METHOD_BEGIN;

    assert(subfield);
    const querydb_map_type::const_iterator& iter = _queryDBs.find(subfield);
    if (iter == _queryDBs.end()) {
        std::ostringstream msg;
        msg << "Could not find spatial database for subfield '" << subfield << "'." << std::endl;
        throw std::logic_error(msg.str());
    } // if

    PYLITH_METHOD_RETURN(iter->second);
} // queryDB


// ----------------------------------------------------------------------
// Get array of query functions.
pylith::topology::FieldQuery::queryfn_type*
pylith::topology::FieldQuery::functions(void) const
{ // functions
    return _functions;
} // functions


// ----------------------------------------------------------------------
// Get array of pointers to contexts.
const pylith::topology::FieldQuery::DBQueryContext* const*
pylith::topology::FieldQuery::contextPtrs(void) const
{ // contextPtrs
    return _contextPtrs;
} // contextPtrs

// ----------------------------------------------------------------------
// Query spatial database to set values in field.
void
pylith::topology::FieldQuery::openDB(spatialdata::spatialdb::SpatialDB* db,
                                     const PylithReal lengthScale)
{ // openDB
    PYLITH_METHOD_BEGIN;

    assert(db);

    // Create contexts and funcs. Need to put contexts into an array of
    // pointers, since Petsc function doesn't know the size of the
    // context.
    const Field::subfields_type& subfields = _field._subfields;
    const unsigned size = subfields.size();
    assert(_queryFns.size() == size);
    delete[] _functions; _functions = (size > 0) ? new queryfn_type[size] : NULL;
    delete[] _contexts; _contexts = (size > 0) ? new DBQueryContext[size] : NULL;
    delete[] _contextPtrs; _contextPtrs = (size > 0) ? new DBQueryContext*[size] : NULL;

    int i = 0;
    for (Field::subfields_type::const_iterator iter = subfields.begin(); iter != subfields.end(); ++iter, ++i) {
        const std::string& name = iter->first;
        if (_queryFns.find(name) == _queryFns.end()) { // if
            std::ostringstream msg;
            msg << "FieldQuery for field '" << _field.label() << "' missing query function for subfield '" << name << "'";
            throw std::logic_error(msg.str());
        } // if/else
        const PylithInt index = iter->second.index;
        assert(size_t(index) < subfields.size());

        _functions[index] = _queryFns[name];

        _contexts[index].db = (_queryDBs[name]) ? _queryDBs[name] : db;
        _contexts[index].cs = _field.mesh().coordsys();
        _contexts[index].lengthScale = lengthScale;

        const pylith::topology::Field::Description& description = iter->second.description;
        _contexts[index].valueScale = description.scale;
        _contexts[index].description = description.label;
        _contexts[index].componentNames = description.componentNames;
        _contexts[index].validator = description.validator;

        _contextPtrs[index] = &_contexts[index];
    } // for

    // Open spatial database.
    db->open();

    PYLITH_METHOD_END;
} // openDB


// ----------------------------------------------------------------------
// Query spatial database to set values in field.
void
pylith::topology::FieldQuery::queryDB(void)
{ // queryDB
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = 0;
    PetscReal dummyTime = 0.0;
    err = DMProjectFunctionLocal(_field.dmMesh(), dummyTime, _functions, (void**)_contextPtrs, INSERT_ALL_VALUES, _field.localVector()); PYLITH_CHECK_ERROR(err);
    //err = PetscObjectCompose((PetscObject) dm, "A", (PetscObject) fieldVec);CHKERRQ(ierr); // :MATT: Which dm is this? Do we need this?

    PYLITH_METHOD_END;
} // queryDB


// ----------------------------------------------------------------------
// Query spatial database to set values in field.
void
pylith::topology::FieldQuery::closeDB(spatialdata::spatialdb::SpatialDB* db)
{ // queryDB
    PYLITH_METHOD_BEGIN;

    assert(db);

    delete[] _functions; _functions = NULL;
    delete[] _contexts; _contexts = NULL;
    delete[] _contextPtrs; _contextPtrs = NULL;

    // Close spatial database.
    db->close();

    PYLITH_METHOD_END;
} // queryDB


// ----------------------------------------------------------------------
// Generic query of values from spatial database.
PetscErrorCode
pylith::topology::FieldQuery::dbQueryGeneric(PylithInt dim,
                                             PylithReal t,
                                             const PylithReal x[],
                                             PylithInt nvalues,
                                             PylithScalar* values,
                                             void* context)
{ // dbQueryPositive
    PYLITH_METHOD_BEGIN;

    assert(x);
    assert(values);
    assert(context);

    const pylith::topology::FieldQuery::DBQueryContext* queryctx = (pylith::topology::FieldQuery::DBQueryContext*)context; assert(queryctx);

    const int numQueryValues = queryctx->componentNames.size();
    assert(numQueryValues == nvalues);

    // Tell database which values we want.
    const char** queryValueNames = (numQueryValues > 0) ? new const char*[numQueryValues] : NULL;
    for (int i = 0; i < numQueryValues; ++i) {
        queryValueNames[i] = queryctx->componentNames[i].c_str();
    } // for
    try {
        queryctx->db->queryVals(queryValueNames, numQueryValues);
    } catch (const std::runtime_error& err) {
        delete[] queryValueNames; queryValueNames = NULL;
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, err.what());
        PYLITH_METHOD_RETURN(1);
    } // try/catch

    // Dimensionalize query location coordinates.
    assert(queryctx->lengthScale > 0);
    double xDim[3];
    for (int i = 0; i < dim; ++i) {
        xDim[i] = x[i] * queryctx->lengthScale;
    } // for

    // Query database.
    assert(queryctx->cs);
    const int err = queryctx->db->query(values, nvalues, xDim, dim, queryctx->cs);

    if (err) {
        delete[] queryValueNames; queryValueNames = NULL;
        std::ostringstream msg;
        msg << "Could not find " << queryctx->description << " at (";
        for (int i = 0; i < dim; ++i)
            msg << "  " << xDim[i];
        msg << ") in spatial database '" << queryctx->db->label() << "'.";
        PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
        PYLITH_METHOD_RETURN(1);
    } // if

    // Validate if validator function was specified.
    if (queryctx->validator) {
        for (int i = 0; i < nvalues; ++i) {
            const std::string& invalidMsg = queryctx->validator(values[i]);
            if (invalidMsg.length() > 0) {
                delete[] queryValueNames; queryValueNames = NULL;

                std::ostringstream msg;
                msg << "Found invalid " << queryValueNames[i] << " (" << values[i] << ") at location (";
                for (int i = 0; i < dim; ++i)
                    msg << "  " << xDim[i];
                msg << ") in spatial database '" << queryctx->db->label() << "'. ";
                msg << invalidMsg;
                PYLITH_SET_ERROR(PETSC_COMM_SELF, PETSC_ERR_LIB, msg.str().c_str());
                PYLITH_METHOD_RETURN(1);
            } // if
        } // for
    } // if

    // Nondimensionalize values
    assert(queryctx->valueScale > 0);
    for (int i = 0; i < nvalues; ++i) {
        values[i] /= queryctx->valueScale;
    } // for

    delete[] queryValueNames; queryValueNames = NULL;

    PYLITH_METHOD_RETURN(0);
} // dbQueryGeneric


// ----------------------------------------------------------------------
const char*
pylith::topology::FieldQuery::validatorPositive(const PylithReal value)
{ // validatorPositive
    return (value > 0.0) ? "" : "Value must be positive.";
} // validatorPositive


// End of file
