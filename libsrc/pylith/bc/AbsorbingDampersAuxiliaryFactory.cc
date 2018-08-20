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

#include "AbsorbingDampersAuxiliaryFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // HOLDSA AuxiliaryField
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::AbsorbingDampersAuxiliaryFactory::AbsorbingDampersAuxiliaryFactory(void) {
    GenericComponent::name("absorbingdampersauxiliaryfactory");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::AbsorbingDampersAuxiliaryFactory::~AbsorbingDampersAuxiliaryFactory(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add density field to auxiliary fields.
void
pylith::bc::AbsorbingDampersAuxiliaryFactory::addDensity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addDensity(void)");

    const char* fieldName = "density";
    const PylithReal densityScale = _normalizer->densityScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = densityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addDensity


// ---------------------------------------------------------------------------------------------------------------------
// Add shear wave speed field to auxiliary fields.
void
pylith::bc::AbsorbingDampersAuxiliaryFactory::addVs(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addVs(void)");

    const char* fieldName = "vs";
    const PylithReal velocityScale = _normalizer->lengthScale() / _normalizer->timeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = velocityScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addVs


// ---------------------------------------------------------------------------------------------------------------------
// Add dilatational wave speed field to auxiliary fields.
void
pylith::bc::AbsorbingDampersAuxiliaryFactory::addVp(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addVp(void)");

    const char* fieldName = "vp";
    const PylithReal velocityScale = _normalizer->lengthScale() / _normalizer->timeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = velocityScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addVp


// End of file
