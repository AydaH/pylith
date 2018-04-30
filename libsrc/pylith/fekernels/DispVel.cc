/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University of Chicago
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2015 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/DispVel.hh"

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for time evolution equation with displacement and velocity
 * solution fields.
 *
 * Solution fields = [disp(dim), vel(dim, optional)]
 * Auxiliary fields = None
 *
 * \int_V \vec{\phi}_v \cdot \left( \frac{\partial \vec{u}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_v \cdot \vec{v}(t) \, dV.
 *
 * ======================================================================
 */

// ----------------------------------------------------------------------
// f0 entry function for disp/velocity equation.
void
pylith::fekernels::DispVel::f0u(const PylithInt dim,
                                const PylithInt numS,
                                const PylithInt numA,
                                const PylithInt sOff[],
                                const PylithInt sOff_x[],
                                const PylithScalar s[],
                                const PylithScalar s_t[],
                                const PylithScalar s_x[],
                                const PylithInt aOff[],
                                const PylithInt aOff_x[],
                                const PylithScalar a[],
                                const PylithScalar a_t[],
                                const PylithScalar a_x[],
                                const PylithReal t,
                                const PylithScalar x[],
                                const PylithInt numConstants,
                                const PylithScalar constants[],
                                PylithScalar f0[]) {
    const PylithInt _numS = 2;
    const PylithInt i_disp = 0;
    const PylithScalar* disp_t = &s_t[sOff[i_disp]];

    PylithInt i;

    assert(_numS == numS);
    assert(sOff);
    assert(s);
    assert(s_t);
    assert(f0);

    for (i = 0; i < dim; ++i) {
        f0[i] += disp_t[i];
    } // for
} // f0u


// ----------------------------------------------------------------------
// g0 function for disp/velocity equation.
void
pylith::fekernels::DispVel::g0u(const PylithInt dim,
                                const PylithInt numS,
                                const PylithInt numA,
                                const PylithInt sOff[],
                                const PylithInt sOff_x[],
                                const PylithScalar s[],
                                const PylithScalar s_t[],
                                const PylithScalar s_x[],
                                const PylithInt aOff[],
                                const PylithInt aOff_x[],
                                const PylithScalar a[],
                                const PylithScalar a_t[],
                                const PylithScalar a_x[],
                                const PylithReal t,
                                const PylithScalar x[],
                                const PylithInt numConstants,
                                const PylithScalar constants[],
                                PylithScalar g0[]) {
    const PylithInt _numS = 2;
    const PylithInt i_vel = 1;
    const PylithScalar* vel = &s[sOff[i_vel]];

    PylithInt i;

    assert(_numS == numS);
    assert(sOff);
    assert(s);
    assert(g0);

    for (i = 0; i < dim; ++i) {
        g0[i] += vel[i];
    } // for
} // g0u


// ----------------------------------------------------------------------
// Jf0 function for disp/velocity equation with zero values on diagonal.
void
pylith::fekernels::DispVel::Jf0uu_zero(const PylithInt dim,
                                       const PylithInt numS,
                                       const PylithInt numA,
                                       const PylithInt sOff[],
                                       const PylithInt sOff_x[],
                                       const PylithScalar s[],
                                       const PylithScalar s_t[],
                                       const PylithScalar s_x[],
                                       const PylithInt aOff[],
                                       const PylithInt aOff_x[],
                                       const PylithScalar a[],
                                       const PylithScalar a_t[],
                                       const PylithScalar a_x[],
                                       const PylithReal t,
                                       const PylithReal utshift,
                                       const PylithScalar x[],
                                       const PylithInt numConstants,
                                       const PylithScalar constants[],
                                       PylithScalar Jf0[]) {
    // No work to do for zero values.
} //


// ----------------------------------------------------------------------
// Jf0 function for disp/velocity equation with implicit time-stepping.
void
pylith::fekernels::DispVel::Jf0uu(const PylithInt dim,
                                  const PylithInt numS,
                                  const PylithInt numA,
                                  const PylithInt sOff[],
                                  const PylithInt sOff_x[],
                                  const PylithScalar s[],
                                  const PylithScalar s_t[],
                                  const PylithScalar s_x[],
                                  const PylithInt aOff[],
                                  const PylithInt aOff_x[],
                                  const PylithScalar a[],
                                  const PylithScalar a_t[],
                                  const PylithScalar a_x[],
                                  const PylithReal t,
                                  const PylithReal utshift,
                                  const PylithScalar x[],
                                  const PylithInt numConstants,
                                  const PylithScalar constants[],
                                  PylithScalar Jf0[]) {
    const PylithInt _numS = 2;

    PylithInt i;

    assert(_numS == numS);

    for (i = 0; i < dim; ++i) {
        Jf0[i*dim+i] += utshift;
    } // for
} // Jf0uu


// ----------------------------------------------------------------------
/* Jg0 function for disp/velocity equation.
 *
 * Solution fields = [disp(dim), vel(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::DispVel::Jg0uv(const PylithInt dim,
                                  const PylithInt numS,
                                  const PylithInt numA,
                                  const PylithInt sOff[],
                                  const PylithInt sOff_x[],
                                  const PylithScalar s[],
                                  const PylithScalar s_t[],
                                  const PylithScalar s_x[],
                                  const PylithInt aOff[],
                                  const PylithInt aOff_x[],
                                  const PylithScalar a[],
                                  const PylithScalar a_t[],
                                  const PylithScalar a_x[],
                                  const PylithReal t,
                                  const PylithReal utshift,
                                  const PylithScalar x[],
                                  const PylithInt numConstants,
                                  const PylithScalar constants[],
                                  PylithScalar Jg0[]) {
    const PylithInt _numS = 2;

    PylithInt i;

    assert(_numS == numS);

    for (i = 0; i < dim; ++i) {
        Jg0[i*dim+i] += 1.0;
    } // for
} // Jg0uv

// End of file
