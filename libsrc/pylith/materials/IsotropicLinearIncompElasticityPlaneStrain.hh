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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/IsotropicLinearIncompElasticityPlaneStrain.hh
 *
 * @brief C++ class for isotropic linear incompressible elastic plane strain material.
 */

#if !defined(pylith_materials_isotropiclinearincompelasticityplanestrain_hh)
#define pylith_materials_isotropiclinearincompelasticityplanestrain_hh

// Include directives ---------------------------------------------------
#include "materialsfwd.hh" // forward declarations

#include "pylith/materials/Material.hh" // ISA Material

// Material -------------------------------------------------------------
/** @brief C++ class for isotropic linear incompressible elastic plane strain material.
 */

class pylith::materials::IsotropicLinearIncompElasticityPlaneStrain : public pylith::materials::Material { // class IsotropicLinearIncompElasticityPlaneStrain
    friend class TestIsotropicLinearIncompElasticityPlaneStrain; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    IsotropicLinearIncompElasticityPlaneStrain(void);

    /// Destructor.
    ~IsotropicLinearIncompElasticityPlaneStrain(void);

    /** Include inertia?
     *
     * @param[in] value Flag indicating to include inertial term.
     */
    void useInertia(const bool value);

    /** Include inertia?
     *
     * @returns True if including inertial term, false otherwise.
     */
    bool useInertia(void) const;

    /** Include body force?
     *
     * @param[in] value Flag indicating to include body force term.
     */
    void useBodyForce(const bool value);

    /** Include body force?
     *
     * @returns True if including body force term, false otherwise.
     */
    bool useBodyForce(void) const;

    /** Include reference stress/strain?
     *
     * @param value Flag indicating to include reference stress/strain.
     */
    void useReferenceState(const bool value);

    /** Use reference stress and strain in computation of stress and
     * strain?
     *
     * @returns True if using reference stress and strain, false otherwise.
     */
    bool useReferenceState(void) const;

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /// Setup auxiliary subfields (discretization and query fns).
    void _auxFieldSetup(void);

    /** Set kernels for RHS residual G(t,u).
     *
     * @param[in] solution Solution field.
     */
    void _setFEKernelsRHSResidual(const topology::Field& solution) const;

    /** Set kernels for RHS Jacobian G(t,u).
     *
     * @param[in] solution Solution field.
     */
    void _setFEKernelsRHSJacobian(const topology::Field& solution) const;

    /** Set kernels for LHS residual F(t,u,\dot{u}).
     *
     * @param[in] solution Solution field.
     */
    void _setFEKernelsLHSResidual(const topology::Field& solution) const;

    /** Set kernels for LHS Jacobian F(t,u,\dot{u}).
     *
     * @param[in] solution Solution field.
     */
    void _setFEKernelsLHSJacobian(const topology::Field& solution) const;


    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    bool _useInertia; ///< Flag to include inertial term.
    bool _useBodyForce; ///< Flag to include body force term.
    bool _useReferenceState; ///< Flag to use reference stress and strain.

    static const char* _pyreComponent; ///< Name of Pyre component.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    IsotropicLinearIncompElasticityPlaneStrain(const IsotropicLinearIncompElasticityPlaneStrain&); ///< Not implemented.
    const IsotropicLinearIncompElasticityPlaneStrain& operator=(const IsotropicLinearIncompElasticityPlaneStrain&); ///< Not implemented

}; // class IsotropicLinearIncompElasticityPlaneStrain

#endif // pylith_materials_isotropiclinearincompelasticityplanestrain_hh


// End of file
