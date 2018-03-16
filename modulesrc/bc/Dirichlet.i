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

/** @file modulesrc/bc/DirichletBoundary.i
 *
 * @brief Python interface to C++ DirichletBoundary object.
 */

namespace pylith {
    namespace bc {

        class Dirichlet :
	    public pylith::bc::BoundaryCondition,
	    public pylith::feassemble::ConstraintPointwise {
	    
	    // PUBLIC METHODS /////////////////////////////////////////////////
	public:
	    
	    /// Default constructor.
	    Dirichlet(void);
	    
	    /// Destructor.
	    ~Dirichlet(void);
	    
	    /// Deallocate PETSc and local data structures.
	    virtual
	    void deallocate(void);
	    
	    /** Verify configuration is acceptable.
	     *
	     * @param[in] solution Solution field.
	     */
	    void verifyConfiguration(const pylith::topology::Field& solution) const;
	    
	    /** Initialize boundary condition.
	     *
	     * @param[in] solution Solution field.
	     */
	    void initialize(const pylith::topology::Field& solution);
	    
	    /** Set constrained values in solution field.
	     *
	     * @param[out] solution Solution field.
	     * @param[in] t Current time.
	     */
	    void setSolution(pylith::topology::Field* solution,
			     const double t);
	    
	    // PROTECTED METHODS //////////////////////////////////////////////////
	protected:
	    
	    /** Setup auxiliary subfields (discretization and query fns).
	     *
	     * Create subfields in auxiliary fields (includes name of the field,
	     * vector field type, discretization, and scale for
	     * nondimensionalization) and set query functions for filling them
	     * from a spatial database.
	     *
	     * @attention The order of the calls to subfieldAdd() must match the
	     * order of the auxiliary fields in the FE kernels.
	     */
	    virtual
	    void _auxFieldsSetup(void) = 0;
	    
	    /** Set kernels for RHS residual G(t,s).
	     *
	     * Potentially, there are g0 and g1 kernels for each equation. If no
	     * kernel is needed, then set the kernel function to NULL.
	     *
	     * @param solution Solution field.
	     */
	    virtual
	    void _setFEKernelsConstraint(const topology::Field& solution) = 0;
	    
        }; // class Dirichlet
	
    } // bc
} // pylith


// End of file
