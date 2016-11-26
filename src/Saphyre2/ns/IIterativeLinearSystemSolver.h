#ifndef S2_NS_I_ITERATIVE_LINEAR_SYSTEM_SOLVER_H
#define S2_NS_I_ITERATIVE_LINEAR_SYSTEM_SOLVER_H

#include "Config.h"
#include "MatrixD.h"
#include <boost/function.hpp>

// #define __S2_NS_ENABLE_TRACE_ILSS_RESIDUAL

namespace S2 { namespace ns {

/* Iterative LS solver interface
   - Init() allocates all required temporary memory
   - Solve() solves the LS with specific method

   \todo Most Iterative solvers can be "restarted" to prevent
   numerical drift or, in the case of GMRES, to reduce memory
   cost. PhysBAM has this concept integrated into all iterative
   solvers, consider promoting it from GMRES subclass to the common
   interface.
*/
class IIterativeLinearSystemSolver
{
public:

    typedef boost::function<void (Real*,const Real*)> D_Av; //arg0 = A * arg1
    typedef boost::function<void (Real*)> D_ApplyConstraints;

public:
    inline IIterativeLinearSystemSolver() : m_Size(0), m_NumIterations(0), m_RelResidual(0), m_AbsResidual(0) {}
    virtual ~IIterativeLinearSystemSolver() {}

    virtual void Init( unsigned size ) { NS_ASSERT( size > 0 ); m_Size = size; }

    // Solve Ax=b for x, with initial approximation x = x_0
    virtual void Solve( D_Av d_Av, D_ApplyConstraints d_ApplyConstraints, Real* x, const Real* b, Real rel_epsilon, Real abs_epsilon, unsigned max_iterations ) = 0;

    inline unsigned GetSize() const { return m_Size; }
    inline unsigned GetNumIterations() const { return m_NumIterations; }
    inline Real GetRelResidual() const { return m_RelResidual; }
    inline Real GetAbsResidual() const { return m_AbsResidual; }

protected:
    unsigned m_Size;
    unsigned m_NumIterations;
    Real m_RelResidual;
    Real m_AbsResidual;
};

}} //namespace S2::ns


//----------------------------------------------------------------
// Helper Methods \todo SHOULD BE IMPLEMENTED IN A CPP NOT INLINE!
//----------------------------------------------------------------

#include "real_array.h"
#include <Mal/GRandom.h>

namespace S2 { namespace ns {

/* Assemble matrix A from Av products.

   \todo THIS IS EXTREMELY inefficient, better find exact jacobian
   using single-component formulas or, at least, use the fact that
   d_Av is called with a SINGLE non-zero component... consider
   d_Av_component( tmp1, i).  Also consider block-components (eg:
   per-node blocks => noderows)
*/
inline void Assemble_A_From_Av( IIterativeLinearSystemSolver::D_Av d_Av, unsigned size,
                                Real* tmp0, Real* tmp1,
                                ns::MatrixD& A )
{
    for( unsigned int j=0; j<size; j++ )
    {
        // Init e_j unit vector
        ns::real_array::zero( tmp0, size );
        tmp0[j] = Real(1);
        // Compute A_j column
        d_Av( tmp1, tmp0 );
        // Assemble A_j column
        for( unsigned int i=0; i<size; i++ ) A(i,j) = tmp1[i];
    }
}

// Compute Largest (Real) eigenvalue,eigenvector pair using Av products
inline Real ComputeLargestEigenvalue_Symmetric_Shifted( IIterativeLinearSystemSolver::D_Av d_Av, unsigned size,
                                                        Real* tmp0, Real* tmp1, Real* tmp2,
                                                        unsigned int max_iter, Real epsilon, Real shift_eigen_value,
                                                        Real& eigen_value, Real* eigen_vector = 0 )
{
    Real *v0( tmp1 );
    Real e0(0);
    Real *v1( tmp2 );
    do
    {
        for( unsigned int i=0; i<size; i++ ) v1[i] = mal::RandomF<Real>(-1,1);
    } while( real_array::norm_2( v1, size ) < epsilon );
    Real e1( real_array::norm_2( v1, size ) );
    MAL_ASSERT( e1 > epsilon );
    unsigned int iter(0);
    do
    {
        e0 = e1;
        real_array::assign( v0, v1, size );
        real_array::muleq( v0, mal::Rcp(e1), size );
        d_Av( v1, v0 ); //v1 = A*v0
        real_array::subeq_scaled( v1, shift_eigen_value, v0, size ); //v1' = v1 - shift*Id*v0
        //\todo SHOULD APPLY KNC, which, I think, MUST leave kinematic entries untouched...
        //LS_ZeroKNC( m_CG_vec_q ); //\todo this does NOT SEEM RIGHT... we do this inside CG to avoid changing KNC-node magnitudes, not to turn them into 0
        e1 = real_array::norm_2( v1, size );
        iter++;
    }
    while( iter < max_iter //Limit cost
           && e1 > epsilon //Ensure non-collapsing eigenvector
           && mal::Abs(e1-e0) > epsilon ); //Ensure significant convergence
    // Compute final eigen_vector and its sign
    if( e1 > epsilon )
    {
        if( 0 != eigen_vector )
        {
            real_array::assign( eigen_vector, v1, size );
            real_array::muleq( eigen_vector, mal::Rcp(e1), size );
        }
        // Compute eigen_value sign by checking if v0*v1 < 0 (which means it has inverted when multiplied by M in v1 = M * v0)
        eigen_value = ( real_array::dot(v0,v1,size) >= 0 ) ? e1 : -e1;
    }
    // No non-zero component found in iterated eigen_vector, no valid eigen_value
    eigen_value = 0;
    if( 0 != eigen_vector ) real_array::zero( eigen_vector, size );
    return 0;
}

inline Real ComputeMinMaxEigenvalueRatio_Symmetric( IIterativeLinearSystemSolver::D_Av d_Av, unsigned size,
                                                    Real* tmp0, Real* tmp1, Real* tmp2,
                                                    unsigned int max_iter, Real epsilon )
{
    // Find largest EV of A using power method
    Real lambda_max(0);
    ComputeLargestEigenvalue_Symmetric_Shifted( d_Av, size,
                                                tmp0, tmp1, tmp2,
                                                max_iter, epsilon, 0,
                                                lambda_max );
    // Find smallest EV of A using shifted power method (http://math.stackexchange.com/questions/271864/power-iteration-smallest-eigenvalue)
    Real lambda_min(0);
    ComputeLargestEigenvalue_Symmetric_Shifted( d_Av, size,
                                                tmp0, tmp1, tmp2,
                                                max_iter, epsilon, lambda_max,
                                                lambda_min );
    lambda_min += lambda_max;
    return lambda_min/lambda_max;
}

}} //namespace S2::ns

#endif // S2_NS_I_ITERATIVE_LINEAR_SYSTEM_SOLVER_H
