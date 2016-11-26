#ifndef S2_NS_ILSS_CR_H
#define S2_NS_ILSS_CR_H

#include "IIterativeLinearSystemSolver.h"
#include "real_array.h"

namespace S2 { namespace ns {

/* Conjugate Residuals

   Method based on the thesis "Iterative Methods For Singular Linear
   Equations and Least Squares Problems", Table 2.12, pg 27.

   \todo FAILS MISERABLY even for unconstrained, SPD, systems that CG
   solves without trouble... which SUCKS HARD... implementation seems
   equivalent to source papers and to PhysBAM (independently
   written!)...

   \todo Stopping criteria adapted from CG... may be WRONG

   \todo MOVE CODE TO CPP
*/
class ILSS_CR: public IIterativeLinearSystemSolver
{
public:
    finline ILSS_CR() : m_r(0), m_p(0), m_z(0), m_w(0) {}
    ~ILSS_CR()
        {
            if( m_r ) delete[] m_r;
            if( m_p ) delete[] m_p;
            if( m_z ) delete[] m_z;
            if( m_w ) delete[] m_w;
        }

    void Init( unsigned size )
        {
            IIterativeLinearSystemSolver::Init(size);
            m_r = new Real[size];
            m_p = new Real[size];
            m_z = new Real[size];
            m_w = new Real[size];
        }

    void Solve( D_Av d_Av, D_ApplyConstraints d_ApplyConstraints,
                Real* x, const Real* b,
                Real rel_epsilon, Real abs_epsilon, unsigned max_iterations )
        {
            NS_ASSERT( x != 0 && b != 0 );
            NS_ASSERT( !ns::real_array::is_nan( x, m_Size ) );
            NS_ASSERT( !ns::real_array::is_nan( b, m_Size ) );
            NS_ASSERT( max_iterations > 0 );
            NS_ASSERT( m_Size > 0 && m_r != 0 && m_p != 0 && m_z != 0 && m_w != 0 );

            Real abs_epsilon_sq( mal::Sq(abs_epsilon) );

            /*TEMP: PhysBAM applies boundary conditions to X here, we
              assume it's already been done before calling Solve()
              d_ApplyConstraints( x );
            */

            // r <== b - A*x //KNC r <== S * (b - A*x)
            d_Av( m_r, x ); // r = A*x
            ns::real_array::muleq_addeq( m_r, -1, b, m_Size ); // r = -1*r + b
            d_ApplyConstraints( m_r ); // r = S*r
            NS_ASSERT( !ns::real_array::is_nan( m_r, m_Size ) );

            // p <== r
            ns::real_array::assign( m_p, m_r, m_Size );

            // z <== Ar
            d_Av( m_z, m_r );
            d_ApplyConstraints( m_z ); // z = S*z
            NS_ASSERT( !ns::real_array::is_nan( m_z, m_Size ) );

            // w <== Ap == Ar == z
            ns::real_array::assign( m_w, m_z, m_Size );

            // delta <== r^T * r
            Real delta_0( ns::real_array::sq_norm_2( m_r, m_Size ) ); //\todo the paper OTMCRMICS paper inits this differently to account for the KNC
            Real delta_new( delta_0 );
            Real delta_old( delta_0 );

            // mu = r^T * z
            Real mu_new( ns::real_array::dot( m_r, m_z, m_Size ) );
            Real mu_old( mu_new );

            unsigned int k( 0 );
            Real delta_target( mal::Sq(rel_epsilon)*delta_0 ); //target delta magnitude
            /*\todo Since the convergence criteria is relative, it may ask
              for an arbitrarily small delta_new. If the initial y is
              "almost" exact, the error delta_0 is very small, and we
              require it to be reduced by a relative ls_prec factor, which
              may be impossible due to finite precision. Also, dot_d_q has
              been seen to be near 0 when delta_new is VERY small (10e-42)
              but still above delta_target (10e-48).

              \todo MAYBE the modified delta_0 computation in OTMCRMICS
              would solve this, at the expense of an additional
              multiplication by A...
            */
            while( delta_new > abs_epsilon_sq //absolute precision condition, not in original alg
                   && delta_new > delta_target
                   && k < max_iterations )
            {
                // norm_sq_w
                Real norm_sq_w( ns::real_array::sq_norm_2( m_w, m_Size ) );
                if( norm_sq_w < abs_epsilon_sq ) break; //\todo Can this happen??

                // alpha = mu / norm_sq_w
                Real alpha( mu_old / norm_sq_w );

                // x <== x + alpha * p
                ns::real_array::addeq_scaled( x, alpha, m_p, m_Size );
                NS_ASSERT( !ns::real_array::is_nan( x, m_Size ) );

                // \todo Consider restart to avoid residual error accumulation (see CG paper)
                if( false )//restart )
                {
                    // r <== b - A*x //KNC r <== S * (b - A*x)
                    d_Av( m_r, x ); // r = A*x
                    ns::real_array::muleq_addeq( m_r, -1, b, m_Size ); // r = -1*r + b
                    d_ApplyConstraints( m_r ); // r = S*r
                    NS_ASSERT( !ns::real_array::is_nan( m_r, m_Size ) );
                }
                else
                {
                    // r <== r - alpha * w
                    ns::real_array::subeq_scaled( m_r, alpha, m_w, m_Size );
                    NS_ASSERT( !ns::real_array::is_nan( m_r, m_Size ) );
                }

                // delta_new <== r^T * r;
                delta_old = delta_new;
                delta_new = ns::real_array::sq_norm_2( m_r, m_Size );

                // z <== A*r
                d_Av( m_z, m_r );
                d_ApplyConstraints( m_z );
                NS_ASSERT( !ns::real_array::is_nan( m_z, m_Size ) );

                // mu = r^T * z
                mu_old = mu_new;
                mu_new = ns::real_array::dot( m_r, m_z, m_Size );

                // beta <== mu_new/mu_old
                if( mal::Abs(mu_old) < abs_epsilon_sq ) break; //\todo This may happen if initial if we asked for too much precision...
                Real beta( mu_new / mu_old );

                // p <== r + beta*p //KNC d <== S * (r + beta*p) \todo KNC are ok here???
                ns::real_array::muleq_addeq( m_p, beta, m_r, m_Size );
                d_ApplyConstraints( m_p );
                NS_ASSERT( !ns::real_array::is_nan( m_p, m_Size ) );

                // w <== z + beta*w //KNC w <== S * (z + beta*w) \todo KNC are ok here???
                ns::real_array::muleq_addeq( m_w, beta, m_z, m_Size );
                d_ApplyConstraints( m_w );
                NS_ASSERT( !ns::real_array::is_nan( m_w, m_Size ) );

                //DEBUG
                {
                    NS_LOG_WARNING("CR[%u] mu = %f, delta = %f", k, mu_new, delta_new );
                }

                k++;
            }
            NS_ASSERT( !ns::real_array::is_nan( x, m_Size ) );
            // Results
            m_NumIterations = k;
            m_RelResidual = (delta_0>0) ? mal::Sqrt( delta_new / delta_0 ) : 0;
            m_AbsResidual = mal::Sqrt( delta_new );
            // NS_LOG_WARNING("CR ends with %u iter, abs_residual = %f", m_NumIterations, m_AbsResidual );
        }

private:
    Real* m_r;
    Real* m_p;
    Real* m_z;
    Real* m_w;
};

}} //namespace S2::ns

#endif // S2_NS_ILSS_CR_H
