#ifndef S2_NS_ILSS_SYMMQMR_H
#define S2_NS_ILSS_SYMMQMR_H

#include "IIterativeLinearSystemSolver.h"
#include "real_array.h"

namespace S2 { namespace ns {

/* Symmetric-QMR

   Method based on the paper "A New Krylov Subspace Method For
   Symmetric Indefinite Linear Systems", Algorithm 1, pg 5.

   \note No preconditioning supported, see original algorithm if
   required. Same alg is implemented in PhysBAM with different
   preconditioning strategy. According to original paper above, with
   no preconditioning SYMMQMR is equivalent to MINRES.

   \todo Stopping criteria adapted from CG... may be WRONG

   \todo MOVE CODE TO CPP
*/
class ILSS_SYMMQMR: public IIterativeLinearSystemSolver
{
public:
    finline ILSS_SYMMQMR() : m_r(0), m_t(0), m_q(0), m_d(0) {}
    ~ILSS_SYMMQMR()
        {
            if( m_r ) delete[] m_r;
            if( m_t ) delete[] m_t;
            if( m_q ) delete[] m_q;
            if( m_d ) delete[] m_d;
        }

    void Init( unsigned size )
        {
            IIterativeLinearSystemSolver::Init(size);
            m_r = new Real[size];
            m_t = new Real[size];
            m_q = new Real[size];
            m_d = new Real[size];
        }

    void Solve( D_Av d_Av, D_ApplyConstraints d_ApplyConstraints,
                Real* x, const Real* b,
                Real rel_epsilon, Real abs_epsilon, unsigned max_iterations )
        {
            NS_ASSERT( x != 0 && b != 0 );
            NS_ASSERT( !ns::real_array::is_nan( x, m_Size ) );
            NS_ASSERT( !ns::real_array::is_nan( b, m_Size ) );
            NS_ASSERT( max_iterations > 0 );
            NS_ASSERT( m_Size > 0 && m_r != 0 && m_t != 0 && m_q != 0 && m_d != 0 );

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

            // t <== r ( if precond t = M_1^-1 r)
            ns::real_array::assign( m_t, m_r, m_Size );

            // q << t ( if precond q = M_2^-1 t)
            ns::real_array::assign( m_q, m_t, m_Size );

            Real tau_old( ns::real_array::norm_2( m_t, m_Size ) );
            Real nu_old( 0 );
            Real rho_old( ns::real_array::dot( m_r, m_q, m_Size ) );

            Real delta_0( ns::real_array::sq_norm_2( m_r, m_Size ) ); //\todo the paper OTMSYMMQMRMICS paper inits this differently to account for the KNC
            Real delta_new( delta_0 );
            Real delta_old( delta_0 );

#ifdef __S2_NS_ENABLE_TRACE_ILSS_RESIDUAL
            NS_LOG("SYMMQMR[%4d] delta = %f", -1, mal::Sqrt(delta_new) );
#endif

            unsigned int k( 0 );
            Real delta_target( mal::Sq(rel_epsilon)*delta_0 ); //target delta magnitude
            while( delta_new > abs_epsilon_sq //absolute precision condition, not in original alg
                   && delta_new > delta_target
                   && k < max_iterations )
            {
                // t <== A*q \todo t overwritten at k==0 => NO NEED TO INIT IT??
                d_Av( m_t, m_q );
                d_ApplyConstraints( m_t );
                NS_ASSERT( !ns::real_array::is_nan( m_t, m_Size ) );

                // theta = q^T * t
                Real theta( ns::real_array::dot( m_q, m_t, m_Size ) );
                if( mal::Abs(theta) < abs_epsilon_sq ) break; //\todo Can this happen??

                // alpha = rho_old / theta
                Real alpha( rho_old / theta );

                // r <== r - alpha * t
                ns::real_array::subeq_scaled( m_r, alpha, m_t, m_Size );
                NS_ASSERT( !ns::real_array::is_nan( m_r, m_Size ) );

                // t <== r (\todo if precond t = M_1^-1 r)
                //ns::real_array::assign( m_t, m_r, m_Size );

                // delta_new <== r^T * r;
                delta_old = delta_new;
                delta_new = ns::real_array::sq_norm_2( m_r, m_Size );
                Real nu_new( mal::Sqrt(delta_new) / tau_old );
                Real c( 1.0/mal::Sqrt(1+mal::Sq(nu_new)) );
                Real tau_new( tau_old*nu_new*c );

                // d <== c^2 nu^2 d + c^2 alpha q;
                ns::real_array::muleq( m_d, mal::Sq(c*nu_old), m_Size );
                ns::real_array::addeq_scaled( m_d, mal::Sq(c)*alpha, m_q, m_Size );

                // x <== x + d
                ns::real_array::addeq( x, m_d, m_Size );
                NS_ASSERT( !ns::real_array::is_nan( x, m_Size ) );

                // t <== r ( if precond t = M_1^-1 r)
                //ns::real_array::assign( m_t, m_r, m_Size );

                // rho = r^T * t ( == delta_new if no precond)
                // Real rho_new( ns::real_array::dot( m_r, m_r, m_Size ) );
                Real rho_new( delta_new );

                // beta <== rho_new/rho_old
                if( mal::Abs(rho_old) < abs_epsilon_sq ) break; //\todo This may happen if initially we asked for too much precision...

                Real beta( rho_new / rho_old );

                // q <== u + beta*q //KNC d <== S * (r + beta*p) \todo KNC are ok here???
                // where u = t = r if no preconditioning
                ns::real_array::muleq_addeq( m_q, beta, m_r, m_Size );
                d_ApplyConstraints( m_q );
                NS_ASSERT( !ns::real_array::is_nan( m_q, m_Size ) );

#ifdef __S2_NS_ENABLE_TRACE_ILSS_RESIDUAL
                NS_LOG("SYMMQMR[%4u] delta = %f", k, mal::Sqrt(delta_new) );
#endif
                nu_old = nu_new;
                tau_old = tau_new;
                rho_old = rho_new;

                k++;
            }
            NS_ASSERT( !ns::real_array::is_nan( x, m_Size ) );
            // Results
            m_NumIterations = k;
            m_RelResidual = (delta_0>0) ? mal::Sqrt( delta_new / delta_0 ) : 0;
            m_AbsResidual = mal::Sqrt( delta_new );
            // NS_LOG_WARNING("SYMMQMR ends with %u iter, abs_residual = %f", m_NumIterations, m_AbsResidual );
        }

private:
    Real* m_r;
    Real* m_t;
    Real* m_q;
    Real* m_d;
};

}} //namespace S2::ns

#endif // S2_NS_ILSS_SYMMQMR_H
