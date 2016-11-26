#ifndef S2_NS_ILSS_MINRES_H
#define S2_NS_ILSS_MINRES_H

#include "IIterativeLinearSystemSolver.h"
#include "real_array.h"

namespace S2 { namespace ns {

/* MINRES

   Method based on the paper "Solution of Sparse Indefinite Systems of
   Linear Equations" C.C.Paige and M.A.Saunders, 1975

   Solves Ax = b for symmetric A (definite/indefinite), with
   "something-like" least-squares solution if A is singular.

   \note Code adapted from
   Eigen... http://eigen.tuxfamily.org/dox/unsupported/MINRES_8h_source.html

   \todo MOVE CODE TO CPP
*/
class ILSS_MINRES: public IIterativeLinearSystemSolver
{
public:
    finline ILSS_MINRES() : m_v(0), m_v0(0), m_v1(0), m_w(0), m_w1(0), m_p(0), m_p0(0), m_p0_tmp(0) {}
    ~ILSS_MINRES()
        {
            if( m_v ) delete[] m_v;
            if( m_v0 ) delete[] m_v0;
            if( m_v1 ) delete[] m_v1;
            if( m_w ) delete[] m_w;
            if( m_w1 ) delete[] m_w1;
            if( m_p ) delete[] m_p;
            if( m_p0 ) delete[] m_p0;
            if( m_p0_tmp ) delete[] m_p0_tmp;
        }

    void Init( unsigned size )
        {
            IIterativeLinearSystemSolver::Init(size);
            m_v = new Real[size];
            m_v0 = new Real[size];
            m_v1 = new Real[size];
            m_w = new Real[size];
            m_w1 = new Real[size];
            m_p = new Real[size];
            m_p0 = new Real[size];
            m_p0_tmp = new Real[size];
        }

    void Solve( D_Av d_Av, D_ApplyConstraints d_ApplyConstraints,
                Real* x, const Real* b,
                Real rel_epsilon, Real abs_epsilon, unsigned max_iterations )
        {
            NS_ASSERT( x != 0 && b != 0 );
            NS_ASSERT( !ns::real_array::is_nan( x, m_Size ) );
            NS_ASSERT( !ns::real_array::is_nan( b, m_Size ) );
            NS_ASSERT( max_iterations > 0 );
            NS_ASSERT( m_Size > 0
                       && m_v != 0 && m_v0 != 0 && m_v1 != 0
                       && m_w != 0 && m_w1 != 0
                       && m_p != 0 && m_p0 != 0 && m_p0_tmp != 0 );

            Real abs_epsilon_sq( mal::Sq(abs_epsilon) );

            // v <== 0
            ns::real_array::zero( m_v, m_Size );
            // r <== b - A*x //KNC r <== S * (b - A*x)
            d_Av( m_v1, x ); // r = A*x
            ns::real_array::muleq_addeq( m_v1, -1, b, m_Size ); // r = -1*r + b
            d_ApplyConstraints( m_v1 ); // r = S*r
            NS_ASSERT( !ns::real_array::is_nan( m_v1, m_Size ) );

            Real norm_r_sq( ns::real_array::sq_norm_2( m_v1, m_Size ) );

            // w1 <== Precond * v1 (== v1 if no preconditioner)
            ns::real_array::assign( m_w1, m_v1, m_Size );
            Real beta_new_sq( ns::real_array::dot( m_v1, m_w1, m_Size ) );
            NS_ASSERT( beta_new_sq >= 0 );
            Real beta_new( mal::Sqrt(beta_new_sq) );
            Real beta_one( beta_new );

            // v1 <== v1/beta_new
            ns::real_array::muleq( m_v1, 1.0/beta_new, m_Size );
            // w1 <== w1/beta_new
            ns::real_array::muleq( m_w1, 1.0/beta_new, m_Size );

            // Initialize other variables
            Real c(1.0); // the cosine of the Givens rotation
            Real c_old(1.0);
            Real s(0.0); // the sine of the Givens rotation
            Real s_old(0.0); // the sine of the Givens rotation
            ns::real_array::zero( m_p0, m_Size );
            ns::real_array::zero( m_p, m_Size );
            Real eta(1.0);

#ifdef __S2_NS_ENABLE_TRACE_ILSS_RESIDUAL
            NS_LOG_WARNING("MINRES[%d] delta = %f", -1, mal::Sqrt(norm_r_sq) );
#endif

            unsigned int k( 0 );
            Real norm_b_sq( ns::real_array::sq_norm_2( b, m_Size ) );
            Real delta_target( mal::Sq(rel_epsilon)*norm_b_sq ); //target delta magnitude
            while( norm_r_sq > abs_epsilon_sq //absolute precision condition, not in original alg
                   && norm_r_sq > delta_target
                   && k < max_iterations )
            {
                Real beta(beta_new);
                // v0 <== v
                ns::real_array::assign( m_v0, m_v, m_Size );
                // v <== v1
                ns::real_array::assign( m_v, m_v1, m_Size );
                // w <== w1
                ns::real_array::assign( m_w, m_w1, m_Size );
                //\todo CONSIDER swaps instead of assigns for some vectors (v0,v,v1,w,w1) to avoid redundant copies...
                // v1 <== A*w - beta*v0
                d_Av( m_v1, m_w );
                d_ApplyConstraints( m_v1 );
                NS_ASSERT( !ns::real_array::is_nan( m_v1, m_Size ) );
                ns::real_array::addeq_scaled( m_v1, -beta, m_v0, m_Size );
                Real alpha( ns::real_array::dot( m_v1, m_w, m_Size ) );
                // v1 <== v1 - alpha*v
                ns::real_array::addeq_scaled( m_v1, -alpha, m_v, m_Size );
                // w1 <== Precond * v1 (== v1 if no preconditioner)
                ns::real_array::assign( m_w1, m_v1, m_Size );
                beta_new_sq = ns::real_array::dot( m_v1, m_w1, m_Size ); //== |v1|^2 as v1==w2 if no preconditioner
                NS_ASSERT( beta_new_sq >= 0 );
                beta_new = mal::Sqrt(beta_new_sq);
                // v1 <== v1/beta_new
                ns::real_array::muleq( m_v1, 1.0/beta_new, m_Size );
                // w1 <== w1/beta_new
                ns::real_array::muleq( m_w1, 1.0/beta_new, m_Size );

                // Givens rotation
                const Real r2( s*alpha+c*c_old*beta ); // s, s_old, c and c_old are still from previous iteration
                const Real r3( s_old*beta ); // s, s_old, c and c_old are still from previous iteration
                const Real r1_hat( c*alpha-c_old*s*beta );
                const Real r1( mal::Sqrt( mal::Sq(r1_hat) + beta_new_sq ) );
                c_old = c; // store for next iteration
                s_old = s; // store for next iteration
                c = r1_hat/r1; // new cosine
                s = beta_new/r1; // new sine

                // Update solution
                // p0_tmp <== p0, p0 <== p
                ns::real_array::assign( m_p0_tmp, m_p0, m_Size );
                ns::real_array::assign( m_p0, m_p, m_Size ); //\todo See if temporal p0_tmp could reuse some existing vector instead that will be written-to in the next iter
                // p <== (w - r2*p - r3*p0_tmp) / r1
                ns::real_array::muleq_addeq( m_p, -r2, m_w, m_Size );
                ns::real_array::addeq_scaled( m_p, -r3, m_p0_tmp, m_Size );
                ns::real_array::muleq( m_p, 1.0/r1, m_Size ); //\todo Consider more complex ops to avoid 3 passes...
                // x <== x + beta_one*c*eta*p
                ns::real_array::addeq_scaled( x, beta_one*c*eta, m_p, m_Size );
                norm_r_sq *= mal::Sq(s);
                eta = -s*eta;

#ifdef __S2_NS_ENABLE_TRACE_ILSS_RESIDUAL
                NS_LOG_WARNING("MINRES[%u] delta = %f", k, mal::Sqrt(norm_r_sq) );
#endif
                k++;
            }
            NS_ASSERT( !ns::real_array::is_nan( x, m_Size ) );
            // Results
            m_NumIterations = k;
            m_RelResidual = mal::Sqrt( norm_r_sq / norm_b_sq );
            m_AbsResidual = mal::Sqrt( norm_r_sq );
            // NS_LOG_WARNING("MINRES ends with %u iter, abs_residual = %f", m_NumIterations, m_AbsResidual );
        }
private:
    Real* m_v;
    Real* m_v0;
    Real* m_v1;
    Real* m_w;
    Real* m_w1;
    Real* m_p;
    Real* m_p0;
    Real* m_p0_tmp;
};

}} //namespace S2::ns

#endif // S2_NS_ILSS_MINRES_H
