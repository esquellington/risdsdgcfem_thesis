#ifndef S2_NS_ILSS_CGS_H
#define S2_NS_ILSS_CGS_H

#include "IIterativeLinearSystemSolver.h"
#include "real_array.h"

namespace S2 { namespace ns {

/* CG-Squared algorithm with a constraint projection matrix S
   Implementation from the paper:
     "HOW FAST ARE NONSYMMETRIC MATRIX ITERATIONS?*"
     SIAM J. MATRIX ANAL. APPL.
     Vol. 13, No. 3, pp. 778-795, July 1992

   \todo Stopping criteria taken from CG... seems reasonable, but I
   should check the original CGS paper...

   //\todo MOVE CODE TO CPP
*/
class ILSS_CGS: public IIterativeLinearSystemSolver
{
public:
    finline ILSS_CGS()
    : m_r(0)
    , m_d(0)
    , m_q(0)
    , m_tilde_r0(0)
    , m_u(0)
    , m_v(0)
    , m_w(0)
    , m_A_times_w(0)
    , m_p(0)
        {}
    ~ILSS_CGS()
        {
            if( m_r ) delete[] m_r;
            if( m_d ) delete[] m_d;
            if( m_q ) delete[] m_q;
            if( m_tilde_r0 ) delete[] m_tilde_r0;
            if( m_u ) delete[] m_u;
            if( m_v ) delete[] m_v;
            if( m_w ) delete[] m_w;
            if( m_A_times_w ) delete[] m_A_times_w;
            if( m_p ) delete[] m_p;
        }

    void Init( unsigned size )
        {
            IIterativeLinearSystemSolver::Init(size);
            m_r = new Real[size];
            m_d = new Real[size];
            m_q = new Real[size];
            m_tilde_r0 = new Real[size];
            m_u = new Real[size];
            m_v = new Real[size];
            m_w = new Real[size];
            m_A_times_w = new Real[size];
            m_p = new Real[size];
        }

    void Solve( D_Av d_Av, D_ApplyConstraints d_ApplyConstraints,
                Real* x, const Real* b,
                Real rel_epsilon, Real abs_epsilon, unsigned max_iterations )
        {
            NS_ASSERT( x != 0 && b != 0 );
            NS_ASSERT( !ns::real_array::is_nan( x, m_Size ) );
            NS_ASSERT( !ns::real_array::is_nan( b, m_Size ) );
            NS_ASSERT( max_iterations > 0 );
            NS_ASSERT( m_Size > 0 && m_r != 0 && m_d != 0 && m_q != 0
                       && m_tilde_r0 != 0 && m_u != 0 && m_v != 0
                       && m_w != 0 && m_A_times_w != 0 && m_p != 0 );

            Real abs_epsilon_sq( mal::Sq(abs_epsilon) );

            // r <== b - A*y //KNC r <== S * (b - A*x)
            d_Av( m_r, x ); // r = A*x
            ns::real_array::muleq_addeq( m_r, -1, b, m_Size ); // r = -1*r + b
            d_ApplyConstraints( m_r );
            NS_ASSERT( !ns::real_array::is_nan( m_r, m_Size ) );

            // \tilde r0 <== r
            ns::real_array::assign( m_tilde_r0, m_r, m_Size );

            // q = p = 0
            ns::real_array::zero( m_p, m_Size );
            ns::real_array::zero( m_q, m_Size );

            //\todo not sure if CGS-rho is equivalent to CG-delta...
            Real rho_0 = ns::real_array::sq_norm_2( m_r, m_Size ); //\todo the paper OTMCGMICS paper inits this differently to account for the KNC
            // \rho = 1
            Real rho_old( 1 );
            Real rho_new( rho_old );
            //
            /*\todo TESTING SAME INIT AS PURE CG... this crashes...
              Real rho_old( rho_0 );
              Real rho_new( rho_0 );
            */

            unsigned k = 0;
            Real rho_target( mal::Sq(rel_epsilon)*rho_0 ); //target rho magnitude
            /*\todo Since the convergence criteria is relative, it may ask
              for an arbitrarily small rho_new. If the initial y is
              "almost" exact, the error rho_0 is very small, and we
              require it to be reduced by a relative ls_prec factor, which
              may be impossible due to finite precision. Also, dot_d_q has
              been seen to be near 0 when rho_new is VERY small (10e-42)
              but still above rho_target (10e-48).

              \todo MAYBE the modified rho_0 computation in OTMCGMICS
              would solve this, at the expense of an additional
              multiplication by A...
            */
            //S2_DS_TRACE_LOCAL( "rho_0", rho_0 );
            //S2_DS_TRACE_LOCAL( "rho_target", rho_target );
            while( /*rho_new > ls_abs_prec_sq //absolute precision condition, not in original alg
                     && rho_new > rho_target
                     &&
                   */
                k < max_iterations )
            {
                rho_old = rho_new;
                rho_new = ns::real_array::dot( m_tilde_r0, m_r, m_Size );
                //S2_DS_TRACE_LOCAL( "rho_new", rho_new );
                Real beta = (k==0) ? 0 : rho_new / rho_old;
                //S2_DS_TRACE_LOCAL( "beta", beta );
                // u = r + beta*q
                ns::real_array::assign( m_u, m_r, m_Size );             //u = r
                ns::real_array::addeq_scaled( m_u, beta, m_q, m_Size ); //u += beta*q
                //S2_DS_TRACE_LOCAL_ARRAY( "u", m_u, m_Size );
                // p = u + beta(q + beta*p)
                ns::real_array::muleq( m_p, beta*beta, m_Size );      //p = beta*beta*p
                ns::real_array::addeq_scaled( m_p, beta, m_q, m_Size ); //p += beta*q
                ns::real_array::addeq( m_p, m_u, m_Size );              //p += u
                //S2_DS_TRACE_LOCAL_ARRAY( "p", p, m_Size );
                // v <== A*p //KNC v <== S*A*p
                d_Av( m_v, m_p );
                d_ApplyConstraints( m_v );
                NS_ASSERT( !ns::real_array::is_nan( m_v, m_Size ) );
                //S2_DS_TRACE_LOCAL_ARRAY( "v", m_v, m_Size );
                // \sigma = tilde_r0*v
                Real sigma = ns::real_array::dot( m_tilde_r0, m_v, m_Size );
                //S2_DS_TRACE_LOCAL( "sigma", sigma );
                if( mal::Abs(sigma) < abs_epsilon_sq ) break; //\todo This may happen if initial if we asked for too much precision...
                Real alpha = rho_new / sigma;
                //S2_DS_TRACE_LOCAL( "alpha", alpha );
                // q = u - alpha*v
                ns::real_array::assign( m_q, m_u, m_Size );              //q = u
                ns::real_array::subeq_scaled( m_q, alpha, m_v, m_Size ); //q -= alpha*v
                //S2_DS_TRACE_LOCAL_ARRAY( "q", m_q, m_Size );
                // w = u+q
                ns::real_array::assign( m_w, m_u, m_Size );
                ns::real_array::addeq( m_w, m_q, m_Size );
                // r = r - alpha*A*(u+q)
                {
                    d_Av( m_A_times_w, m_w );
                    d_ApplyConstraints( m_A_times_w );
                    NS_ASSERT( !ns::real_array::is_nan( m_A_times_w, m_Size ) );
                    ns::real_array::subeq_scaled( m_r, alpha, m_A_times_w, m_Size );
                    //S2_DS_TRACE_LOCAL_ARRAY( "r", m_r, m_Size );
                }
                // x = x + alpha*(u+q)
                d_ApplyConstraints( m_w ); //TEMP: added "by intuition"...
                ns::real_array::addeq_scaled( x, alpha, m_w, m_Size );
                k++;
            }
            NS_ASSERT( !ns::real_array::is_nan( x, m_Size ) );
            // results \todo THIS IS WRONG, rho_new CAN BE NEGATIVE, not proper prec... maybe Abs(rho_new)??
            m_NumIterations = k;
            m_RelResidual = (rho_0>0) ? mal::Sqrt( rho_new / rho_0 ) : 0;
            m_AbsResidual = mal::Sqrt( rho_new );
            /*
              S2_DS_TRACE_LOCAL( "ls_iter", ls_iter );
              S2_DS_TRACE_LOCAL( "ls_rel_prec", ls_rel_prec );
              S2_DS_TRACE_LOCAL( "ls_abs_prec", ls_abs_prec );
            */
        }

private:
    Real* m_r;
    Real* m_d;
    Real* m_q;
    Real* m_tilde_r0;
    Real* m_u;
    Real* m_v;
    Real* m_w;
    Real* m_A_times_w;
    Real* m_p;
};

}} //namespace S2::ns

#endif // S2_NS_ILSS_CGS_H
