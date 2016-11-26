#ifndef S2_NS_ILSS_CG_H
#define S2_NS_ILSS_CG_H

#include "IIterativeLinearSystemSolver.h"
#include "real_array.h"

namespace S2 { namespace ns {

/* Conjugate Gradients

   Method based on the paper "On the modified conjugate gradient
   method in cloth simulation", without Preconditioning.

   The method uses a projection S that selects dynamic nodes only,
   here implemented by the d_ApplyConstraints() method.

   \todo MOVE CODE TO CPP
*/
class ILSS_CG: public IIterativeLinearSystemSolver
{
public:
    finline ILSS_CG() : m_r(0), m_d(0), m_q(0) {}
    ~ILSS_CG()
        {
            if( m_r ) delete[] m_r;
            if( m_d ) delete[] m_d;
            if( m_q ) delete[] m_q;
        }

    void Init( unsigned size )
        {
            IIterativeLinearSystemSolver::Init(size);
            m_r = new Real[size];
            m_d = new Real[size];
            m_q = new Real[size];
        }

    void Solve( D_Av d_Av, D_ApplyConstraints d_ApplyConstraints,
                Real* x, const Real* b,
                Real rel_epsilon, Real abs_epsilon, unsigned max_iterations )
        {
            NS_ASSERT( x != 0 && b != 0 );
            NS_ASSERT( !ns::real_array::is_nan( x, m_Size ) );
            NS_ASSERT( !ns::real_array::is_nan( b, m_Size ) );
            NS_ASSERT( max_iterations > 0 );
            NS_ASSERT( m_Size > 0 && m_r != 0 && m_d != 0 && m_q != 0 );

            Real abs_epsilon_sq( mal::Sq(abs_epsilon) );

            //\todo 2003_OnTheModifiedConjugateGradientMethodInClothSimulation suggests
            // x = Sy + (I-S)z, where y = x0 (guess) and z = prescribed/constrained values

            // r <== b - A*x //KNC r <== S * (b - A*x)
            d_Av( m_r, x ); // r = A*x
            ns::real_array::muleq_addeq( m_r, -1, b, m_Size ); // r = -1*r + b
            d_ApplyConstraints( m_r ); // r = S*r
            NS_ASSERT( !ns::real_array::is_nan( m_r, m_Size ) );

            // d <== r
            ns::real_array::assign( m_d, m_r, m_Size );
            //TEMP: Unnecessary d_ApplyConstraints( d ); if no precod, otherwise d = S*P^-1*r

            // delta <== r^T * r
            Real delta_0 = ns::real_array::sq_norm_2( m_r, m_Size ); //\todo the paper OTMCGMICS paper inits this differently to account for the KNC
            Real delta_new( delta_0 );
            Real delta_old( delta_0 );

#ifdef __S2_NS_ENABLE_TRACE_ILSS_RESIDUAL
            NS_LOG("CG[%4d] delta = %f", -1, mal::Sqrt(delta_new) );
#endif

            unsigned int k( 0 );
            Real delta_target( mal::Sq(rel_epsilon)*delta_0 ); //target delta magnitude
            /*\todo Since the convergence criteria is relative, it may ask
              for an arbitrarily small delta_new. If the initial y is
              "almost" exact, the error delta_0 is very small, and we
              require it to be reduced by a relative ls_prec factor, which
              may be impossible due to finite precision. Also, dot_d_q has
              been seen to be near 0 when delta_new is VERY small (10e-42)
              but still above delta_target (10e-48).

              \todo MAYBE the modified delta_0 computation in OTMCGMICS
              would solve this, at the expense of an additional
              multiplication by A...
            */
            while( delta_new > abs_epsilon_sq //absolute precision condition, not in original alg
                   && delta_new > delta_target
                   && k < max_iterations )
            {
                // q <== A*d //KNC q <== S*A*d
                d_Av( m_q, m_d );
                d_ApplyConstraints( m_q );
                NS_ASSERT( !ns::real_array::is_nan( m_q, m_Size ) );

                // alpha <== delta_new/dot(d,q)
                Real dot_d_q( ns::real_array::dot( m_d, m_q, m_Size ) );
                if( mal::Abs(dot_d_q) < abs_epsilon_sq ) break; //\todo This may happen if initially we asked for too much precision...
                Real alpha( delta_new / dot_d_q );

                // x <== x + alpha * d
                ns::real_array::addeq_scaled( x, alpha, m_d, m_Size );
                NS_ASSERT( !ns::real_array::is_nan( x, m_Size ) );

                // r <== r - alpha * q
                ns::real_array::subeq_scaled( m_r, alpha, m_q, m_Size );
                NS_ASSERT( !ns::real_array::is_nan( m_r, m_Size ) );

                /* \todo Consider restart to avoid residual error accumulation (see CG paper)
                   if( 0==k%50 )
                     r = b - A*x;
                   else
                     r = r - alpha * q
                */

                // delta_new <== r^T * r;
                delta_old = delta_new;
                delta_new = ns::real_array::sq_norm_2( m_r, m_Size );

                /*\todo HERE we could early-out on convergence, as any
                  further computations are ONLY used if an additional
                  iteration is necessary
                */

                /* \todo Consider restarting to avoid cancellation error (see CG paper)
                   if( delta_new <= delta_target )
                   {
                     r = b - A*x;
                     delta_new = r.SqNorm2();
                   }
                */

                // beta <== delta_new / delta_old
                Real beta( delta_new / delta_old );

                // d <== r + beta*d //KNC d <== S * (r + beta*d)
                ns::real_array::muleq_addeq( m_d, beta, m_r, m_Size );
                d_ApplyConstraints( m_d );
                NS_ASSERT( !ns::real_array::is_nan( m_d, m_Size ) );

#ifdef __S2_NS_ENABLE_TRACE_ILSS_RESIDUAL
                NS_LOG("CG[%4u] delta = %f", k, mal::Sqrt(delta_new) );
#endif

                k++;
            }
            NS_ASSERT( !ns::real_array::is_nan( x, m_Size ) );
            // Results
            m_NumIterations = k;
            m_RelResidual = (delta_0>0) ? mal::Sqrt( delta_new / delta_0 ) : 0;
            m_AbsResidual = mal::Sqrt( delta_new );
        }

private:
    Real* m_r;
    Real* m_d;
    Real* m_q;
};

}} //namespace S2::ns

#endif // S2_NS_ILSS_CG_H
