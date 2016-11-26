#ifndef S2_NS_ILSS_RI_H
#define S2_NS_ILSS_RI_H

#include "IIterativeLinearSystemSolver.h"
#include "real_array.h"

namespace S2 { namespace ns {

/* Richardson Iteration

   From:
   http://en.wikipedia.org/wiki/Modified_Richardson_iteration
   http://math.nist.gov/iml++/ir.h.txt
   Requires A to be Positive or Negative-definite

   \todo SEEMS to fail, and original version in
   LeafDSH_Solid2D_FEM.cpp did not... but it was mostly useless
   anyway.
*/
class ILSS_RI: public IIterativeLinearSystemSolver
{
public:
    inline ILSS_RI()
    : m_r(0)
    , m_tmp0(0)
    , m_tmp1(0)
    , m_tmp2(0)
        {}
    ~ILSS_RI()
        {
            if( m_r ) delete[] m_r;
            if( m_tmp0 ) delete[] m_tmp0;
            if( m_tmp1 ) delete[] m_tmp1;
            if( m_tmp2 ) delete[] m_tmp2;
        }

    void Init( unsigned size )
        {
            IIterativeLinearSystemSolver::Init(size);
            m_r = new Real[size];
            m_tmp0 = new Real[size];
            m_tmp1 = new Real[size];
            m_tmp2 = new Real[size];
        }

    void Solve( D_Av d_Av, D_ApplyConstraints d_ApplyConstraints,
                Real* x, const Real* b,
                Real rel_epsilon, Real abs_epsilon, unsigned max_iterations )
        {
            NS_ASSERT( x != 0 && b != 0 );
            NS_ASSERT( !real_array::is_nan( x, m_Size ) );
            NS_ASSERT( !real_array::is_nan( b, m_Size ) );
            NS_ASSERT( max_iterations > 0 );
            NS_ASSERT( m_Size > 0 && m_r != 0 );

            Real lambda_max(0);
            ComputeLargestEigenvalue_Symmetric_Shifted( d_Av, m_Size,
                                                        m_tmp0, m_tmp1, m_tmp2,
                                                        25, 0.01, 0,
                                                        lambda_max );
            Real omega( 2.0f / lambda_max );

            // r <== b - A*x //KNC r <== S * (b - A*x)
            d_Av( m_r, x ); // r = A*x
            real_array::muleq_addeq( m_r, -1, b, m_Size ); // r = -1*r + b
            d_ApplyConstraints( m_r );
            NS_ASSERT( !real_array::is_nan( m_r, m_Size ) );

            Real norm_b = real_array::norm_2( b, m_Size );
            Real norm_r = real_array::norm_2( m_r, m_Size );
            if (norm_b == 0.0) norm_b = 1;
            Real rel_prec = norm_r / norm_b;

            unsigned k = 0;
            while( rel_prec > rel_epsilon
                   && k < max_iterations )
            {
                // r <== b - A*x //KNC r <== S * (b - A*x)
                d_Av( m_r, x ); // r = A*x
                real_array::muleq_addeq( m_r, -1, b, m_Size ); // r = -1*r + b
                d_ApplyConstraints( m_r );
                // y += omega * ( b - Ax )
                real_array::addeq_scaled( x, omega, m_r, m_Size );
                norm_r = real_array::norm_2( m_r, m_Size );
                rel_prec = norm_r / norm_b;
                k++;
            }
            NS_ASSERT( !real_array::is_nan( x, m_Size ) );
            // Results
            m_NumIterations = k;
            m_RelResidual = rel_prec;
            m_AbsResidual = norm_r;
        }

private:
    Real* m_r;
    Real* m_tmp0;
    Real* m_tmp1;
    Real* m_tmp2;
};

}} //namespace S2::ns

#endif // S2_NS_ILSS_RI_H
