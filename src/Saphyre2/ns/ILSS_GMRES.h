#ifndef S2_NS_ILSS_GMRES_H
#define S2_NS_ILSS_GMRES_H

#include "IIterativeLinearSystemSolver.h"
#include "real_array.h"
#include "MatrixD.h"

namespace S2 { namespace ns {

/* GMRES

   From the paper "GMRES: A GENERALIZED MINIMAL RESIDUAL ALGORITHM FOR
   SOLVING NONSYMMETRIC LINEAR SYSTEMS*, Y.Saad & M.H.Schultz, 1986

   \note Code adapted from http://math.nist.gov/iml++/gmres.h.txt

   \todo MOVE CODE TO CPP
*/
class ILSS_GMRES: public IIterativeLinearSystemSolver
{
public:
    finline ILSS_GMRES()
    : m_MaxRestarts(0)
    , m_r(0)
    , m_d(0)
    , m_q(0)
    , m_tmp(0)
    , m_w(0)
    , m_s(0)
    , m_cs(0)
    , m_sn(0)
    , m_V(0)
        {}
    ~ILSS_GMRES()
        {
            if( m_r ) delete[] m_r;
            if( m_d ) delete[] m_d;
            if( m_q ) delete[] m_q;
            if( m_tmp ) delete[] m_tmp;
            if( m_w ) delete[] m_w;
            if( m_s ) delete[] m_s;
            if( m_cs ) delete[] m_cs;
            if( m_sn ) delete[] m_sn;
            if( m_V )
            {
                for( unsigned int it_v=0; it_v<m_MaxRestarts+1; it_v++ ) delete[] m_V[it_v];
                delete[] m_V;
            }
        }

    void Init( unsigned size, unsigned max_restarts )
        {
            IIterativeLinearSystemSolver::Init(size);
            m_r = new Real[size];
            m_d = new Real[size];
            m_q = new Real[size];
            m_tmp = new Real[size];
            m_w = new Real[size];
            SetMaxRestarts( max_restarts );
        }

    void SetMaxRestarts( unsigned max_restarts )
        {
            //no point in max_restarts >= size, it seems from GMRES paper
            unsigned new_max_restarts = mal::Min( max_restarts, m_Size-1 );
            if( new_max_restarts > m_MaxRestarts )
            {
                if( m_MaxRestarts > 0 )
                {
                    delete[] m_s;
                    delete[] m_cs;
                    delete[] m_sn;
                    for( unsigned int it_v=0; it_v<m_MaxRestarts+1; it_v++ ) delete[] m_V[it_v];
                    delete[] m_V;
                }
                // Vector s(m+1), cs(m+1), sn(m+1)
                m_s = new Real[new_max_restarts+1];
                m_cs = new Real[new_max_restarts+1];
                m_sn = new Real[new_max_restarts+1];
                // Alloc V and H
                m_V = new Real*[new_max_restarts+1];
                for( unsigned int it_v=0; it_v<new_max_restarts+1; it_v++ ) m_V[it_v] = new Real[m_Size];
                //\todo Matrix H is upper-triangular, it seems... consider better allocation
                m_H.Resize( new_max_restarts+1, new_max_restarts+1 );
            }
            m_MaxRestarts = new_max_restarts;
            //NS_LOG_WARNING( "SetMaxRestarts = %d", m_MaxRestarts );
        }

    void Update_H( Real* x,
                   Real* tmp, //intermediate
                   int k,
                   const ns::MatrixD& H, const Real* s,
                   Real** V, //\todo should be const
                   unsigned int size )
        {
            ns::real_array::assign( tmp, s, k+1 );
            // Backsolve:
            for( int i = k; i >= 0; i-- )
            {
                tmp[i] /= H(i,i);
                for( int j = i - 1; j >= 0; j-- )
                    tmp[j] -= H(j,i) * tmp[i];
            }
            for( int j = 0; j <= k; j++ )
                ns::real_array::addeq_scaled( x, tmp[j], V[j], size ); //x += V[j] * tmp[j];
        }


    inline void GeneratePlaneRotation( const double& dx, const double& dy, Real& cs, Real& sn )
        {
            if (dy == 0.0) {
                cs = 1.0;
                sn = 0.0;
            } else if (mal::Abs(dy) > mal::Abs(dx)) {
                Real temp = dx / dy;
                sn = 1.0 / mal::Sqrt( 1.0 + temp*temp );
                cs = temp * sn;
            } else {
                Real temp = dy / dx;
                cs = 1.0 / mal::Sqrt( 1.0 + temp*temp );
                sn = temp * cs;
            }
        }

    inline void ApplyPlaneRotation( Real& dx, Real& dy, const Real& cs, const Real& sn )
        {
            Real temp = cs * dx + sn * dy;
            dy = -sn * dx + cs * dy;
            dx = temp;
        }

    inline void ApplyPlaneRotation( double& dx, double& dy, const Real& cs, const Real& sn )
        {
            Real temp = cs * dx + sn * dy;
            dy = -sn * dx + cs * dy;
            dx = temp;
        }

    void Solve( D_Av d_Av, D_ApplyConstraints d_ApplyConstraints,
                Real* x, const Real* b,
                Real rel_epsilon, Real abs_epsilon, unsigned max_iterations )
        {
            NS_ASSERT( x != 0 && b != 0 );
            NS_ASSERT( !ns::real_array::is_nan( x, m_Size ) );
            NS_ASSERT( !ns::real_array::is_nan( b, m_Size ) );
            NS_ASSERT( max_iterations > 0 );
            NS_ASSERT( m_Size > 0 && m_MaxRestarts > 0
                       && m_r != 0 && m_d != 0 && m_q != 0
                       && m_tmp != 0 && m_w != 0 && m_s != 0
                       && m_cs != 0 && m_sn != 0 && m_V != 0 );

            // r <== b - A*y //KNC r <== S * (b - A*y)
            d_Av( m_r, x ); // r = A*x
            ns::real_array::muleq_addeq( m_r, -1, b, m_Size ); // r = -1*r + b
            NS_ASSERT( !ns::real_array::is_nan( m_r, m_Size ) );
            d_ApplyConstraints( m_r );
            NS_ASSERT( !ns::real_array::is_nan( m_r, m_Size ) );

            Real norm_b = ns::real_array::norm_2( b, m_Size );
            Real norm_r = ns::real_array::norm_2( m_r, m_Size ); //==beta
            if( norm_b == 0.0 ) norm_b = 1;
            Real rel_residual = norm_r / norm_b; //==resid

#ifdef __S2_NS_ENABLE_TRACE_ILSS_RESIDUAL
            NS_LOG("GMRES[%4d] delta = %f", -1, norm_r );
#endif

            unsigned num_iter = 0;
            while( rel_residual > rel_epsilon
                   && norm_r > abs_epsilon //absolute precision condition, not in original alg
                   && num_iter < max_iterations )
            {
                ns::real_array::assign_scaled( m_V[0], mal::Rcp(norm_r), m_r, m_Size ); //v[0] = r * (1.0 / beta);    // ??? r / beta
                ns::real_array::zero( m_s, m_MaxRestarts+1 ); //s = 0.0;
                m_s[0] = norm_r; //s(0) = beta;
                unsigned int i;
                for( i = 0; i < m_MaxRestarts && num_iter < max_iterations; i++, num_iter++ )
                {
                    //w = M.solve(A * v[i]);
                    d_Av( m_w, m_V[i] ); // w = A*v[i]
                    d_ApplyConstraints( m_w ); //\todo NOT SURE
                    NS_ASSERT( !ns::real_array::is_nan( m_w, m_Size ) );
                    for( unsigned int k = 0; k <= i; k++ )
                    {
                        m_H(k,i) = ns::real_array::dot( m_w, m_V[k], m_Size );
                        ns::real_array::subeq_scaled( m_w, m_H(k,i), m_V[k], m_Size ); //w -= H(k,i) * V[k];
                    }
                    m_H(i+1,i) = ns::real_array::norm_2( m_w, m_Size );
                    ns::real_array::assign_scaled( m_V[i+1], mal::Rcp( m_H(i+1,i) ), m_w, m_Size ); //v[i+1] = w * (1.0 / H(i+1, i));

                    for( unsigned int k = 0; k < i; k++ )
                        ApplyPlaneRotation( m_H(k,i), m_H(k+1,i), m_cs[k], m_sn[k] );

                    GeneratePlaneRotation( m_H(i,i), m_H(i+1,i), m_cs[i], m_sn[i] );
                    ApplyPlaneRotation( m_H(i,i), m_H(i+1,i), m_cs[i], m_sn[i] );
                    ApplyPlaneRotation( m_s[i], m_s[i+1], m_cs[i], m_sn[i] );

                    rel_residual = mal::Abs(m_s[i+1]) / norm_b;

#ifdef __S2_NS_ENABLE_TRACE_ILSS_RESIDUAL
                    NS_LOG("GMRES[%4u] delta = %f", num_iter, rel_residual * norm_b );
#endif

                    Real abs_residual(rel_residual * norm_b);
                    if( rel_residual < rel_epsilon
                        || abs_residual < abs_epsilon ) //absolute precision condition, not in original alg
                    {
                        Update_H( x, m_tmp, i, m_H, m_s, m_V, m_Size );
                        m_RelResidual = rel_residual;
                        m_AbsResidual = abs_residual; //\todo Not sure... but norm_r may not be up-to-date here
                        m_NumIterations = num_iter;
                        return; //converged
                    }
                }
                Update_H( x, m_tmp, i - 1, m_H, m_s, m_V, m_Size );
                // r <== b - A*y //KNC r <== S * (b - A*y)
                d_Av( m_r, x ); // r = A*x
                ns::real_array::muleq_addeq( m_r, -1, b, m_Size ); // r = -1*r + b
                d_ApplyConstraints( m_r );
                NS_ASSERT( !ns::real_array::is_nan( m_r, m_Size ) );
                norm_r = ns::real_array::norm_2( m_r, m_Size );
                rel_residual = norm_r / norm_b;

#ifdef __S2_NS_ENABLE_TRACE_ILSS_RESIDUAL
                NS_LOG("GMRES[%4u] delta = %f", num_iter, norm_r );
#endif

                if( rel_residual < rel_epsilon
                    || norm_r < abs_epsilon ) //absolute precision condition, not in original alg
                {
                    m_RelResidual = rel_residual;
                    m_AbsResidual = norm_r;
                    m_NumIterations = num_iter;
                    return; //converged
                }
            }
            m_RelResidual = rel_residual;
            m_AbsResidual = norm_r;
            m_NumIterations = num_iter;
            return; //not converged
        }

private:
    unsigned m_MaxRestarts;
    Real* m_r;
    Real* m_d;
    Real* m_q;
    Real* m_tmp;
    Real* m_w;
    Real* m_s;
    Real* m_cs;
    Real* m_sn;
    Real** m_V;
    MatrixD m_H;
};

}} //namespace S2::ns

#endif // S2_NS_ILSS_GMRES_H
