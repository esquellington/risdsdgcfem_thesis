#include "MatrixD.h"
#include <Mal/GRandom.h>

//#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

#define S2_NS_SOLVECG_DEFAULT_MAX_ITER 50
#define S2_NS_SOLVECG_DEFAULT_PRECISION 1.0e-5

namespace S2 { namespace ns {

using namespace boost::numeric::ublas;

//---- Initializers
void MatrixD::Resize( unsigned int rows, unsigned int cols )
{
    m_ubm.resize( rows, cols, false );
}
void MatrixD::Zero() { m_ubm.assign( zero_matrix<RealD>( Size1(), Size2() ) ); }
void MatrixD::Identity()
{
    NS_ASSERT(Size1() == Size2());
    m_ubm.assign( identity_matrix<RealD>( Size1(), Size2() ) );
}

void MatrixD::Set( const VectorD &v, bool bTransposed )
{
    if( bTransposed )
    {
        NS_ASSERT( Size2() == v.Size() );
        for( int i=0; i<v.Size(); i++ )
            m_ubm(0,i) = v[i];
    }
    else
    {
        NS_ASSERT( Size1() == v.Size() );
        for( int i=0; i<v.Size(); i++ )
            m_ubm(i,0) = v[i];
    }
}

//---- Operators
MatrixD MatrixD::operator+( const MatrixD &m ) const
{
    NS_ASSERT( Size1() == m.Size1() && Size2() == m.Size2() );
    return MatrixD( m_ubm + m.m_ubm );
}
MatrixD MatrixD::operator-( const MatrixD &m ) const
{
    NS_ASSERT( Size1() == m.Size1() && Size2() == m.Size2() );
    return MatrixD( m_ubm - m.m_ubm );
}
MatrixD MatrixD::operator*( const MatrixD &m ) const
{
    NS_ASSERT( Size2() == m.Size1() );
    return MatrixD( prod( m_ubm, m.m_ubm ) );
}
MatrixD MatrixD::operator-() const
{
    return MatrixD( -m_ubm );
}

void MatrixD::operator+=( const MatrixD &m )
{
    NS_ASSERT( Size1() == m.Size1() && Size2() == m.Size2() );
    m_ubm += m.m_ubm;
}

void MatrixD::operator-=( const MatrixD &m )
{
    NS_ASSERT( Size1() == m.Size1() && Size2() == m.Size2() );
    m_ubm -= m.m_ubm;
}

void MatrixD::operator*=( RealD real )
{
    m_ubm *= real;
}

MatrixD operator*( RealD real, const MatrixD &m )
{
    return MatrixD( real * m.m_ubm );
}

//---- Operations
MatrixD MatrixD::Transposed() const
{
    return MatrixD( trans( m_ubm ) );
}

/*! Linear System solve using LU Factorization

  \warning Fails if /fp:fast is enabled!!
*/
void SolveLU( const MatrixD &A, VectorD &x, const VectorD &b, RealD &prec )
{
    ublas_matrix LU( A.m_ubm );
    permutation_matrix<std::size_t> P( A.Size1() );
    int res = lu_factorize(LU,P);
    NS_ASSERT( 0 == res );
    ublas_vector vx(b.Size());
    for( int i=0; i<b.Size(); i++)
        vx[i] = b[i];
    lu_substitute( LU, P, vx );
    for( int i=0; i<b.Size(); i++)
        x[i] = vx[i];
    prec = VectorD(A*x-b).Norm2(); //\todo THIS IS VERY INEFFICIENT!!
}

/*! Linear System solve using Conjugate Gradient
  Solves:
     J*W*Jt*x = b
  Where:
     J*W*Jt is Symmetric and Positive Definite

  Uses only Matrix x Vector products to exploit sparsity of J and W.

  \todo Optimize with TransposedMultiply() instead of J.Transposed()
  and 3-param multiplication functions instead of overloaded operator*()!

  \note Loosely based in "An Introduction to the Conjugate Gradient
  Method Without the Agonyzing Pain", Chapter B2, Page 50.
*/
void SolveCG_JWJt( const MatrixD &J, const MatrixD &W, VectorD &x, const VectorD &b,
                   RealD &prec, int &iter )
{
    NS_ASSERT( x.Size() == b.Size() );
    //NS_ASSERT( (J*W*Jt).IsSPD() );
    for ( int i = 0; i < b.Size(); ++i ) NS_ASSERT( !mal::IsNaN(b[i]) );
    for ( int i = 0; i < J.Size1(); ++i ) for ( int j = 0; j < J.Size2(); ++j ) NS_ASSERT( !mal::IsNaN(J(i,j)) );
    for ( int i = 0; i < W.Size1(); ++i ) for ( int j = 0; j < W.Size2(); ++j ) NS_ASSERT( !mal::IsNaN(W(i,j)) );

    int max_iter = (iter > 0) ? iter : S2_NS_SOLVECG_DEFAULT_MAX_ITER;
    RealD epsilon = (prec > 0) ? prec : S2_NS_SOLVECG_DEFAULT_PRECISION;

    MatrixD Jt( J.Transposed() );
    VectorD r( b - J*( W*( Jt*x ) ) );
    VectorD d( r );
    VectorD q( b.Size() );

    RealD psi_0 = r.SqNorm2(); //psi <== r^T * r;
    RealD psi = psi_0;
    RealD psi_ant = psi_0;
    int k = 0;
    while( psi > epsilon*epsilon*psi_0 && k < max_iter )
    {
        // q <== A*d
        q = J*( W*( Jt*d ) ); //!< \todo OPTIMIZE!!
        RealD alpha = psi/Dot( d, q );
        x += alpha*d;

        // r <== r - alpha*q
        q *= alpha;
        r -= q;

        /* To Avoid residual error accumulation (see CG paper)
        if( 0==k%50 )
            r = b - J*( W*( Jt*x ) );
        */

        psi_ant = psi;
        psi = r.SqNorm2(); // psi <== r^T * r;

        /* To Avoid cancellation error (see CG paper)
        if( psi <= epsilon*epsilon*psi_0 )
        {
            r = b - J*( W*( Jt*x ) );
            psi = r.SqNorm2();
        }
        */

        RealD beta = (psi/psi_ant);

        // d <== r + beta*d
        d *= beta;
        d += r;

        k++;
    }

    // write results
    iter = k;
    prec = mal::Sqrt(psi);
}

void SolveCG( const MatrixD &A, VectorD &x, const VectorD &b,
              RealD &prec, int &iter )
{
    NS_ASSERT( A.IsSquare() );
    NS_ASSERT( A.Size2() == b.Size() );
    NS_ASSERT( x.Size() == b.Size() );

    // b_paranoid mode
    //NS_ASSERT( IsSPD() );
    for ( int i = 0; i < b.Size(); ++i ) NS_ASSERT( !mal::IsNaN(b[i]) );
    for ( int i = 0; i < A.Size1(); ++i ) for ( int j = 0; j < A.Size2(); ++j ) NS_ASSERT( !mal::IsNaN(A(i,j)) );

    int max_iter = (iter > 0) ? iter : S2_NS_SOLVECG_DEFAULT_MAX_ITER;
    RealD epsilon = (prec > 0) ? prec : S2_NS_SOLVECG_DEFAULT_PRECISION;

    VectorD r( b - A*x );
    VectorD d( r );
    VectorD q( b.Size() );

    RealD psi_0 = r.SqNorm2(); //psi <== r^T * r;
    RealD psi = psi_0;
    RealD psi_ant = psi_0;
    int k = 0;
    while( psi > epsilon*epsilon*psi_0 && k < max_iter )
    {
        // q <== A*d
        q = A*d;
        RealD alpha = psi/Dot( d, q );
        x += alpha*d;

        // r <== r - alpha*q
        q *= alpha;
        r -= q;

        /* To Avoid residual error accumulation (see CG paper)
        if( 0==k%50 )
            r = b - A*x;
        */

        psi_ant = psi;
        psi = r.SqNorm2(); // psi <== r^T * r;

        /* To Avoid cancellation error (see CG paper)
        if( psi <= epsilon*epsilon*psi_0 )
        {
            r = b - A*x;
            psi = r.SqNorm2();
        }
        */

        RealD beta = (psi/psi_ant);

        // d <== r + beta*d
        d *= beta;
        d += r;

        k++;
    }

    // write results
    iter = k;
    prec = mal::Sqrt(psi);
}

/*! Linear System solve using Gauss-Seidel
TEMPORAL IMPLEMENTATION... this comes from OpenTissue and is copy-pasted... thus ugly, just for testing
\sa http://image.diku.dk/svn/OpenTissue/archieve/sunegn/OpenTissue/math/boost_matrix_solvers/
*/
void SolveGS( const MatrixD &A, VectorD &x, const VectorD &b,
              RealD &prec, int &iter )
{
    NS_ASSERT( A.IsSquare() );
    NS_ASSERT( A.Size2() == b.Size() );
    NS_ASSERT( x.Size() == b.Size() );

    // b_paranoid mode
    for ( int i = 0; i < b.Size(); ++i ) NS_ASSERT( !mal::IsNaN(b[i]) );
    for ( int i = 0; i < A.Size1(); ++i ) for ( int j = 0; j < A.Size2(); ++j ) NS_ASSERT( !mal::IsNaN(A(i,j)) );
    //NS_ASSERT( IsSPD() || IsDiagonallyDominant() );

    int max_iter = (iter > 0) ? iter : S2_NS_SOLVECG_DEFAULT_MAX_ITER;
    RealD epsilon = (prec > 0) ? prec : S2_NS_SOLVECG_DEFAULT_PRECISION;

    int n = x.Size();
    int k = 0;
    for ( int i = 0; i < x.Size(); ++i ) x[i] = 1.0f/x.Size(); //\todo should receive approx solution
    while ( k < max_iter )
    {
        for ( int i = 0; i < n; ++i )
        {
            x[i] = b[i];
            //NS_ASSERT( !mal::IsNaN(x[i]) );
            for ( int j = 0; j < i; ++j ) x[i] -= A( i, j ) * x[ j ];
            //NS_ASSERT( !mal::IsNaN(x[i]) );
            for ( int j = i + 1; j < n; ++j ) x[i] -= A( i, j ) * x[ j ];
            //NS_ASSERT( !mal::IsNaN(x[i]) );
            if( mal::Abs(A( i, i )) > 0.00001f ) x[i] /= A( i, i );
            //NS_ASSERT( mal::Abs(A( i, i )) > 0.000001f ); //If this happens, we're doomed
            //NS_ASSERT( !mal::IsNaN(x[i]) );
        }
        k++;
    }
    iter = k;
    prec = VectorD(A*x-b).Norm2(); //\todo THIS IS VERY INEFFICIENT!!
}

void SolveJacobi( const MatrixD &A, VectorD &x, const VectorD &b,
                  RealD &prec, int &iter )
{
    NS_ASSERT( A.IsSquare() );
    NS_ASSERT( A.Size2() == b.Size() );
    NS_ASSERT( x.Size() == b.Size() );

    // b_paranoid mode
    for ( int i = 0; i < b.Size(); ++i ) NS_ASSERT( !mal::IsNaN(b[i]) );
    for ( int i = 0; i < A.Size1(); ++i ) for ( int j = 0; j < A.Size2(); ++j ) NS_ASSERT( !mal::IsNaN(A(i,j)) );
    //NS_ASSERT( IsSPD() || IsDiagonallyDominant() );

    int max_iter = (iter > 0) ? iter : S2_NS_SOLVECG_DEFAULT_MAX_ITER;
    RealD epsilon = (prec > 0) ? prec : S2_NS_SOLVECG_DEFAULT_PRECISION;

    int n = x.Size();
    int k = 0;
    for ( int i = 0; i < x.Size(); ++i ) x[i] = 1.0f/x.Size(); //\todo should receive approx solution
    VectorD x_old( x );
    while ( k < max_iter )
    {
        for ( int i = 0; i < n; ++i )
        {
            x[i] = b[i];
            //NS_ASSERT( !mal::IsNaN(x[i]) );
            for ( int j = 0; j < i; ++j ) x[i] -= A( i, j ) * x_old[ j ];
            //NS_ASSERT( !mal::IsNaN(x[i]) );
            for ( int j = i + 1; j < n; ++j ) x[i] -= A( i, j ) * x_old[ j ];
            //NS_ASSERT( !mal::IsNaN(x[i]) );
            if( mal::Abs(A( i, i )) > 0.00001f ) x[i] = x_old[i] / A( i, i );
            //NS_ASSERT( mal::Abs(A( i, i )) > 0.000001f ); //If this happens, we're doomed
            NS_ASSERT( !mal::IsNaN(x[i]) );
        }
        x_old = x;
        k++;
    }
    iter = k;
    prec = VectorD(A*x-b).Norm2(); //\todo THIS IS VERY INEFFICIENT!!
}

/*! Compute largest (eigen_value,eigen_vector) pair of a real
    symmetric matrix.

    Returns absolute error of the last two eigen_value approximations
    computed, which will be ideally < epsilon unless max_iter is hit
    or the algorithm fails to converge.

    Returns < 0 if no significant result can be found (may be due to
    random seed vector)

    \note May return 0 eigenvalue/eigenvector pair for M = 0 AND for
    random v1 perpendicular to all M rows, which CAN happen if M is not
    full-rank.
*/
RealD ComputeLargestEigenvalue_Symmetric( const MatrixD &m,
                                          RealD &eigen_value, VectorD &eigen_vector,
                                          unsigned max_iter, RealD epsilon )
{
    // b_paranoid mode
    NS_ASSERT( m.IsSquare() );
    //NS_ASSERT( m.IsSymmetric() );
    for ( int i = 0; i < m.Size1(); ++i ) for ( int j = 0; j < m.Size2(); ++j ) NS_ASSERT( !mal::IsNaN(m(i,j)) );

    VectorD v0( m.Size1() );
    RealD e0(0);
    // Gen random non-null v1
    VectorD v1( m.Size1() );
    do
    {
        for( int i=0; i<v0.Size(); i++ )
            v1[i] = mal::RandomF<RealD>( -1, 1 );
    }
    while( v1.Norm2() < epsilon );
    RealD e1( v1.Norm2() );
    unsigned int iter(0);
    do
    {
        e0 = e1;
        //if( e1 <= epsilon ) std::cout << "M=" << m << "v1=" << v1 << "e1=" << e1 << std::endl;
        v0 = v1; v0 *= mal::Rcp(e1);
        v1 = m * v0;
        e1 = v1.Norm2();
        iter++;
    }
    while( iter < max_iter
           && e1 > epsilon //Ensure non-collapsing eigenvector
           && mal::Abs(e1-e0) > epsilon );
    // Compute final eigen_vector //\todo ENSURE THIS WORKS FOR m = 0!!
    if( e1 > epsilon )
    {
        eigen_vector = v1; eigen_vector *= mal::Rcp(e1);
        // Compute eigen_value sign and return
        for( int i=0; i<v0.Size(); i++ )
            if( v0[i] != RealD(0) && v1[i] != RealD(0) ) //\todo should use eigen_value instead of v1, but sign will be the same
            {
                eigen_value = ( v0[i] * v1[i] > RealD(0) ) ? e1 : -e1;
                return mal::Abs(e1-e0);
            }
        // No non-zero component found in iterated eigen_vector, no valid eigen_value
        eigen_value = 0;
        return RealD(0);
    }
    else
    {
        // eigen_vector is almost zero, assume 0 largest eigen_value
        eigen_vector.Zero();
        eigen_value = 0;
        return RealD(0);
    }
}

RealD ComputeEigenValues_Symmetric( const MatrixD &m,
                                    RealD *vec_eigen_value, VectorD *vec_eigen_vector, unsigned count,
                                    unsigned max_iter, RealD epsilon )
{
    //NS_ASSERT( m.IsSymmetric() );
    MatrixD A(m);
    RealD max_error(0);
    for( unsigned int c=0; c<count; c++ )
    {
        // Compute largest
        RealD eps_error = ComputeLargestEigenvalue_Symmetric( A, vec_eigen_value[c], vec_eigen_vector[c], max_iter, epsilon );
        if( eps_error < RealD(0) ) return eps_error;
        else if( max_error < eps_error ) max_error = eps_error;
        // Deflate largest
        if( c < count-1 )
            for( int i=0; i<m.Size1(); i++ )
                for( int j=0; j<m.Size2(); j++ )
                    A(i,j) -= vec_eigen_value[c] * vec_eigen_vector[c][i] * vec_eigen_vector[c][j];
    }
    return max_error;
}

RealD ComputeLargestAssymetry( const MatrixD &m )
{
    Real las(0);
    for( int i=0; i<m.Size1(); i++ )
        for( int j=i; j<m.Size2(); j++ )
            las = mal::Max<Real>( las, mal::Abs( m(i,j) - m(j,i) ) );
    return las;
}

}} //namespace S2::ns
