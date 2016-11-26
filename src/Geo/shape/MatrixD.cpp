
#include "MatrixD.h"

//#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

#define S2_LA_SOLVECG_DEFAULT_MAX_ITER 50
#define S2_LA_SOLVECG_DEFAULT_PRECISION 1.0e-5

namespace S2 { namespace LA {

using namespace boost::numeric::ublas;

//---- Initializers
void MatrixD::Resize( unsigned int rows, unsigned int cols )
{
    m_ubm.resize( rows, cols, false );
}
void MatrixD::Zero() { m_ubm.assign( zero_matrix<Real>( Size1(), Size2() ) ); }
void MatrixD::Identity()
{
    assert(Size1() == Size2());
    m_ubm.assign( identity_matrix<Real>( Size1(), Size2() ) );
}

void MatrixD::Set( const VectorD &v, bool bTransposed )
{
    if( bTransposed )
    {
        assert( Size2() == v.Size() );
        for( int i=0; i<v.Size(); i++ )
            m_ubm(0,i) = v[i];
    }
    else
    {
        assert( Size1() == v.Size() );
        for( int i=0; i<v.Size(); i++ )
            m_ubm(i,0) = v[i];                         
    }
}

//---- Operators
MatrixD MatrixD::operator+( const MatrixD &m ) const
{
    assert( Size1() == m.Size1() && Size2() == m.Size2() );
    return MatrixD( m_ubm + m.m_ubm );
}
MatrixD MatrixD::operator-( const MatrixD &m ) const
{
    assert( Size1() == m.Size1() && Size2() == m.Size2() );
    return MatrixD( m_ubm - m.m_ubm );    
}
MatrixD MatrixD::operator*( const MatrixD &m ) const
{
    assert( Size2() == m.Size1() );
    return MatrixD( prod( m_ubm, m.m_ubm ) );
}
MatrixD MatrixD::operator-() const
{
    return MatrixD( -m_ubm );
}

void MatrixD::operator+=( const MatrixD &m )
{
    assert( Size1() == m.Size1() && Size2() == m.Size2() );
    m_ubm += m.m_ubm;
}

void MatrixD::operator-=( const MatrixD &m )
{
    assert( Size1() == m.Size1() && Size2() == m.Size2() );
    m_ubm -= m.m_ubm;
}
    
void MatrixD::operator*=( Real real )
{
    m_ubm *= real;
}

MatrixD operator*( Real real, const MatrixD &m )
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
void SolveLU( const MatrixD &A, VectorD &x, const VectorD &b )
{    
    ublas_matrix LU( A.m_ubm );
    permutation_matrix<std::size_t> P( A.Size1() );
    int res = lu_factorize(LU,P);
    assert( 0 == res );
    ublas_vector vx(b.Size());
    for( int i=0; i<b.Size(); i++)
        vx[i] = b[i];
    lu_substitute( LU, P, vx );
    for( int i=0; i<b.Size(); i++)
        x[i] = vx[i];
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
void SolveCG_JWJt( const MatrixD &J, const MatrixD &W, const VectorD &b,
                   VectorD &x,
                   Real &prec, int &iter )
{    
    assert( x.Size() == b.Size() );

    int max_iter = (iter > 0) ? iter : S2_LA_SOLVECG_DEFAULT_MAX_ITER;
    Real epsilon = (prec > 0) ? prec : S2_LA_SOLVECG_DEFAULT_PRECISION;

    MatrixD Jt( J.Transposed() );
    VectorD r( b - J*( W*( Jt*x ) ) );
    VectorD d( r );
    VectorD q( b.Size() );

    Real psi_0 = r.SqNorm2(); //psi <== r^T * r;
    Real psi = psi_0;
    Real psi_ant = psi_0;
    int k = 0;
    while( psi > epsilon*epsilon*psi_0 && k < max_iter )
    {
        // q <== A*d
        q = J*( W*( Jt*d ) ); //!< \todo OPTIMIZE!!
        Real alpha = psi/Dot( d, q );
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
        
        Real beta = (psi/psi_ant);
        
        // d <== r + beta*d
        d *= beta;
        d += r;

        k++;
    }

    // write results
    iter = k;
    prec = mal::Sqrt(psi);
}

void SolveCG( const MatrixD &A, const VectorD &b,
              VectorD &x,
              Real &prec, int &iter )
{    
    assert( x.Size() == b.Size() );

    int max_iter = (iter > 0) ? iter : S2_LA_SOLVECG_DEFAULT_MAX_ITER;
    Real epsilon = (prec > 0) ? prec : S2_LA_SOLVECG_DEFAULT_PRECISION;

    VectorD r( b - A*x );
    VectorD d( r );
    VectorD q( b.Size() );

    Real psi_0 = r.SqNorm2(); //psi <== r^T * r;
    Real psi = psi_0;
    Real psi_ant = psi_0;
    int k = 0;
    while( psi > epsilon*epsilon*psi_0 && k < max_iter )
    {
        // q <== A*d
        q = A*d;
        Real alpha = psi/Dot( d, q );
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
        
        Real beta = (psi/psi_ant);
        
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
void SolveGS( const MatrixD &A, const VectorD &b,
              VectorD &x, Real &prec, int &iter )
{
    // b_paranoid mode
    for ( int i = 0; i < b.Size(); ++i ) GEO_ASSERT( !mal::IsNaN(b[i]) );
    for ( int i = 0; i < A.Size1(); ++i ) for ( int j = 0; j < A.Size2(); ++j ) GEO_ASSERT( !mal::IsNaN(A(i,j)) );
    //GEO_ASSERT( IsSPD() || IsDiagonallyDominant() );
    
    int max_iter = (iter > 0) ? iter : S2_LA_SOLVECG_DEFAULT_MAX_ITER;
    Real epsilon = (prec > 0) ? prec : S2_LA_SOLVECG_DEFAULT_PRECISION;
    
    int n = x.Size();
    int k = 0;
    for ( int i = 0; i < x.Size(); ++i ) x[i] = 1.0f/x.Size(); //\todo should receive approx solution
    while ( k < max_iter )
    {
        for ( int i = 0; i < n; ++i )
        {
            x[i] = b[i];
            //GEO_ASSERT( !mal::IsNaN(x[i]) );
            for ( int j = 0; j < i; ++j ) x[i] -= A( i, j ) * x[ j ];
            //GEO_ASSERT( !mal::IsNaN(x[i]) );
            for ( int j = i + 1; j < n; ++j ) x[i] -= A( i, j ) * x[ j ];
            //GEO_ASSERT( !mal::IsNaN(x[i]) );
            if( mal::Abs(A( i, i )) > 0.00001f ) x[i] /= A( i, i );            
            //GEO_ASSERT( mal::Abs(A( i, i )) > 0.000001f ); //If this happens, we're doomed
            //GEO_ASSERT( !mal::IsNaN(x[i]) );
        }
        k++;
    }
    iter = k;
    prec = -1;
}

void SolveJacobi( const MatrixD &A, const VectorD &b,
                  VectorD &x, Real &prec, int &iter )
{
    // b_paranoid mode
    for ( int i = 0; i < b.Size(); ++i ) GEO_ASSERT( !mal::IsNaN(b[i]) );
    for ( int i = 0; i < A.Size1(); ++i ) for ( int j = 0; j < A.Size2(); ++j ) GEO_ASSERT( !mal::IsNaN(A(i,j)) );
    //GEO_ASSERT( IsSPD() || IsDiagonallyDominant() );
    
    int max_iter = (iter > 0) ? iter : S2_LA_SOLVECG_DEFAULT_MAX_ITER;
    Real epsilon = (prec > 0) ? prec : S2_LA_SOLVECG_DEFAULT_PRECISION;
    
    int n = x.Size();
    int k = 0;
    for ( int i = 0; i < x.Size(); ++i ) x[i] = 1.0f/x.Size(); //\todo should receive approx solution
    VectorD x_old( x );
    while ( k < max_iter )
    {
        for ( int i = 0; i < n; ++i )
        {
            x[i] = b[i];
            //GEO_ASSERT( !mal::IsNaN(x[i]) );
            for ( int j = 0; j < i; ++j ) x[i] -= A( i, j ) * x_old[ j ];
            //GEO_ASSERT( !mal::IsNaN(x[i]) );
            for ( int j = i + 1; j < n; ++j ) x[i] -= A( i, j ) * x_old[ j ];
            //GEO_ASSERT( !mal::IsNaN(x[i]) );
            if( mal::Abs(A( i, i )) > 0.00001f ) x[i] = x_old[i] / A( i, i );
            //GEO_ASSERT( mal::Abs(A( i, i )) > 0.000001f ); //If this happens, we're doomed
            GEO_ASSERT( !mal::IsNaN(x[i]) );
        }
        x_old = x;
        k++;
    }
    iter = k;
    prec = -1;
}


}} //namespace S2::LA
