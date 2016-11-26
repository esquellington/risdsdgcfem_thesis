#ifndef S2_NS_MATRIX_D_H
#define S2_NS_MATRIX_D_H

#include <Saphyre2/ns/Config.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace S2 { namespace ns {

typedef double RealD;

typedef boost::numeric::ublas::matrix<RealD,boost::numeric::ublas::row_major> ublas_matrix;
typedef boost::numeric::ublas::vector<RealD> ublas_vector;

class MatrixD;
class VectorD;

//! Dense Variable-Sized Matrix
/*! Dense Variable-Sized Matrix

Implemented as a wrapper of ublas::matrix by now.

\todo Make template of class <RealD>

\notes:
- Operates only with other MatrixD and VectorD, not with MatrixB
  (block-matrix) nor MatrixS (sparse-matrix)
*/
class MatrixD
{
public:

    //! Default Constructor
    MatrixD() {}

    //! Size-init constructor
    MatrixD( unsigned int rows, unsigned int cols )
    : m_ubm( rows, cols )
    {}

    //! Copy constructor
    MatrixD( const MatrixD &m )
    : m_ubm( m.m_ubm )
    {}

    //! \name Bookkeeping
    //@{
    inline int Size1() const { return m_ubm.size1(); }
    inline int Size2() const { return m_ubm.size2(); }
    inline bool IsValidIndex(int i, int j) const { return ( i>=0 && i<Size1() && j>=0 && j<Size2() ); }
    inline bool IsSquare() const { return Size1() == Size2(); }
    /*
    inline bool IsSymmetric() const { return true; }
    inline bool IsSPD() const { return true; }
    */
    //@}

    //! \name Initializers
    //@{
    void Resize( unsigned int rows, unsigned int cols );
    void Zero();
    void Identity();
    void Set( const VectorD &v, bool bTransposed );
    //@}

    //! \name Accessors
    //@{
    inline RealD &operator()(int i, int j)
    {
        assert( IsValidIndex(i,j) );
        return m_ubm(i,j);
    }
    inline const RealD &operator()(int i, int j) const
    {
        assert( IsValidIndex(i,j) );
        return m_ubm(i,j);
    }
    //@}

    //! \name Operators
    //@{
    MatrixD operator+( const MatrixD &m ) const;
    MatrixD operator-( const MatrixD &m ) const;
    MatrixD operator*( const MatrixD &m ) const;
    MatrixD operator-() const;

    void operator+=( const MatrixD &m );
    void operator-=( const MatrixD &m );
    void operator*=( RealD real );
    //@}

    //! \name Operations
    //@{
    MatrixD Transposed() const;
    //@}

    //! \name Linear system solvers
    //! Linear System solve using LU Factorization, returns prec = residual norm
    friend void SolveLU( const MatrixD &A, VectorD &x, const VectorD &b, RealD &prec );

    friend MatrixD operator*( RealD real, const MatrixD &m );

protected:
    MatrixD( const ublas_matrix &ubm ) { m_ubm = ubm; }

public:
    ublas_matrix m_ubm;
    friend class VectorD;
};


/*! Column Vector wrapper of a MatrixD. */
class VectorD: public MatrixD
{
public:
    //! Default constructor
    VectorD() {}

    //! Size-init constructor
    VectorD( unsigned int size )
    : MatrixD( size, 1 )
    {}

    //! Copy constructor from VectorD
    VectorD( const VectorD &v )
    : MatrixD( v.m_ubm )
    {}

    //! Copy constructor from MatrixD
    VectorD( const MatrixD &m )
    : MatrixD( m.m_ubm )
    { NS_ASSERT(m.m_ubm.size2() == 1); }

    //! Copy constructor from a VectorD range
    VectorD( const VectorD &v, unsigned int begin, unsigned int size )
    : MatrixD( size, 1 )
    {
        for( unsigned int i=0; i<size; i++ )
            m_ubm(i,0) = v[begin+i];
    }

    void Resize( unsigned int size ) { MatrixD::Resize(size,1); }

    //! \name Accessors
    //@{
    inline RealD &operator[](int i)
    {
        assert( IsValidIndex(i,0) );
        return m_ubm(i,0);
    }
    inline const RealD &operator[](int i) const
    {
        assert( IsValidIndex(i,0) );
        return m_ubm(i,0);
    }
    //@}

    inline int Size() const { return Size1(); }

    //!< Squared 2-norm ( == Dot(v,v) )
    inline RealD SqNorm2() const
    {
        RealD res = 0.0f;
        for(int i=0; i<Size1(); i++ )
            res += m_ubm(i,0)*m_ubm(i,0);
        return res;
    }

    //! 2-Norm
    inline RealD Norm2() const
    {
        return mal::Sqrt( SqNorm2() );
    }
};

//! Transposition as a function for readability
inline MatrixD Transpose( const MatrixD &m ) { return m.Transposed(); }

//! Dot product
inline RealD Dot( const VectorD &v1, const VectorD &v2 )
{
    RealD res = 0.0f;
    for( int i=0; i<v1.Size(); i++ )
        res += v1[i] * v2[i];
    return res;
}

void SolveCG_JWJt( const MatrixD &J, const MatrixD &W, VectorD &x, const VectorD &b,
                   RealD &prec, int &iter );

void SolveCG( const MatrixD &A, VectorD &x, const VectorD &b,
              RealD &prec, int &iter );

void SolveGS( const MatrixD &A, VectorD &x, const VectorD &b,
              RealD &prec, int &iter );

void SolveJacobi( const MatrixD &A, VectorD &x, const VectorD &b,
                  RealD &prec, int &iter );

RealD ComputeLargestEigenvalue_Symmetric( const MatrixD &m,
                                          RealD &eigen_value, VectorD &eigen_vector,
                                          unsigned max_iter, RealD epsilon );

RealD ComputeEigenValues_Symmetric( const MatrixD &m,
                                    RealD *vec_eigen_value, VectorD *vec_eigen_vector, unsigned count,
                                    unsigned max_iter, RealD epsilon );

RealD ComputeLargestAssymetry( const MatrixD &m );

template <typename OStreamT>
inline OStreamT &operator<<( OStreamT &o_stream, MatrixD &m )
{
    /*TEMP: Old format
    for( int i=0; i<m.Size1(); i++ )
    {
        o_stream << "[ ";
        for( int j=0; j<m.Size2(); j++ ) o_stream << m(i,j) << " ";
        o_stream << "]\n";
    }
    */
    // Ouput in Octave format for easy copy-paste
    o_stream << "[ ";
    for( int i=0; i<m.Size1(); i++ )
    {
        for( int j=0; j<m.Size2(); j++ )
            if( j==m.Size2()-1 )
                o_stream << m(i,j);
            else
                o_stream << m(i,j) << ", ";
        if( i==m.Size1()-1 )
            o_stream << " ]\n";
        else
            o_stream << " ;\n";
    }
    return o_stream;
}

}} //namespace S2::ns

#endif // S2_NS_MATRIX_D_H
