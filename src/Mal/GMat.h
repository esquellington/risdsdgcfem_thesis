#ifndef MAL_GMAT_H
#define MAL_GMAT_H

#include <Mal/Config.h>
#include <Mal/RealUtils.h>
#include <algorithm>

namespace mal
{

template<typename T, int NR, int NC>
class GMat
{
public:
    typedef T real_type;
    static const unsigned cRows = NR;
    static const unsigned cCols = NC;
    static const int cNumElements = NR*NC;
    static const unsigned size_in_reals = NR*NC;

    finline static GMat Zero() { return GMat(T(0)); }
    finline static GMat Identity() { GMat identity(T(0));
                                     for(unsigned int i=0; i<cRows; i++) identity(i,i) = T(1);
                                     return identity; }
public:
    //!\name Construction
    //@{
    finline GMat() {}
    finline explicit GMat( T val ) { for(int k=0; k<cNumElements; k++) GetPtr()[k] = val; }
    template<typename S>
    finline GMat( const GMat<S,NR,NC> &m ) { for(int i=0; i<NR; i++) for(int j=0; j<NC; j++) (*this)(i,j) = T(m(i,j)); }
    finline explicit GMat( const T* real_array ) { for(int k=0; k<cNumElements; k++) GetPtr()[k] = real_array[k]; }
    //@}


    //\name Square matrix construction
    //@{
    // 2x2
    template <typename S00, typename S01,
              typename S10, typename S11>
    finline explicit GMat( S00 m00, S01 m01,
                           S10 m10, S11 m11 ) { MAL_STATIC_ASSERT(NC==2 && NR==2);
                                                data[0][0] = T(m00); data[0][1] = T(m01);
                                                data[1][0] = T(m10); data[1][1] = T(m11); }
    // 3x3
    template <typename S00, typename S01, typename S02,
              typename S10, typename S11, typename S12,
              typename S20, typename S21, typename S22>
    finline explicit GMat( S00 m00, S01 m01, S02 m02,
                           S10 m10, S11 m11, S12 m12,
                           S20 m20, S21 m21, S22 m22 ) { MAL_STATIC_ASSERT(NC==3 && NR==3);
                                                         data[0][0] = T(m00); data[0][1] = T(m01); data[0][2] = T(m02);
                                                         data[1][0] = T(m10); data[1][1] = T(m11); data[1][2] = T(m12);
                                                         data[2][0] = T(m20); data[2][1] = T(m21); data[2][2] = T(m22); }
    // 4x4
    template <typename S00, typename S01, typename S02, typename S03,
              typename S10, typename S11, typename S12, typename S13,
              typename S20, typename S21, typename S22, typename S23,
              typename S30, typename S31, typename S32, typename S33>
    finline explicit GMat( S00 m00, S01 m01, S02 m02, S03 m03,
                           S10 m10, S11 m11, S12 m12, S13 m13,
                           S20 m20, S21 m21, S22 m22, S23 m23,
                           S30 m30, S31 m31, S32 m32, S33 m33 ) { MAL_STATIC_ASSERT(NC==4 && NR==4);
                                                                  data[0][0] = T(m00); data[0][1] = T(m01); data[0][2] = T(m02); data[0][3] = T(m03);
                                                                  data[1][0] = T(m10); data[1][1] = T(m11); data[1][2] = T(m12); data[1][3] = T(m13);
                                                                  data[2][0] = T(m20); data[2][1] = T(m21); data[2][2] = T(m22); data[2][3] = T(m23);
                                                                  data[3][0] = T(m30); data[3][1] = T(m31); data[3][2] = T(m32); data[3][3] = T(m33); }
    //@}

    //!\name Access
    //@{
    finline const T &operator()( int i, int j ) const { MAL_ASSERT( IsInRangeCO(i,0,NR) && IsInRangeCO(j,0,NC) ); return data[i][j]; }
    finline T &operator()( int i, int j ) { MAL_ASSERT( IsInRangeCO(i,0,NR) && IsInRangeCO(j,0,NC) ); return data[i][j]; }
    //@}

    //!\name Add/Sub
    //@{
    finline GMat operator+( const GMat &m ) const { GMat sum;
                                                    for(int k=0; k<cNumElements; k++)
                                                        sum.GetPtr()[k] = GetPtr()[k] + m.GetPtr()[k];
                                                    return sum; }
    finline GMat operator-( const GMat &m ) const { GMat sub;
                                                    for(int k=0; k<cNumElements; k++)
                                                        sub.GetPtr()[k] = GetPtr()[k] - m.GetPtr()[k];
                                                    return sub; }
    finline GMat operator+=( const GMat &m ) { for(int k=0; k<cNumElements; k++) GetPtr()[k] += m.GetPtr()[k]; return *this; }
    finline GMat operator-=( const GMat &m ) { for(int k=0; k<cNumElements; k++) GetPtr()[k] -= m.GetPtr()[k]; return *this; }
    //@}

    //!\name Products
    //@{
    finline GMat operator*( T val ) const { GMat prod;
                                            for(int k=0; k<cNumElements; k++)
                                                prod.GetPtr()[k] = val*GetPtr()[k];
                                            return prod; }
    finline GMat operator*=( T val ) { for(int k=0; k<cNumElements; k++) GetPtr()[k] *= val; return *this; }
    friend finline GMat operator*( T val, const GMat &m ) { return m*val; }
    //@}

    //!\name Unary
    //@{
    finline GMat operator-() const { GMat neg; for(int k=0; k<cNumElements; k++) neg.GetPtr()[k] = -(GetPtr()[k]); return neg; }
    /*\todo Consider removing, free funcs are available
    finline GMat Inverse() const { MAL_STATIC_ASSERT(NR==NC); MAL_ASSERT(false); return Zero(); }
    finline void Invert() { MAL_STATIC_ASSERT(NR==NC); MAL_ASSERT(false); }
    */
    finline GMat<T,NC,NR> Transposed() const { GMat<T,NC,NR> transp;
                                               for(int i=0; i<NC; i++)
                                                   for(int j=0; j<NR; j++)
                                                       transp(i,j) = (*this)(j,i);
                                               return transp; }
    finline void Transpose() { MAL_STATIC_ASSERT( NR == NC );
                               for(int i=0; i<NR; i++)
                                   for(int j=0; j < i; j++)
                                       std::swap( (*this)(i,j), (*this)(j,i) ); }
    //@}

    //!\name Norms
    //@{
    // Frobenius norm
    finline T NormSqF() const { T norm_sq(0); for(int k=0; k<cNumElements; k++) norm_sq += Sq(GetPtr()[k]); return norm_sq; }
    finline T NormF() const { return Sqrt( NormSqF() ); }
    //@}

    //!\name Utility
    //@{
    finline void ToArray( T *real_array ) const { for(int k=0; k<cNumElements; k++) real_array[k] = GetPtr()[k]; }
    finline void FromArray( T *real_array ) { for(int k=0; k<cNumElements; k++) GetPtr()[k] = real_array[k]; }
    //@}

public:
    //!\name Product of matrices of compatible sizes
    //@{
    template <int NR2, int NC2>
    finline GMat<T,NR,NC2> operator*( const GMat<T,NR2,NC2> &m ) const { MAL_STATIC_ASSERT( NC == NR2 );
                                                                         GMat<T,NR,NC2> prod;
                                                                         for( int i=0; i<NR; i++ )
                                                                             for( int j=0; j<NC2; j++ )
                                                                             {
                                                                                 T acc(0);
                                                                                 for( int k=0; k<NC; k++ )
                                                                                     acc += (*this)(i,k) * m(k,j);
                                                                                 prod(i,j) = acc;
                                                                             }
                                                                         return prod; }
    //! *= NOT AVAILABLE because it is ambiguous: this = m * this Vs this = this * m
    //finline void operator*=( const GMat &m );
    //@}

private:
    finline T *GetPtr() { return reinterpret_cast<T*>(data); }
    finline const T *GetPtr() const { return reinterpret_cast<const T*>(data); }

//private:
public:
    T data[NR][NC];
};

//---- Misc
template<typename T, int NR, int NC>
finline bool IsNaN( const GMat<T,NR,NC> &mat )
{
    for( int i=0; i<NR; i++ )
        for(int j=0; j<NC; j++)
            if( IsNaN( mat(i,j) ) )
                return true;
    return false;
}

template<typename T, int NR, int NC>
finline GMat<T,NC,NR> Transposed( const GMat<T,NR,NC> &m ) { return m.Transposed(); }

//! Norm
template < typename T, int NR, int NC >
finline T NormF( const GMat<T,NR,NC>& mat ) { return mat.NormF(); }

//! NormSq
template < typename T, int NR, int NC >
finline T NormSqF( const GMat<T,NR,NC>& mat ) { return mat.NormSqF(); }

//! Abs
template<typename T, int NR, int NC>
finline GMat<T,NR,NC> Abs( const GMat<T,NR,NC> &mat )
{
    GMat<T,NR,NC> res;
    for( int i=0; i<NR; i++ )
        for(int j=0; j<NC; j++)
            res(i,j) = Abs( mat(i,j) );
    return res;
}

//! Max
template<typename T, int NR, int NC>
finline T Max( const GMat<T,NR,NC> &mat )
{
    T max_entry( -Infinity<T>() );
    for( int i=0; i<NR; i++ )
        for(int j=0; j<NC; j++)
            if( max_entry < mat(i,j) )
                max_entry = mat(i,j);
    return max_entry;
}

//! NormInf
template < typename T, int NR, int NC >
finline T NormInf( const GMat<T,NR,NC>& mat ) { return Max( Abs( mat ) ); }

//! SignedNormInf
template < typename T, int NR, int NC >
finline T SignedNormInf( const GMat<T,NR,NC>& mat )
{
    T max_abs_entry( 0 );
    for( int i=0; i<NR; i++ )
        for(int j=0; j<NC; j++)
            if( mal::Abs(max_abs_entry) < mal::Abs(mat(i,j)) )
                max_abs_entry = mat(i,j);
    return max_abs_entry;
}

//---- Square-specific matrix ops
template <typename T, int N>
finline T Trace( const GMat<T,N,N> &mat )
{
    T trace(0);
    for( int i=0; i<N; i++ ) trace += mat(i,i);
    return trace;
}

template <typename T>
finline T Det( const GMat<T,2,2> &mat )
{
    return mat(0,0) * mat(1,1) - mat(1,0) * mat(0,1);
}

template <typename T>
finline T Det( const GMat<T,3,3> &mat )
{
    return mat(0,0) * mat(1,1) * mat(2,2)
        - mat(0,0) * mat(2,1) * mat(1,2)
        + mat(1,0) * mat(2,1) * mat(0,2)
        - mat(1,0) * mat(0,1) * mat(2,2)
        + mat(2,0) * mat(0,1) * mat(1,2)
        - mat(2,0) * mat(1,1) * mat(0,2);
}

// template <typename T>
// finline T Det( const GMat<T,4,4> &mat )
// {
// \todo \see Vertices3_In_Tetrahedrons3_Barycentric_Embedding.cpp for a temporary implementation, if needed globally reconsider SIMD version
// }

//----------------------------------------------------------------
// Inverse
//----------------------------------------------------------------
template <typename T> finline T Epsilon_Matrix_Inverse() { return T(0); }
template<> finline float Epsilon_Matrix_Inverse() { return 1e-6; }
template<> finline double Epsilon_Matrix_Inverse() { return 1E-12; }

template <typename T>
finline GMat<T,2,2> Inverse( const GMat<T,2,2> &mat, T determinant )
{
    GMat<T,2,2> res;
    //if( Abs(determinant) < 0.000001f ) MAL_LOG_WARNING("Inverse(m) with Det(m) = %f", determinant);
    MAL_LOG_ERROR_IF( Abs(determinant) < Epsilon_Matrix_Inverse<T>(), "Inverse(m) with Det(m) = %e", determinant );
    T inv_det( mal::Rcp(determinant) );
    res(0,0) = inv_det * mat(1,1);
    res(0,1) = -inv_det * mat(0,1);
    res(1,0) = -inv_det * mat(1,0);
    res(1,1) = inv_det * mat(0,0);
    return res;
}

template <typename T>
finline GMat<T,2,2> Adjugate( const GMat<T,2,2> &mat )
{
    GMat<T,2,2> res;
    res(0,0) = mat(1,1);
    res(0,1) = -mat(0,1);
    res(1,0) = -mat(1,0);
    res(1,1) = mat(0,0);
    return res;
}

template <typename T>
finline GMat<T,2,2> Inverse( const GMat<T,2,2> &mat )
{
    return Inverse( mat, Det(mat) );
}

template <typename T>
finline GMat<T,3,3> Inverse( const GMat<T,3,3> &mat, T determinant )
{
    GMat<T,3,3> res;
    // if( Abs(determinant) < 0.000001f ) MAL_LOG_WARNING("Inverse(m) with Det(m) = %f", determinant);
    MAL_LOG_ERROR_IF( Abs(determinant) < Epsilon_Matrix_Inverse<T>(), "Inverse(m) with Det(m) = %e", determinant );
    T inv_det( mal::Rcp(determinant) );
    res(0,0) = inv_det * ( mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1) );
    res(0,1) = inv_det * ( mat(0,2) * mat(2,1) - mat(0,1) * mat(2,2) );
    res(0,2) = inv_det * ( mat(0,1) * mat(1,2) - mat(0,2) * mat(1,1) );
    res(1,0) = inv_det * ( mat(1,2) * mat(2,0) - mat(1,0) * mat(2,2) );
    res(1,1) = inv_det * ( mat(0,0) * mat(2,2) - mat(0,2) * mat(2,0) );
    res(1,2) = inv_det * ( mat(0,2) * mat(1,0) - mat(0,0) * mat(1,2) );
    res(2,0) = inv_det * ( mat(1,0) * mat(2,1) - mat(1,1) * mat(2,0) );
    res(2,1) = inv_det * ( mat(0,1) * mat(2,0) - mat(0,0) * mat(2,1) );
    res(2,2) = inv_det * ( mat(0,0) * mat(1,1) - mat(0,1) * mat(1,0) );
    return res;
}

template <typename T>
finline GMat<T,3,3> Inverse( const GMat<T,3,3> &mat )
{
    return Inverse( mat, Det(mat) );
}

template <typename T>
finline GMat<T,3,3> Adjugate( const GMat<T,3,3> &mat )
{
    GMat<T,3,3> res;
    res(0,0) = ( mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1) );
    res(0,1) = ( mat(0,2) * mat(2,1) - mat(0,1) * mat(2,2) );
    res(0,2) = ( mat(0,1) * mat(1,2) - mat(0,2) * mat(1,1) );
    res(1,0) = ( mat(1,2) * mat(2,0) - mat(1,0) * mat(2,2) );
    res(1,1) = ( mat(0,0) * mat(2,2) - mat(0,2) * mat(2,0) );
    res(1,2) = ( mat(0,2) * mat(1,0) - mat(0,0) * mat(1,2) );
    res(2,0) = ( mat(1,0) * mat(2,1) - mat(1,1) * mat(2,0) );
    res(2,1) = ( mat(0,1) * mat(2,0) - mat(0,0) * mat(2,1) );
    res(2,2) = ( mat(0,0) * mat(1,1) - mat(0,1) * mat(1,0) );
    return res;
}

// template <typename T>
// finline GMat<T,4,4> Inverse( const GMat<T,4,4> &mat, T determinant )
// {
// }

// template <typename T>
// finline GMat<T,4,4> Inverse( const GMat<T,4,4> &mat )
// {
//     return Inverse( mat, Det(mat) );
// }

// template <typename T>
// finline GMat<T,4,4> Adjugate( const GMat<T,4,4> &mat )
// {
// }

/* 4x4 matrix inversion
   - Scalar code from http://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix (MESA-glu)
   - By now We force double computation, if ever becomes time-critical consider T-only version
   - \todo Consider SIMD version \see
     https://github.com/LiraNuna/glsl-sse2/blob/master/source/mat4.h#L324
     and the original by Intel (OLD!)
     http://download.intel.com/design/PentiumIII/sml/24504301.pdf
*/
template <typename T>
inline GMat<T,4,4> Inverse( const GMat<T,4,4>& mat )
{
    double m[16], inv[16], det;
    GMat<double,4,4> dmat( mat );
    dmat.ToArray(m);
    inv[0] = m[5]  * m[10] * m[15] -
             m[5]  * m[11] * m[14] -
             m[9]  * m[6]  * m[15] +
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] -
             m[13] * m[7]  * m[10];
    inv[4] = -m[4]  * m[10] * m[15] +
              m[4]  * m[11] * m[14] +
              m[8]  * m[6]  * m[15] -
              m[8]  * m[7]  * m[14] -
              m[12] * m[6]  * m[11] +
              m[12] * m[7]  * m[10];
    inv[8] = m[4]  * m[9] * m[15] -
             m[4]  * m[11] * m[13] -
             m[8]  * m[5] * m[15] +
             m[8]  * m[7] * m[13] +
             m[12] * m[5] * m[11] -
             m[12] * m[7] * m[9];
    inv[12] = -m[4]  * m[9] * m[14] +
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] -
               m[8]  * m[6] * m[13] -
               m[12] * m[5] * m[10] +
               m[12] * m[6] * m[9];
    inv[1] = -m[1]  * m[10] * m[15] +
              m[1]  * m[11] * m[14] +
              m[9]  * m[2] * m[15] -
              m[9]  * m[3] * m[14] -
              m[13] * m[2] * m[11] +
              m[13] * m[3] * m[10];
    inv[5] = m[0]  * m[10] * m[15] -
             m[0]  * m[11] * m[14] -
             m[8]  * m[2] * m[15] +
             m[8]  * m[3] * m[14] +
             m[12] * m[2] * m[11] -
             m[12] * m[3] * m[10];
    inv[9] = -m[0]  * m[9] * m[15] +
              m[0]  * m[11] * m[13] +
              m[8]  * m[1] * m[15] -
              m[8]  * m[3] * m[13] -
              m[12] * m[1] * m[11] +
              m[12] * m[3] * m[9];
    inv[13] = m[0]  * m[9] * m[14] -
              m[0]  * m[10] * m[13] -
              m[8]  * m[1] * m[14] +
              m[8]  * m[2] * m[13] +
              m[12] * m[1] * m[10] -
              m[12] * m[2] * m[9];
    inv[2] = m[1]  * m[6] * m[15] -
             m[1]  * m[7] * m[14] -
             m[5]  * m[2] * m[15] +
             m[5]  * m[3] * m[14] +
             m[13] * m[2] * m[7] -
             m[13] * m[3] * m[6];
    inv[6] = -m[0]  * m[6] * m[15] +
              m[0]  * m[7] * m[14] +
              m[4]  * m[2] * m[15] -
              m[4]  * m[3] * m[14] -
              m[12] * m[2] * m[7] +
              m[12] * m[3] * m[6];
    inv[10] = m[0]  * m[5] * m[15] -
              m[0]  * m[7] * m[13] -
              m[4]  * m[1] * m[15] +
              m[4]  * m[3] * m[13] +
              m[12] * m[1] * m[7] -
              m[12] * m[3] * m[5];
    inv[14] = -m[0]  * m[5] * m[14] +
               m[0]  * m[6] * m[13] +
               m[4]  * m[1] * m[14] -
               m[4]  * m[2] * m[13] -
               m[12] * m[1] * m[6] +
               m[12] * m[2] * m[5];
    inv[3] = -m[1] * m[6] * m[11] +
              m[1] * m[7] * m[10] +
              m[5] * m[2] * m[11] -
              m[5] * m[3] * m[10] -
              m[9] * m[2] * m[7] +
              m[9] * m[3] * m[6];
    inv[7] = m[0] * m[6] * m[11] -
             m[0] * m[7] * m[10] -
             m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] +
             m[8] * m[2] * m[7] -
             m[8] * m[3] * m[6];
    inv[11] = -m[0] * m[5] * m[11] +
               m[0] * m[7] * m[9] +
               m[4] * m[1] * m[11] -
               m[4] * m[3] * m[9] -
               m[8] * m[1] * m[7] +
               m[8] * m[3] * m[5];
    inv[15] = m[0] * m[5] * m[10] -
              m[0] * m[6] * m[9] -
              m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] +
              m[8] * m[1] * m[6] -
              m[8] * m[2] * m[5];
    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
    if (det == 0) return GMat<T,4,4>::Identity();
    det = 1.0 / det;
    for (int i = 0; i < 16; i++) inv[i] *= det;
    dmat.FromArray( inv );
    return GMat<T,4,4>(dmat);
}

} //namespace mal

#endif
