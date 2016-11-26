#ifndef MAL_GMAT_UTILS_H
#define MAL_GMAT_UTILS_H

#include <Mal/Config.h>
#include <Mal/GMat.h>
#include <Mal/GVec.h>

namespace mal
{

//! Mat-Vector product
template <typename T, int NR, int NC>
inline GVec<T,NR> operator*( const GMat<T,NR,NC> &m, const GVec<T,NC> &v )
{
    GVec<T,NR> res;
    for(int i=0; i<NR; i++)
    {
        T acc(0);
        for(int j=0; j<NC; j++)
            acc += m(i,j) * v[j];
        res[i] = acc;
    }
    return res;
}

//! Mat row vector extraction
template <typename T, int NR, int NC>
inline GVec<T,NC> Row( int row, const GMat<T,NR,NC> &m )
{
    MAL_ASSERT( row < NR );
    GVec<T,NC> res;
    for(int j=0; j<NC; j++) res[j] = m(row,j);
    return res;
}

//! Mat column vector extraction
template <typename T, int NR, int NC>
inline GVec<T,NR> Column( int col, const GMat<T,NR,NC> &m )
{
    MAL_ASSERT( col < NC );
    GVec<T,NR> res;
    for(int i=0; i<NR; i++) res[i] = m(i,col);
    return res;
}

//! Mat row vector extraction
template <unsigned I, typename T, int NR, int NC>
inline GVec<T,NC> GRow( const GMat<T,NR,NC> &m )
{
    MAL_STATIC_ASSERT( I < NR );
    GVec<T,NC> res;
    for(int j=0; j<NC; j++) res[j] = m(I,j);
    return res;
}

//! Mat column vector extraction
template <unsigned J, typename T, int NR, int NC>
inline GVec<T,NR> GColumn( const GMat<T,NR,NC> &m )
{
    MAL_STATIC_ASSERT( J < NC );
    GVec<T,NR> res;
    for(int i=0; i<NR; i++) res[i] = m(i,J);
    return res;
}

//! Mat row vector write
template <unsigned I, typename T, int NR, int NC>
inline void GSetRow( GMat<T,NR,NC> &m, const GVec<T,NC> &v )
{
    MAL_STATIC_ASSERT( I < NR );
    for(int j=0; j<NC; j++) m(I,j) = v[j];
}

//! Mat column vector write
template <unsigned J, typename T, int NR, int NC>
inline void GSetColumn( GMat<T,NR,NC> &m, const GVec<T,NR> &v )
{
    MAL_STATIC_ASSERT( J < NC );
    for(int i=0; i<NR; i++) m(i,J) = v[i];
}

//! Mat row vector write
template <typename T, int NR, int NC>
inline void SetRow( GMat<T,NR,NC> &m, int row, const GVec<T,NC> &v )
{
    for(int j=0; j<NC; j++) m(row,j) = v[j];
}

//! Mat column vector write
template <typename T, int NR, int NC>
inline void SetColumn( GMat<T,NR,NC> &m, int col, const GVec<T,NR> &v )
{
    for(int i=0; i<NR; i++) m(i,col) = v[i];
}

template < int I0, int J0, int I1, int J1,
           typename T, int NR, int NC >
inline GMat<T,I1-I0+1,J1-J0+1> GRange( const GMat<T,NR,NC> &m )
{
    MAL_STATIC_ASSERT( I0 >= 0 && I0 < I1 && I1 < NR && J0 >= 0 && J0 < J1 && J1 < NC );
    GMat<T,I1-I0+1,J1-J0+1> r;
    for(int i=I0; i<=I1; i++)
        for(int j=J0; j<=J1; j++)
            r(i-I0,j-J0) = m(i,j);
    return r;
}

template < int I0, int J0, int I1, int J1,
           typename T, int NR, int NC >
inline void GSetRange( GMat<T,NR,NC> &m, const GMat<T,I1-I0+1,J1-J0+1> &rm )
{
    MAL_STATIC_ASSERT( I0 >= 0 && I0 < I1 && I1 < NR && J0 >= 0 && J0 < J1 && J1 < NC );
    for(int i=I0; i<=I1; i++)
        for(int j=J0; j<=J1; j++)
            m(i,j) = rm(i-I0,j-J0);
}

template <typename T>
inline GMat<T,2,2> GMat2x2_From_Columns( const GVec<T,2> &e0, const GVec<T,2> &e1 )
{
    return GMat<T,2,2>( e0[0], e1[0],
                        e0[1], e1[1] );
}

template <typename T>
inline GMat<T,2,2> GMat2x2_From_Rows( const GVec<T,2> &r0, const GVec<T,2> &r1 )
{
    return GMat<T,2,2>( r0[0], r0[1],
                        r1[0], r1[1] );
}

template <typename T>
inline GMat<T,3,3> GMat3x3_From_Columns( const GVec<T,3> &e0, const GVec<T,3> &e1, const GVec<T,3> &e2 )
{
    return GMat<T,3,3>( e0[0], e1[0], e2[0],
                        e0[1], e1[1], e2[1],
                        e0[2], e1[2], e2[2] );
}

template <typename T>
inline GMat<T,3,3> GMat3x3_From_Rows( const GVec<T,3> &e0, const GVec<T,3> &e1, const GVec<T,3> &e2 )
{
    return GMat<T,3,3>( e0[0], e0[1], e0[2],
                        e1[0], e1[1], e1[2],
                        e2[0], e2[1], e2[2] );
}

template <typename T>
inline GMat<T,3,3> GMat3x3_From_Skew( const GVec<T,3> &v )
{
    return GMat<T,3,3>(     0, -v[2],  v[1],
                         v[2],     0, -v[0],
                        -v[1],  v[0],     0  );
}

template <typename T, int N>
inline GMat<T,N,N> GMatNxN_From_Diagonal( const GVec<T,N> &diag )
{
    GMat<T,N,N> res(T(0));
    for( int i=0; i<N; i++ ) res(i,i) = diag[i];
    return res;
}

/*! Build a reflection matrix from a given (unitary) axis
  \note Its based on:
  R' = Reflection( Mat2x2 R, Vec2 a )
  {
    Vec2 e0( R(0,0),
             R(1,0) );
    Vec2 e1( R(0,1),
             R(1,1) );
    Vec2 re0( e0 - (2 * mal::Dot(e0,a)) * a );
    Vec2 re1( e1 - (2 * mal::Dot(e1,a)) * a );
    R'(0,0) = re0[0]; R'(0,1) = re1[0];
    R'(1,0) = re0[1]; R'(1,1) = re1[1];
    return R';
  }
  Which explicitly *reflects* each R column vector against a.
  This linear op can be encoded into a reflection matrix A, so that:
    A = GReflection2x2_From_Axis(a);
    R' = A * R;
*/
template <typename T>
inline GMat<T,2,2> GReflection2x2_From_Axis( const GVec<T,2> &axis )
{
    T minus_two_axis_x_axis_y( - T(2)*axis[0]*axis[1] );
    return GMat<T,2,2>(  T(1) - T(2)*mal::Sq(axis[0]) ,      minus_two_axis_x_axis_y ,
                         minus_two_axis_x_axis_y      , T(1) - T(2)*mal::Sq(axis[1]) );
}

/*! Build a reflection matrix from a given (unitary) axis
  \note Its based on:
  R' = Reflection( Mat3x3 R, Vec3 a )
  {
    Vec3 e0( R(0,0),
             R(1,0),
             R(2,0) );
    Vec3 e1( R(0,1),
             R(1,1),
             R(2,1));
    Vec3 e2( R(0,2),
             R(1,2),
             R(2,2) );
    Vec3 re0( e0 - (2 * mal::Dot(e0,a)) * a );
    Vec3 re1( e1 - (2 * mal::Dot(e1,a)) * a );
    Vec3 re2( e2 - (2 * mal::Dot(e2,a)) * a );
    R'(0,0) = re0[0]; R'(0,1) = re1[0]; R'(0,2) = re2[0];
    R'(1,0) = re0[1]; R'(1,1) = re1[1]; R'(1,2) = re2[1];
    R'(2,0) = re0[2]; R'(2,1) = re1[2]; R'(2,2) = re2[2];
    return R';
  }
  Which explicitly *reflects* each R column vector against a.
  This linear op can be encoded into a reflection matrix A, so that:
    A = GReflection3x3_From_Axis(a);
    R' = A * R;
*/
template <typename T>
inline GMat<T,3,3> GReflection3x3_From_Axis( const GVec<T,3> &axis )
{
    T minus_two_axis_x_axis_y( - T(2)*axis[0]*axis[1] );
    T minus_two_axis_x_axis_z( - T(2)*axis[0]*axis[2] );
    T minus_two_axis_y_axis_z( - T(2)*axis[1]*axis[2] );
    return GMat<T,3,3>( T(1) - T(2)*Sq(axis[0]) , minus_two_axis_x_axis_y , minus_two_axis_x_axis_z ,
                        minus_two_axis_x_axis_y , T(1) - T(2)*Sq(axis[1]) , minus_two_axis_y_axis_z ,
                        minus_two_axis_x_axis_z , minus_two_axis_y_axis_z , T(1) - T(2)*Sq(axis[2]) );
}

/* Build a rotation matrix that rotates v0 into v1
  \note Correct, proved on paper that v1 = R * v0
*/
template <typename T>
inline GMat<T,2,2> GRotation2x2_Vec2Vec( const GVec<T,2> &v0, const GVec<T,2> &v1 )
{
    T cos_theta( mal::Dot(v0,v1) );
    T sin_theta( mal::Dot(v0, GVec<T,2>(v1[1],-v1[0]) ) ); //(v1[1],-v1[0]) == -mal::Perpendicular(v1), but with explicit orientation
    return GMat<T,2,2>( cos_theta, -sin_theta,
                        sin_theta, cos_theta );
}

/* Relative error between two matrices, defined as the largest
   relative error of any entry.
*/
template < typename T, int NR, int NC >
inline T GRelativeError( const GMat<T,NR,NC>& m1, const GMat<T,NR,NC>& m2 )
{
    T max_rel_error(0);
    for(int i=0; i<NC; i++)
        for(int j=0; j<NR; j++)
        {
            T error = Abs( m1(i,j) - m2(i,j) );
            T rel_error = ( Abs(m1(i,j)) > Epsilon<T>() )
                          ? error / Abs(m1(i,j))
                          : ( Abs(m2(i,j)) > Epsilon<T>() )
                            ? error / Abs(m2(i,j))
                            : 0;
            if( rel_error > max_rel_error ) max_rel_error = rel_error;
        }
    return max_rel_error;
}

template <typename T, int N>
inline GMat<T,N,N> GMatNxN_Symmetric_Part( const GMat<T,N,N>& m )
{
    return T(0.5) * (m + Transposed(m));
}

} // namespace mal

#endif // MAL_GMAT_UTILS_H
