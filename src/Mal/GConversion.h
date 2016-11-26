/*! \file GConversion.h
\brief Generic conversions between Mal types.

\note Some UNTESTED functions fire an assert(false) to FORCE checking
their results on the first use }-D
*/

#ifndef MAL_GCONVERSION_H
#define MAL_GCONVERSION_H

#include <Mal/Config.h>
#include <Mal/RealUtils.h>
#include <Mal/GVec.h>
#include <Mal/GQuat.h>
#include <Mal/GMat.h>

namespace mal
{

// --------------------------------------------------------
// ---- GQuat conversion from several other representations
// --------------------------------------------------------

template<typename T>
finline GQuat<T> GQuat_From( const GVec<T,3> &axis, T angle )
{
    T axis_norm = axis.Norm(); //\todo INEFFICIENT if it was already normalized!
    MAL_ASSERT( axis_norm != T(0) );
    T half_angle = T(0.5) * angle;
    T tmp = Sin( half_angle ) / axis_norm;
    return GQuat<T>( tmp * axis.x(), tmp * axis.y(), tmp * axis.z(), Cos( half_angle ));
}

/*! Mat->Quat: Codi basat en ShoemakeII i comprovat a DLE.
  \todo CONVERTIR i PROVAR
  \note El codi de RotationIssues.pdf està MALAMENT!
*/
template<typename T>
finline GQuat<T> GQuat_From( const GMat<T,3,3> &mat3 )
{
    GQuat<T> q;
    T trace = mal::Trace(mat3);
    if( trace > Epsilon<T>() ) //Critical precision point
    {
        q.w = T(0.5)*Sqrt( trace + T(1) );
        T rcp_s( T(0.25) / q.w );
        q.x = ( mat3(2,1) - mat3(1,2) ) * rcp_s;
        q.y = ( mat3(0,2) - mat3(2,0) ) * rcp_s;
        q.z = ( mat3(1,0) - mat3(0,1) ) * rcp_s;
    }
    else
    {
        if ( mat3(0,0) > mat3(1,1) && mat3(0,0) > mat3(2,2) )
        {
            q.x = T(0.5)*Sqrt( T(1) + mat3(0,0) - mat3(1,1) - mat3(2,2) );
            T rcp_s( T(0.25) / q.x );
            q.y = (mat3(0,1) + mat3(1,0) ) * rcp_s;
            q.z = (mat3(2,0) + mat3(0,2) ) * rcp_s;
            q.w = (mat3(2,1) - mat3(1,2) ) * rcp_s;
        }
        else if (mat3(1,1) > mat3(2,2))
        {
            // Shoemake-tested
            q.y = T(0.5)*Sqrt( T(1) + mat3(1,1) - mat3(0,0) - mat3(2,2) );
            T rcp_s( T(0.25) / q.y );
            q.x = (mat3(0,1) + mat3(1,0) ) * rcp_s;
            q.z = (mat3(1,2) + mat3(2,1) ) * rcp_s;
            q.w = (mat3(0,2) - mat3(2,0) ) * rcp_s;
        }
        else
        {
            q.z = T(0.5)*Sqrt( T(1) + mat3(2,2) - mat3(0,0) - mat3(1,1) );
            T rcp_s( T(0.25) / q.z );
            q.x = (mat3(2,0) + mat3(0,2) ) * rcp_s;
            q.y = (mat3(1,2) + mat3(2,1) ) * rcp_s;
            q.w = (mat3(1,0) - mat3(0,1) ) * rcp_s;
        }
    }
    return q.Normalized();
}

template<typename T>
finline GQuat<T> GQuat_From( const GMat<T,2,2> &mat2 )
{
    return GQuat_From( GVec<T,3>(0,0,1), Angle_From(mat2) );
}

template<typename T>
finline GQuat<T> GQuat_From( T yaw, T pitch, T roll )
{
    MAL_ASSERT( false );
    T cos_yaw = Cos( T(0.5) * yaw );
    T sin_yaw = Sin( T(0.5) * yaw );
    T cos_pitch = Cos( T(0.5) * pitch );
    T sin_pitch = Sin( T(0.5) * pitch );
    T cos_roll = Cos( T(0.5) * roll );
    T sin_roll = Sin( T(0.5) * roll );
    return GQuat<T>( cos_roll * sin_pitch * cos_yaw + sin_roll * cos_pitch * sin_yaw,
                     cos_roll * cos_pitch * sin_yaw - sin_roll * sin_pitch * cos_yaw,
                     sin_roll * cos_pitch * cos_yaw - cos_roll * sin_pitch * sin_yaw,
                     cos_roll * cos_pitch * cos_yaw + sin_roll * sin_pitch * sin_yaw );
}

//!< Quaternion from a Vec3
template<typename T>
finline GQuat<T> GQuat_From( const GVec<T,3> &vec )
{
    return GQuat<T>( vec.x(), vec.y(), vec.z(), 0 );
}


// ----------------------------------------------------------
// ---- GMat conversion from several other representations
// ----------------------------------------------------------

//! 2D Rotation matrix on XY plane (around +Z axis)
template<typename T>
finline GMat<T,2,2> GRotation2x2_From( T angle )
{
    GMat<T,2,2> res;
    T s = Sin( angle );
    T c = Cos( angle );
    res(0,0) = c;
    res(0,1) = -s;
    res(1,0) = s;
    res(1,1) = c;
    return res;
}

template<typename T>
finline GMat<T,3,3> GRotation3x3_From( const GVec<T,3> &axis, T angle )
{
    MAL_ASSERT( axis.NormSq() != T(0) );

    GMat<T,3,3> res;

    T s = Sin( angle );
    T c = Cos( angle );
    T t = T(1) - c;

    res(0,0) = t * axis[0] * axis[0] + c;
    res(0,1) = t * axis[0] * axis[1] - s * axis[2];
    res(0,2) = t * axis[0] * axis[2] + s * axis[1];
    res(1,0) = t * axis[0] * axis[1] + s * axis[2];
    res(1,1) = t * axis[1] * axis[1] + c;
    res(1,2) = t * axis[1] * axis[2] - s * axis[0];
    res(2,0) = t * axis[0] * axis[2] - s * axis[1];
    res(2,1) = t * axis[1] * axis[2] + s * axis[0];
    res(2,2) = t * axis[2] * axis[2] + c;

    return res;
}

template<typename T>
finline GMat<T,4,4> GMat4x4_From( const GVec<T,3> &axis, T angle )
{
    MAL_ASSERT( axis.NormSq() != T(0) );

    GMat<T,4,4> res;

    T s = Sin( angle );
    T c = Cos( angle );
    T t = T(1) - c;

    res(0,0) = t * axis[0] * axis[0] + c;
    res(0,1) = t * axis[0] * axis[1] - s * axis[2];
    res(0,2) = t * axis[0] * axis[2] + s * axis[1];
    res(0,3) = 0;
    res(1,0) = t * axis[0] * axis[1] + s * axis[2];
    res(1,1) = t * axis[1] * axis[1] + c;
    res(1,2) = t * axis[1] * axis[2] - s * axis[0];
    res(1,3) = 0;
    res(2,0) = t * axis[0] * axis[2] - s * axis[1];
    res(2,1) = t * axis[1] * axis[2] + s * axis[0];
    res(2,2) = t * axis[2] * axis[2] + c;
    res(2,3) = 0;

    res(3,0) = 0;
    res(3,1) = 0;
    res(3,2) = 0;
    res(3,3) = 1;

    return res;
}

template<typename T>
finline GMat<T,3,3> GRotation3x3_From( const GQuat<T> &q )
{
    T d = q.NormSq();
    MAL_ASSERT( d != 0 );

    T s = T(2) / d;
    T xs = q.x * s,   ys = q.y * s,   zs = q.z * s;
    T wx = q.w * xs,  wy = q.w * ys,  wz = q.w * zs;
    T xx = q.x * xs,  xy = q.x * ys,  xz = q.x * zs;
    T yy = q.y * ys,  yz = q.y * zs,  zz = q.z * zs;

    GMat<T,3,3> res;
    res(0,0) = T(1) - (yy + zz);
    res(0,1) = xy - wz;
    res(0,2) = xz + wy;
    res(1,0) = xy + wz;
    res(1,1) = T(1) - (xx + zz);
    res(1,2) = yz - wx;
    res(2,0) = xz - wy;
    res(2,1) = yz + wx;
    res(2,2) = T(1) - (xx + yy);

    return res;
}

// ----------------------------------------------------------
// ---- Axis/Angle conversion from several other representations
// ----------------------------------------------------------

template<typename T>
finline GVec<T,3> Axis_From( const GQuat<T> &q )
{
    MAL_ASSERT( false );
    T tmp = T(1) - q.w*q.w;
    if( tmp > T(0) )
    {
        T inv_sqrt = T(1)/Sqrt(tmp);
        return GVec<T,3>( q.x, q.y, q.z ) * inv_sqrt;
    }
    else // tmp <= 0
        return GVec<T,3>::Zero();
}


template<typename T>
finline T Angle_From( const GQuat<T> &q )
{
    return q.GetAngle();
}

template<typename T>
finline T Angle_From( const GMat<T,2,2> &mat2 )
{
    // Angle formed by local X axis (=mat2.Column(0)) and global X axis
    return T( mal::ATan2( mat2(1,0), mat2(0,0) ) );
}

// ----------------------------------------------------------
// ---- Vector freaky ops
// ----------------------------------------------------------

//! Cast GVec<T,N> into different dimension GVec<T,D> (either dropping coords or filling with 0)
template <unsigned D, typename T, unsigned N>
finline GVec<T,D> CastDimension( const GVec<T,N> &v )
{
    GVec<T,D> res;
    if( N >= D )
        for(unsigned int i=0; i<D; i++) res[i] = v[i];
    else
    {
        for(unsigned int i=0; i<N; i++) res[i] = v[i];
        for(unsigned int i=N; i<D; i++) res[i] = T(0);
    }
    return res;
}

} // namespace mal

#endif // MAL_GCONVERSION_H
