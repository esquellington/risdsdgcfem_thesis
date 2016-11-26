#ifndef MAL_GQUAT_UTILS_H
#define MAL_GQUAT_UTILS_H

#include <Mal/GQuat.h>
#include <Mal/GVec.h>

namespace mal
{

// Quat-Vector product
template <typename T>
inline GVec<T,3> operator*( const GQuat<T> &q, const GVec<T,3> &v )
{
    GQuat<T> lTmpQuat( q * GQuat<T>(v.x(),v.y(),v.z(),0) * q.Conjugated() );
    return GVec<T,3>( lTmpQuat.x, lTmpQuat.y, lTmpQuat.z );
}

} // namespace mal

#endif // MAL_GQUAT_UTILS_H
