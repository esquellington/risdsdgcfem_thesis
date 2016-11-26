#ifndef MAL_GQUAT_H
#define MAL_GQUAT_H

#include <Mal/Config.h>
#include <Mal/RealUtils.h>

namespace mal
{

template<typename T>
class GQuat
{
public:
    typedef T real_type;
    static const unsigned size_in_reals = 4;
    finline static GQuat Identity() { return GQuat(T(0),T(0),T(0),T(1)); }

public:    
    T x,y,z,w; //!< Components
    
public:

    //!\name Construction
    //@{
    finline GQuat() {}
    finline GQuat( T qx, T qy, T qz, T qw ) { x = qx; y = qy; z = qz; w = qw; }
    finline GQuat( const GQuat &q ) { x = q.x; y = q.y; z = q.z; w = q.w; }
    finline GQuat( const T *real_array ) { x = real_array[0]; y = real_array[1]; z = real_array[2]; w = real_array[3]; }
    //}@}

    //!\name Add/Sub
    //@{
    finline GQuat operator+( const GQuat &q ) const { return GQuat( x + q.x, y + q.y, z + q.z, w + q.w ); }
    finline GQuat operator-( const GQuat &q ) const { return GQuat( x - q.x, y - q.y, z - q.z, w - q.w ); }
    finline GQuat operator+=( const GQuat &q ) { x += q.x; y += q.y; z += q.z; w += q.w; return *this; }
    finline GQuat operator-=( const GQuat &q ) { x -= q.x; y -= q.y; z -= q.z; w -= q.w; return *this; }
    //@}

    //!\name Products
    //@{
    finline GQuat operator*( T scalar ) const { return GQuat( x*scalar, y*scalar, z*scalar, w*scalar ); }
    finline void operator*=( T scalar ) { x *= scalar; y *= scalar; z *= scalar; w *= scalar; }
    friend finline GQuat operator*( T val, const GQuat &q ) { return q*val; }
    
    finline GQuat operator/( T scalar ) const { MAL_ASSERT(scalar != 0); T inv(T(1)/scalar); return ( (*this) * inv ); }
    finline void operator/=( T scalar ) { MAL_ASSERT(scalar != 0); T inv(T(1)/scalar); (*this) *= inv; }
    
    GQuat operator*( const GQuat &q ) const { return GQuat( w*q.x + x*q.w + y*q.z - z*q.y,
                                                            w*q.y + y*q.w + z*q.x - x*q.z,
                                                            w*q.z + z*q.w + x*q.y - y*q.x,
                                                            w*q.w - x*q.x - y*q.y - z*q.z ); }
    //! *= NOT AVAILABLE because it is ambiguous: this = q * this Vs this = this * q
    //finline void operator*=( const GQuat &q );
    //@}
    
    //!\name Unary
    //@{
    finline GQuat operator-() const { return GQuat(-x,-y,-z,-w); }
    finline GQuat Normalized() const { return (*this / Norm()); }
    finline GQuat Conjugated() const { return GQuat( -x, -y, -z, w ); }
    finline GQuat Inverse() const { return Conjugated(); }
    
    finline void Normalize() { *this /= Norm(); }
    finline void Conjugate() { x = -x; y = -y; z = -z; /*w = w;*/ }
    finline void Invert()  { Conjugate(); }
    //@}    

    //!\name Norms
    //@{
    finline T Dot( const GQuat &q ) const { return ( x*q.x + y*q.y + z*q.z + w*q.w ); }
    finline T NormSq() const { return this->Dot(*this); }
    finline T Norm() const { return Sqrt( NormSq() ); }
    //@}
    
    //!\name Utility
    //@{
    finline void ToArray( T *real_array ) const { real_array[0] = x; real_array[1] = y; real_array[2] = z; real_array[3] = w; }
    finline void FromArray( const T *real_array ) { x = real_array[0]; y = real_array[1]; z = real_array[2]; w = real_array[3]; }

    finline T GetAngle() const { MAL_ASSERT(false); return T(2)*ACos(w); } // WARNIIING, may return angles > 180!! check sign!
    finline GQuat Lerp( const GQuat &q, T lambda ) const { if( (*this).Dot(q) > T(0) )
                                                               return (lambda*(*this) + (T(1)-lambda)*q);
                                                           else
                                                               return (lambda*(*this) - (T(1)-lambda)*q); }
    //\todo NLerp, SLerp
    //@}
};

//---- Misc
template< typename T >
finline bool IsNaN( const GQuat<T> &q )
{
    return IsNaN(q.x) || IsNaN(q.y) || IsNaN(q.z) || IsNaN(q.w);
}

} // namespace mal

#endif // MAL_GQUAT_H
