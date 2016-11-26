#ifndef GEO_NP_SIMD_H
#define GEO_NP_SIMD_H

#include <Geo/Config.h>
// #include "Context.h"
// #include "stats.h"

#ifdef __GEO_ENABLE_NP_SIMD
#  include "smmintrin.h"

namespace simd
{

/* Simple SSE4 Vec3f class
   \sa from http://fastcpp.blogspot.com.es/2011/12/simple-V3f-class-with-sse-support.html

   \sa SSE intrinsics guide https://software.intel.com/sites/landingpage/IntrinsicsGuide
*/
class V3f
{
public:
    union
    {
        struct { float x, y, z, w; }; //w must be 0 (vector) or 1 (point)
        __m128 _data;
    };
public:

    // Constructors
    // inline V3f() : _data(_mm_setzero_ps()) {} 0-init
    inline V3f() {} //TEMP testing empty constructor as in GVec.h, SHOULD WORK, WE MUST NOT ASSUME 0-init
    inline V3f( const V3f& v ) : _data(v._data) {}
    inline V3f( float x, float y, float z ) : _data(_mm_set_ps( 0, z, y, x )) {}
    inline V3f( float x, float y, float z, float w ) : _data(_mm_set_ps( w, z, y, x )) {}
    inline V3f( __m128 m ) : _data(m) {}
    inline explicit V3f( const Vec3f& v ) : _data(_mm_set_ps( 0, v.z(), v.y(), v.x() )) {}
    inline explicit V3f( const Vec3f& v, float w ) : _data(_mm_set_ps( w, v.z(), v.y(), v.x() )) {}

    // Casts
    inline Vec3f AsVec3f() const { return Vec3f(x,y,z); }

    // Arithmetic V3f OP V3f
    inline V3f operator+( const V3f& b ) const { return _mm_add_ps(_data, b._data); }
    inline V3f operator-( const V3f& b ) const { return _mm_sub_ps(_data, b._data); }
    // inline V3f operator*( const V3f& b ) const { return _mm_mul_ps(_data, b._data); } \todo ambiguous OP... not necessary
    // inline V3f operator/( const V3f& b ) const { return _mm_div_ps(_data, b._data); } \todo ambiguous OP... not necessary //\todo w div-by-0??

    // Arithmetic V3f OP= V3f
    inline V3f& operator+=( const V3f& b ) { _data = _mm_add_ps(_data, b._data); return *this; }
    inline V3f& operator-=( const V3f& b ) { _data = _mm_sub_ps(_data, b._data); return *this; }
    // inline V3f& operator*=( const V3f& b ) { _data = _mm_mul_ps(_data, b._data); return *this; } \todo ambiguous OP... not necessary
    // inline V3f& operator/=( const V3f& b ) { _data = _mm_div_ps(_data, b._data); return *this; } \todo ambiguous OP... not necessary //\todo w div-by-0??

    // Arithmetic V3f OP float
    // inline V3f operator+( float b ) const { return _mm_add_ps(_data, _mm_set1_ps(b)); } \todo ambiguous OP... not necessary
    // inline V3f operator-( float b ) const { return _mm_sub_ps(_data, _mm_set1_ps(b)); } \todo ambiguous OP... not necessary
    inline V3f operator*( float b ) const { return _mm_mul_ps(_data, _mm_set1_ps(b)); }
    inline V3f operator/( float b ) const { return _mm_div_ps(_data, _mm_set1_ps(b)); }

    // Arithmetic V3f OP= float
    // inline V3f& operator+=( float b ) { _data = _mm_add_ps(_data, _mm_set1_ps(b)); return *this; } \todo ambiguous OP... not necessary
    // inline V3f& operator-=( float b ) { _data = _mm_sub_ps(_data, _mm_set1_ps(b)); return *this; } \todo ambiguous OP... not necessary
    inline V3f& operator*=( float b ) { _data = _mm_mul_ps(_data, _mm_set1_ps(b)); return *this; }
    inline V3f& operator/=( float b ) { _data = _mm_div_ps(_data, _mm_set1_ps(b)); return *this; }

    // Cross
    inline V3f Cross( const V3f& b ) const
        {
            return _mm_sub_ps( _mm_mul_ps(_mm_shuffle_ps(_data, _data, _MM_SHUFFLE(3, 0, 2, 1)), _mm_shuffle_ps(b._data, b._data, _MM_SHUFFLE(3, 1, 0, 2))),
                               _mm_mul_ps(_mm_shuffle_ps(_data, _data, _MM_SHUFFLE(3, 1, 0, 2)), _mm_shuffle_ps(b._data, b._data, _MM_SHUFFLE(3, 0, 2, 1))) );
        }

    // Dot \note requires SSE4.1
    inline float Dot( const V3f& b ) const { return _mm_cvtss_f32(_mm_dp_ps(_data, b._data, 0x71)); } // +(xx,yy,zz)
    inline float Dot4( const V3f& b ) const { return _mm_cvtss_f32(_mm_dp_ps(_data, b._data, 0xF1)); } // +(xx,yy,zz,ww)

    // Neg
    inline V3f operator-() const { return _mm_sub_ps(_mm_set1_ps(0), _data); }

    /*\todo
    // length of the vector
    inline float length() const { return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(_data, _data, 0x71))); }
    // 1/length() of the vector
    inline float rlength() const { return _mm_cvtss_f32(_mm_rsqrt_ss(_mm_dp_ps(_data, _data, 0x71))); }
    // returns the vector scaled to unit length
    inline V3f normalize() const { return _mm_mul_ps(_data, _mm_rsqrt_ps(_mm_dp_ps(_data, _data, 0x7F))); }
    */

    /*\todo
    // overloaded operators that ensure alignment
    inline void* operator new[](size_t x) { return _aligned_malloc(x, 16); }
    inline void operator delete[](void* x) { if (x) _aligned_free(x); }
    */
} __attribute__ ((aligned (16)));

// Free ops
inline V3f operator*( float a, const V3f& b ) { return b * a; }
// inline V3f operator+( float a, const V3f& b ) { return b + a; } \todo ambiguous OP... not necessary
// inline V3f operator-( float a, const V3f& b ) { return V3f(_mm_set1_ps(a)) - b; } \todo ambiguous OP... not necessary
// inline V3f operator/( float a, const V3f& b ) { return V3f(_mm_set1_ps(a)) / b; } \todo ambiguous OP... not necessary
inline float Dot( const V3f& a, const V3f& b ) { return a.Dot(b); }
inline float Dot4( const V3f& a, const V3f& b ) { return a.Dot4(b); }
inline V3f Cross( const V3f& a, const V3f& b ) { return a.Cross(b); }

// \todo batch-transform all
// inline void transform_array_V3f( const Transform3f& m, V3f* vec_v3f )
// {
// }

//\todo simd::M34f to TRANSFORM N>>1 vertices simd::V3f
// class M34f
// {
// public:
// __m128 m_Row0;
// __m128 m_Row1;
// __m128 m_Row2;
// public:
// inline M34f( const Mat4x4f& m ) {}
// inline V3f operator*( const V3f& v ) {}
// };

} //namespace simd

#  endif //__GEO_ENABLE_NP_SIMD

#endif //GEO_NP_SIMD_H
