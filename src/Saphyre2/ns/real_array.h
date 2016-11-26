#ifndef S2_NS_REAL_ARRAY_H
#define S2_NS_REAL_ARRAY_H

#include <Saphyre2/ns/Config.h>

#include <iostream> //req by print()

#ifdef __S2_NS_ENABLE_SIMD
#  include <smmintrin.h>
#endif

namespace S2 { namespace ns {

//@{
namespace real_array
{

//---- Non-vectorizable ops
// IsNan( v[i] )
finline bool is_nan( const Real* v, unsigned int n ) { for( unsigned int i=0; i<n; i++ ) if( mal::IsNaN(v[i]) ) return true; return false; }
// print
finline void print( const char* str, const Real* v, unsigned int n ) { std::cout << str << " = "; for( unsigned int i=0; i<n; i++ ) std::cout << v[i] << ","; std::cout << std::endl; }

#ifdef __S2_NS_ENABLE_SIMD

finline Real* alloc( unsigned int n ) { size_t an(4*((n+3)/4)); return reinterpret_cast<Real*>(_mm_malloc( an*sizeof(Real), 16 )); }
finline void dealloc( Real* v ) { _mm_free( v ); }

// v = 0
finline void zero( Real* v, unsigned int n ) { size_t an(4*((n+3)/4)); memset( v, 0, an*sizeof(Real) ); }
// v0 = v1
finline void assign( Real*__restrict__ v0, const Real*__restrict__ v1, unsigned int n ) { size_t an(4*((n+3)/4)); memcpy( v0, v1, an*sizeof(Real) ); }

finline Real dot( const Real*__restrict__ v0, const Real*__restrict__ v1, unsigned int n ) //{ Real res(0); for( unsigned int i=0; i<n; i++ ) res += v0[i]*v1[i]; return res; }
{
    __m128 acc(_mm_setzero_ps());
    size_t an(4*((n+3)/4));
    for( unsigned int i=0; i<an; i+=4 )
        acc = _mm_add_ps( acc, _mm_dp_ps( _mm_load_ps(&v0[i]), _mm_load_ps(&v1[i]), 0xF1) );
    return _mm_cvtss_f32(acc);
}
finline Real sq_norm_2( const Real* v, unsigned int n ) //{ Real res(0); for( unsigned int i=0; i<n; i++ ) res += mal::Sq(v[i]); return res; }
{
    __m128 acc(_mm_setzero_ps());
    size_t an(4*((n+3)/4));
    for( unsigned int i=0; i<an; i+=4 )
        acc = _mm_add_ps( acc, _mm_dp_ps( _mm_load_ps(&v[i]), _mm_load_ps(&v[i]), 0xF1) );
    return _mm_cvtss_f32(acc);
}
finline Real norm_2( const Real* v, unsigned int n ) { return mal::Sqrt( sq_norm_2(v,n) ); }

// v0 = s*v1
finline void assign_scaled( Real*__restrict__ v0, Real s, const Real*__restrict__ v1, unsigned int n ) //{ for( unsigned int i=0; i<n; i++ ) v0[i] = s*v1[i]; }
{ __m128 ss(_mm_load1_ps(&s)); size_t an(4*((n+3)/4)); for( unsigned int i=0; i<an; i+=4 ) _mm_store_ps( &v0[i], _mm_mul_ps( _mm_load_ps(&v1[i]), ss ) ); }
// v0 += v1
finline void addeq( Real*__restrict__ v0, const Real*__restrict__ v1, unsigned int n ) //{ for( unsigned int i=0; i<n; i++ ) v0[i] += v1[i]; }
{ size_t an(4*((n+3)/4)); for( unsigned int i=0; i<an; i+=4 ) _mm_store_ps( &v0[i], _mm_add_ps( _mm_load_ps(&v0[i]), _mm_load_ps(&v1[i]) ) ); }
// v0 -= v1
finline void subeq( Real*__restrict__ v0, const Real*__restrict__ v1, unsigned int n ) //{ for( unsigned int i=0; i<n; i++ ) v0[i] -= v1[i]; }
{ size_t an(4*((n+3)/4)); for( unsigned int i=0; i<an; i+=4 ) _mm_store_ps( &v0[i], _mm_sub_ps( _mm_load_ps(&v0[i]), _mm_load_ps(&v1[i]) ) ); }
// v *= s
finline void muleq( Real* v, Real s, unsigned int n ) //{ for( unsigned int i=0; i<n; i++ ) v0[i] *= v1[i]; }
{ __m128 ss(_mm_load1_ps(&s)); size_t an(4*((n+3)/4)); for( unsigned int i=0; i<an; i+=4 ) _mm_store_ps( &v[i], _mm_mul_ps( _mm_load_ps(&v[i]), ss ) ); }

// v0 += s*v1
finline void addeq_scaled( Real*__restrict__ v0, Real s, const Real*__restrict__ v1, unsigned int n ) //{ for( unsigned int i=0; i<n; i++ ) v0[i] += s*v1[i]; }
{ __m128 ss(_mm_load1_ps(&s)); size_t an(4*((n+3)/4)); for( unsigned int i=0; i<an; i+=4 ) _mm_store_ps( &v0[i], _mm_add_ps( _mm_load_ps(&v0[i]), _mm_mul_ps( _mm_load_ps(&v1[i]), ss ) ) ); }
// v0 -= s*v1
finline void subeq_scaled( Real*__restrict__ v0, Real s, const Real*__restrict__ v1, unsigned int n ) //{ for( unsigned int i=0; i<n; i++ ) v0[i] -= s*v1[i]; }
{ __m128 ss(_mm_load1_ps(&s)); size_t an(4*((n+3)/4)); for( unsigned int i=0; i<an; i+=4 ) _mm_store_ps( &v0[i], _mm_sub_ps( _mm_load_ps(&v0[i]), _mm_mul_ps( _mm_load_ps(&v1[i]), ss ) ) ); }

// v0 = s*v0 + v1
/*
finline void muleq_addeq( Real*__restrict__ v0, Real s, const Real*__restrict__ v1, unsigned int n ) //{ for( unsigned int i=0; i<n; i++ ) v0[i] = s*v0[i] + v1[i]; }
{ __m128 ss(_mm_load1_ps(&s)); size_t an(4*((n+3)/4)); for( unsigned int i=0; i<an; i+=4 ) _mm_store_ps( &v0[i], _mm_add_ps( _mm_mul_ps( _mm_load_ps(&v0[i]), ss ), _mm_load_ps(&v1[i]) ) ); }
*/
//TEMP: BUG caused by previous SSE function yields inverted elements in armadillo1K, SEEMS correct, may be something unrelated that arises due to different results in SSE computations, using scalar code by now.
finline void muleq_addeq( Real*__restrict__ v0, Real s, const Real*__restrict__ v1, unsigned int n ) { for( unsigned int i=0; i<n; i++ ) v0[i] = s*v0[i] + v1[i]; }

// TEMPORAL: These should be elsewhere...
template <typename T, int D>
finline mal::GVec<T,D>* alloc_GVec( unsigned int n ) { size_t nf(n*D); size_t anf(4*((nf+3)/4)); return reinterpret_cast<mal::GVec<T,D>*>(_mm_malloc( anf*sizeof(Real), 16 )); }
template <typename T, int D>
finline void dealloc_GVec( mal::GVec<T,D>* v ) { _mm_free( v ); }

#else

finline Real* alloc( unsigned int n ) { return reinterpret_cast<Real*>(malloc( n*sizeof(Real) )); }
finline void dealloc( Real* v ) { free( v ); }

// v = 0
finline void zero( Real* v, unsigned int n ) { memset( v, 0, n*sizeof(Real) ); }
// v0 = v1
finline void assign( Real* v0, const Real* v1, unsigned int n ) { memcpy( v0, v1, n*sizeof(Real) ); }

finline Real dot( const Real* v0, const Real* v1, unsigned int n ) { Real res(0); for( unsigned int i=0; i<n; i++ ) res += v0[i]*v1[i]; return res; }
finline Real sq_norm_2( const Real* v, unsigned int n ) { Real res(0); for( unsigned int i=0; i<n; i++ ) res += mal::Sq(v[i]); return res; }
finline Real norm_2( const Real* v, unsigned int n ) { return mal::Sqrt( sq_norm_2(v,n) ); }

// v0 = s*v1
finline void assign_scaled( Real* v0, Real s, const Real* v1, unsigned int n ) { for( unsigned int i=0; i<n; i++ ) v0[i] = s*v1[i]; }
// v0 += v1
finline void addeq( Real* v0, const Real* v1, unsigned int n ) { for( unsigned int i=0; i<n; i++ ) v0[i] += v1[i]; }
// v0 -= v1
finline void subeq( Real* v0, const Real* v1, unsigned int n ) { for( unsigned int i=0; i<n; i++ ) v0[i] -= v1[i]; }
// v *= s
finline void muleq( Real* v, Real s, unsigned int n ) { for( unsigned int i=0; i<n; i++ ) v[i] *= s; }
// v0 += s*v1
finline void addeq_scaled( Real* v0, Real s, const Real* v1, unsigned int n ) { for( unsigned int i=0; i<n; i++ ) v0[i] += s*v1[i]; }
// v0 -= s*v1
finline void subeq_scaled( Real* v0, Real s, const Real* v1, unsigned int n ) { addeq_scaled( v0, -s, v1, n ); }
// v0 = s*v0 + v1
finline void muleq_addeq( Real* v0, Real s, const Real* v1, unsigned int n ) { for( unsigned int i=0; i<n; i++ ) v0[i] = s*v0[i] + v1[i]; }

// TEMPORAL: These should be elsewhere...
template <typename T, int D>
finline mal::GVec<T,D>* alloc_GVec( unsigned int n ) { return reinterpret_cast<mal::GVec<T,D>*>(malloc( n*sizeof(mal::GVec<T,D>) )); }
template <typename T, int D>
finline void dealloc_GVec( mal::GVec<T,D>* v ) { free( v ); }

#endif

}
//@}

}} //namespace S2::ns

#endif //S2_NS_REAL_ARRAY_H
