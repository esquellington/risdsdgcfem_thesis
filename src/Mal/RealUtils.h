#ifndef MAL_REAL_UTILS_H
#define MAL_REAL_UTILS_H

#include <Mal/Config.h>
#include <Mal/NumberUtils.h>
#include <math.h>
#include <float.h>
#include <limits.h>

namespace mal
{

//---- Fundamental Functions
template<typename T>
inline T Floor( T value ) { return T( floor( float(value) ) ); }

template<typename T>
inline T Ceil( T value ) { return T( ceil( float(value) ) ); }

template<typename T>
inline T Fraction( T value ) { double intpart; return T( modf(value,&intpart) ); }

template<typename IntegralT, typename T>
inline IntegralT IntPart( T value ) { double intpart; modf(value,&intpart); return IntegralT(intpart); }

template<typename T>
inline T Sqrt( T value ) { return T( sqrt( float(value) ) ); }

template<typename T>
inline T Sq( T value ) { return T( value*value ); }

template<typename T>
inline T Cube( T value ) { return T( value*value*value ); }

template<typename T> inline T Rcp( T value ) { MAL_ASSERT(false); } //\note You do NOT want to call Rcp non-floating-point types
template<> inline float Rcp( float value ) { return float(1)/value; }
template<> inline double Rcp( double value ) { return double(1)/value; }

template<typename T> inline T Abs( T value ) { return (value>T(0)) ? value : -value; } //\todo SLOW
template<> inline float Abs( float value ) { return fabs(value); }
template<> inline double Abs( double value ) { return fabs(value); }

//---- Trigonometric Functions
template<typename T>
inline T Sin( T value ) { return T( sin( float(value) ) ); }

template<typename T>
inline T Cos( T value ) { return T( cos( float(value) ) ); }

template<typename T>
inline T Tan( T value ) { return T( tan( float(value) ) ); }

template<typename T>
inline T ASin( T value ) { return T( asin( float(value) ) ); }

template<typename T>
inline T ACos( T value ) { return T( acos( float(value) ) ); }

template<typename T>
inline T ATan( T value ) { return T( atan( float(value) ) ); }

template<typename T>
inline T ATan2( T y, T x ) { return T( atan2( float(y), float(x) ) ); }

template<typename T>
inline T Log( T value ) { return T( log( float(value) ) ); }

template<typename T>
inline T Exp( T value ) { return T( exp( float(value) ) ); }

template<typename T>
inline T Log10( T value ) { return T( log10( float(value) ) ); }

template<typename T>
inline T Exp10( T value ) { return T( pow( T(10), float(value) ) ); }

template<typename T>
inline T Pow( T base, T exponent ) { return T( pow( base, exponent ) ); }

//\note we DO NOT provide generic Inf/-Inf implementation, must be specialized for ANY that uses it
template<typename T> inline T Infinity(); //{ return T(3.4e+38f); }
template<typename T> inline T MinusInfinity(); //{ return -T(3.4e+38f); }

template<> inline unsigned int Infinity() { return UINT_MAX; }
template<> inline int Infinity() { return INT_MAX; }
template<> inline int MinusInfinity() { return -INT_MAX; } //\note abs(INT_MIN) = INT_MIN = -(INT_MAX+1)
template<> inline float Infinity() { return float(FLT_MAX); }
template<> inline float MinusInfinity() { return float(-FLT_MAX); }
template<> inline double Infinity() { return double(DBL_MAX); }
template<> inline double MinusInfinity() { return double(-DBL_MAX); }

//---- Rad/Deg Conversions
#define MAL_CONSTANT_PI (3.14159265359)
#define MAL_CONSTANT_RCP_SQRT_2 (1.0/1.41421356237)
#define MAL_CONSTANT_RCP_SQRT_3 (1.0/1.73205080757)
template <typename T> inline T Pi() { return T(MAL_CONSTANT_PI); }
template <typename T> inline T TwoPi() { return T(2.0*MAL_CONSTANT_PI); }
template <typename T> inline T HalfPi() { return T(0.5*MAL_CONSTANT_PI); }
template <typename T> inline T QuarterPi() { return T(0.25*MAL_CONSTANT_PI); }
template <typename T> inline T SixthPi() { return T(MAL_CONSTANT_PI/6.0); }
template <typename T> inline T RcpSqrtOfTwo() { return T(MAL_CONSTANT_RCP_SQRT_2); }
template <typename T> inline T RcpSqrtOfThree() { return T(MAL_CONSTANT_RCP_SQRT_3); }

template <typename T> inline T Deg2Rad( T value_deg ) { return value_deg * T(MAL_CONSTANT_PI) / T(180); }
template <typename T> inline T Rad2Deg( T value_rad ) { return value_rad * T(180) / T(MAL_CONSTANT_PI); }
template <typename T> inline T NormalizeRad( T value_rad )
{ if(value_rad<0) return value_rad + TwoPi<T>(); else if(value_rad >= TwoPi<T>()) return value_rad - TwoPi<T>(); else return value_rad; } //rad => [0,2Pi)

//-- Approx
template <typename T> inline T Epsilon() { return T(0.0001f); }
template<typename T >
inline bool ApproxEq( T value1, T value2, T epsilon = Epsilon<T>() ) { return Abs<T>( value1 - value2 ) <= epsilon ; }
template<typename T>
inline bool IsZero( T value, T epsilon = Epsilon<T>() ) { return Abs<T>( value ) <= epsilon; }

//--- NaN
template <typename T> inline bool IsNaN( T value ) { return value != value; } //!\todo Not sure if this includes #nan,#inf and #den
//\todo template <> inline bool IsNaN<float>( T value ) { return false; }
//\todo template <> inline bool IsNaN<double>( T value ) { return false; }

} //namespace mal

#endif //MAL_REAL_UTILS_H
