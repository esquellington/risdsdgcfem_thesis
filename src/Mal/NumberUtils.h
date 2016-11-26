#ifndef MAL_NUMBER_UTILS_H
#define MAL_NUMBER_UTILS_H

#include <Mal/Config.h>

//! Mal is MiniMal
namespace mal
{

//---- Boundaries
/*
template<typename T0, typename T1>
inline T0 Min( T0 value1, T1 value2 ) { return (value1<T0(value2)) ? value1 : T0(value2); }

template<typename T0, typename T1>
inline T0 Max( T0 value1, T1 value2 ) { return (value1>T0(value2)) ? value1 : T0(value2); }

template<typename T0, typename T1, typename T2>
inline T0 Clamp( T0 value, T1 min, T2 max ) { return Max( T0(min), Min( value, T0(max) )); } //!< Inside [min,max]
*/

/*---- Bounds
  IMPORTANT: We do NOT ACCEPT different operand types T0,T1,T2 to
  avoid Max(0,float) potential confusion, which truncates to int!
  \note std::max/min does NOT accept T0!=T1 either

  \todo Alternatively, use an "arithmetic promotion" as in
  krm::numeric so that <int,float> becomes float regardless of type
  order. \see http://blog.janmr.com/2010/08/cpp-templates-usual-arithmetic-conversions.html
*/
template<typename T>
inline T Min( T value1, T value2 ) { return (value1<value2) ? value1 : value2; }

template<typename T>
inline T Max( T value1, T value2 ) { return (value1>value2) ? value1 : value2; }

template<typename T>
inline T Clamp( T value, T min, T max ) { return Max( min, Min( value, max )); } //!< Inside [min,max]

template<typename T>
inline T Clamp01( T value ) { return Clamp( value, T(0), T(1) ); } //!< Inside [0,1]

/* \note This is SLOW, use standard abs/fabs in RealUtils.h
template<typename T>
inline T Abs( T value ) { return (value>T(0)) ? value : -value; }
*/

template<typename T>
inline T Sign( T value ) { return (value>=T(0)) ? T(1) : T(-1); }

//--- Ranges
template<typename T>
inline bool IsInRangeCO( T value, T min, T max ) { return (value >= min && value < max); } //!< Inside [min,max)

} //namespace mal

#endif //MAL_NUMBER_UTILS_H
