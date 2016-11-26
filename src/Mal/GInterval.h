#ifndef MAL_GINTERVAL_H
#define MAL_GINTERVAL_H

#include <Mal/Config.h>
#include <Mal/RealUtils.h>

/* Simple [min,max] interval representation
   - empty/infinite intervals are represented as [inf,-inf], [-inf,inf]
   - no need to check for empty/infinity in Merge,Intersect or TestOverlap
   - IMPORTANT T can be integral (int) but NOT unsigned (uint) due to empty/infinite interval conventions
*/
namespace mal
{

template<typename T>
class GInterval
{
public:
    typedef T real_type;
    finline static GInterval Empty() { return GInterval(mal::Infinity<T>(),-mal::Infinity<T>()); }
    finline static GInterval Infinity() { return GInterval(-mal::Infinity<T>(),mal::Infinity<T>()); }

public:
    //!\name Construction
    //@{
    finline GInterval() : m_Min(mal::Infinity<T>()), m_Max(-mal::Infinity<T>()) {}
    template <typename S>
    finline explicit GInterval( S val ) : m_Min(T(val)), m_Max(T(val)) {}
    template <typename S0, typename S1>
    finline explicit GInterval( S0 min, S1 max ) : m_Min(T(min)), m_Max(T(max)) {}
    template <typename S>
    finline GInterval( const GInterval<S> &other ) : m_Min(T(other.m_Min)), m_Max(T(other.m_Max)) {}

    template <typename S0>
    finline GInterval& Set( S0 val ) { m_Min = T(val); m_Max = T(val); return *this; }
    template <typename S0, typename S1>
    finline GInterval& SetMinMax( S0 min, S1 max ) { m_Min = T(min); m_Max = T(max); return *this; }
    template <typename S0, typename S1>
    finline GInterval& SetCenterHalfSizes( S0 center, S0 hs ) { MAL_ASSERT(T(hs)>=0); m_Min = T(center)-T(hs); m_Max = T(center)+T(hs); return *this; }
    //@}

    //!\name Access
    //@{
    finline bool IsEmpty() const { return m_Max < m_Min; }
    finline bool IsInfinity() const { return m_Min == -mal::Infinity<T>() && m_Max == mal::Infinity<T>(); }
    finline T &Min() { return m_Min; }
    finline T &Max() { return m_Max; }
    finline const T &Min() const { return m_Min; }
    finline const T &Max() const { return m_Max; }
    finline T Mid() const { return T(0.5) * (m_Max+m_Min); }
    //@}

    //!\name Modifiers \note They do not consider empty/infinite explicitly but WORK if represented as [inf,-inf], [-inf,inf]
    //\todo Extend() BREAKS infinity and may break empty (if val is LARGE)
    finline GInterval& Extend( T val ) { m_Min -= val; m_Max += val; return *this; }
    //\note [a,b] + [c,d] = [min(a,b),max(b,c)] => works for empty [inf,-inf] and infinity [-inf,inf]
    finline GInterval& Merge( T val ) { if(val<m_Min) m_Min = val; if(val>m_Max) m_Max = val; return *this; }
    finline GInterval& Merge( const GInterval &other ) { if(other.m_Min < m_Min) m_Min = other.m_Min; if(other.m_Max > m_Max) m_Max = other.m_Max; return *this; }
    //\note [a,b] & [c,d] = [max(a,c),min(b,d)] => works for empty [inf,-inf] and infinity [-inf,inf]
    finline GInterval& Intersect( const GInterval &other ) { if( m_Min < other.m_Min ) m_Min = other.m_Min; if( m_Max > other.m_Max ) m_Max = other.m_Max; return *this; }
    //@}

    //!\name Norms
    //@{
    finline T Length() const { if(!IsEmpty() && !IsInfinity()) return m_Max-m_Min; else return T(0); }
    //@}

    //!\name Utility
    //@{
    finline T Lerp( T lambda ) const { return lambda*m_Max + (T(1)-lambda)*m_Min; } //\todo WILL FAIL iff IsIntegral<T>
    finline bool TestOverlap( T value ) const { return (value >= m_Min) && (value <= m_Max); }
    finline bool TestOverlap( const GInterval &other ) const { return (m_Min <= other.m_Max) && (other.m_Min <= m_Max); }
    //@}

private:
    T m_Min, m_Max;
};

} //namespace mal

#ifdef __DISABLED_OLD_CODE
/*\todo This is DISABLED because, even if faster for non-empty
  intervals, there are SERIOUS problems and ambiguities between EMPTY
  and 0-LENGTH intervals that caused ugly bugs at DLE. Artificially
  enlarging 0-LENGTH intervals to avoid mistaking them for empty
  min==max ones is A BAD IDEA (caused raycast bugs in
  broadphase). TESTOVERLAP_MUL seems a good idea, but requires
  0-LENGTH intervals to be considered empty. THEREFORE, by now I
  disable these "optimizations" and leave plain-old max<min empty
  interval representation, as well as comparison-based
  TestOverlap. Slower but robust.

  A possible way to improve TestOverlap and Empty interval operation
  would be having a specific GInterval_NonEmpty class, that would have
  optimized TestOverlap (50% faster) and could be promoted to the
  generic but slower GInterval version.

#define __MAL_INTERVAL_TESTOVERLAP_MUL //\todo Avoids comparisons, uses mul to test overlap, requires __MAL_INTERVAL_EMPTY_IMPLIES_MIN_EQUAL_MAX to support empty intervals properly
#define __MAL_INTERVAL_EMPTY_IMPLIES_MIN_EQUAL_MAX //\todo This is required for BREAKS GRayCast_Capsule...

\todo THIS WAS ALL BULLSHIT, everything becomes SO MUCH SIMPLER if
      empty/infinite intervals are represented as [inf,-inf], [-inf,inf], no
      need to check for empty/infinity in Merge,Intersect or TestOverlap AT
      ALL
*/

namespace mal
{

template<typename T>
class GInterval
{
public:
    typedef T real_type;
    finline static GInterval Empty() { return GInterval(); }
    finline static GInterval Infinity() { return GInterval(-mal::Infinity<T>(),mal::Infinity<T>()); }

public:
    //!\name Construction
    //@{
#ifdef __MAL_INTERVAL_EMPTY_IMPLIES_MIN_EQUAL_MAX
    finline GInterval() : m_Min( Infinity<T>() ), m_Max( Infinity<T>() ) {}
#else
    finline GInterval() : m_Min(T(0)), m_Max(T(-1)) {}
#endif

#ifdef __MAL_INTERVAL_EMPTY_IMPLIES_MIN_EQUAL_MAX
    template <typename S>
    finline explicit GInterval( S val ) : m_Min(T(val)-Epsilon<T>()), m_Max(T(val)+Epsilon<T>()) {} //\note single value +- epsilon because (value,value) would be empty \todo THIS IS A BAAAAAD IDEA!
#else
    template <typename S>
    finline explicit GInterval( S val ) : m_Min(T(val)), m_Max(T(val)) {}
#endif
    template <typename S0, typename S1>
    finline explicit GInterval( S0 min, S1 max ) : m_Min(T(min)), m_Max(T(max)) {}
    template <typename S>
    finline GInterval( const GInterval<S> &other ) : m_Min(T(other.m_Min)), m_Max(T(other.m_Max)) {}

#ifdef __MAL_INTERVAL_EMPTY_IMPLIES_MIN_EQUAL_MAX
    template <typename S0>
    finline GInterval& Set( S0 val ) { m_Min = T(val)-Epsilon<T>(); m_Max = T(val)+Epsilon<T>();  return *this; } //\note single value +- epsilon because (value,value) would be empty
#else
    template <typename S0>
    finline GInterval& Set( S0 val ) { m_Min = T(val); m_Max = T(val);  return *this; }
#endif
    template <typename S0, typename S1>
    finline GInterval& SetMinMax( S0 min, S1 max ) { m_Min = T(min); m_Max = T(max);  return *this; }
    template <typename S0, typename S1>
    finline GInterval& SetCenterHalfSizes( S0 center, S0 hs ) { MAL_ASSERT(T(hs)>=0); m_Min = T(center)-T(hs); m_Max = T(center)+T(hs); return *this; }
    //@}

    //!\name Access
    //@{
#ifdef __MAL_INTERVAL_EMPTY_IMPLIES_MIN_EQUAL_MAX
    finline bool IsEmpty() const { return m_Min == Infinity<T>(); } //\note also, m_Min == m_Max, necessary but NOT sufficient
#else
    finline bool IsEmpty() const { return m_Max < m_Min; }
#endif
    finline T &Min() { return m_Min; }
    finline T &Max() { return m_Max; }
    finline const T &Min() const { return m_Min; }
    finline const T &Max() const { return m_Max; }
    finline T Mid() const { return T(0.5) * (m_Max+m_Min); }
    //@}

    //!\name Add/Sub
    //@{
    /*
    finline GInterval operator+( const GInterval &v ) const { GInterval sum; for(int i=0; i<N; i++) sum[i] = data[i] + v[i]; return sum; }
    finline GInterval operator-( const GInterval &v ) const { GInterval sub; for(int i=0; i<N; i++) sub[i] = data[i] - v[i]; return sub; }
    finline GInterval operator+=( const GInterval &v ) { for(int i=0; i<N; i++) data[i] += v[i]; return *this; }
    finline GInterval operator-=( const GInterval &v ) { for(int i=0; i<N; i++) data[i] -= v[i]; return *this; }
    */

    //\todo Extend, Merge, Intersect do NOT CONSIDER empty/infinite intervals explicitly but WORK if empty/infinite interval is represented as [inf,-inf], [-inf,inf]
    // CAN FAIL for both empty [inf,-inf] and infinite [-inf,inf] due to numeric infinity handling... consider using NaN inf?!?!
    finline GInterval& Extend( T val ) { m_Min -= val; m_Max += val; return *this; }
    //\note [a,b] + [c,d] = [min(a,b),max(b,c)] => works for empty [inf,-inf] and infinity [-inf,inf]
    finline GInterval& Merge( T val ) { if(val<m_Min) m_Min = val; else if(val>m_Max) m_Max = val; return *this; }
    finline GInterval& Merge( const GInterval &other ) { if(other.m_Min < m_Min) m_Min = other.m_Min; if(other.m_Max > m_Max) m_Max = other.m_Max; return *this; }
    //\note [a,b] & [c,d] = [max(a,c),min(b,d)] => works for empty [inf,-inf] and infinity [-inf,inf]
    finline GInterval& Intersect( const GInterval &other ) { if( m_Min < other.m_Min ) m_Min = other.m_Min; if( m_Max > other.m_Max ) m_Max = other.m_Max; return *this; }
    //@}

    //!\name Products
    //@{
    /*
    finline GInterval operator*( T val ) const { GInterval prod; for(int i=0; i<N; i++) prod[i] = val*data[i]; return prod; }
    finline GInterval operator*=( T val ) { for(int i=0; i<N; i++) data[i] *= val; return *this; }
    friend finline GInterval operator*( T val, const GInterval &v ) { return v*val; }

    finline GInterval operator/( T val ) const { T inv_val = T(1)/val; return (*this) * inv_val; }
    finline GInterval operator/=( T val ) { T inv_val = T(1)/val; return (*this)*= inv_val; }

    finline T operator*( const GInterval &v ) const { T acc(0); for(int i=0; i<N; i++) acc += data[i] * v[i]; return acc; }
    finline GInterval operator%( const GInterval &v ) const { MAL_STATIC_ASSERT(N==3); return GInterval( data[1]*v[2] - data[2]*v[1],
                                                                                          data[2]*v[0] - data[0]*v[2],
                                                                                          data[0]*v[1] - data[1]*v[0] ); }
    */
    //@}

    //!\name Norms
    //@{
    finline T Length() const { return m_Max-m_Min; }
    //@}

    //!\name Utility
    //@{
    finline GInterval Lerp( T lambda ) const { return lambda*m_Max + (T(1)-lambda)*m_Min; }

#ifdef __MAL_INTERVAL_TESTOVERLAP_MUL
    finline bool TestOverlap( T value ) const { return !IsEmpty() && (m_Min - value) * (value - m_Max) > T(0); } //optimized as OverlapsMul
    finline bool TestOverlap( const GInterval &other ) const { return !IsEmpty() && !other.IsEmpty() && TestOverlapMul_NoEmpty(other); } //_NoEmpty is 50% faster...
#else
    finline bool TestOverlap( T value ) const { return !IsEmpty() && (value >= m_Min) && (value <= m_Max); }
    finline bool TestOverlap( const GInterval &other ) const { return TestOverlapCmp(other); }
#endif
    //@}

public:
    //\name TEMPORAL: Overlaps() efficiency tests, only fastest version will survive
    //@{
    finline bool TestOverlapMul_NoEmpty( const GInterval &other ) const { return (m_Min - other.m_Max) * (other.m_Min - m_Max) >= T(0); } //Fastestest NoEmpty version, 50% faster than empty-aware version

    finline bool TestOverlapMul( const GInterval &other ) const { return (m_Max - m_Min) //=> !IsEmpty()
                                                                         * (other.m_Max - other.m_Min) //=> !other.IsEmpty()
                                                                         * (m_Min - other.m_Max)
                                                                         * (other.m_Min - m_Max) > T(0); } //Fastestest among empty-aware methods, requires IsEmpty() => min == max => max-min == 0, regardless of actual min,max values
    finline bool TestOverlapCmp_NoEmpty( const GInterval &other ) const { return (m_Min <= other.m_Max) && (other.m_Min <= m_Max); } //almost 2x slower than Mul
    finline bool TestOverlapCmp( const GInterval &other ) const { return !IsEmpty() && !other.IsEmpty() && (m_Min <= other.m_Max) && (other.m_Min <= m_Max); }
    finline bool TestOverlapAbs_NoEmpty( const GInterval &other ) const { // 0.5 factor on both sides simplified, approx as fast as Mul (only if mal::Abs uses fabs)
        return mal::Abs( m_Max + m_Min - other.m_Max - other.m_Min )
            <= ( m_Max - m_Min + other.m_Max - other.m_Min ); }
    finline bool TestOverlapAbs( const GInterval &other ) const { // 0.5 factor on both sides simplified,  approx as fast as Mul (only if mal::Abs uses fabs)
        return !IsEmpty() && !other.IsEmpty() && mal::Abs( m_Max + m_Min - other.m_Max - other.m_Min )
            <= ( m_Max - m_Min + other.m_Max - other.m_Min ); }
    //@}

private:
    T m_Min, m_Max;
};

} //namespace mal

#endif //__DISABLED_OLD_CODE

#endif //MAL_GINTERVAL_H
