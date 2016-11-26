#ifndef MAL_G_RANDOM_H
#define MAL_G_RANDOM_H

#include <pla_types.h>
#include <Mal/Config.h>
#include <Mal/GVec.h>
#include <stdlib.h>

namespace mal
{

inline void SetRandomSeed( int seed ) { srand(seed); }

//---- Uniform random distributions
template <typename T>
T Random01()
{
    return T( float(rand()) / float(RAND_MAX) );
}

template <typename T, typename T1, typename T2>
T RandomF( T1 min, T2 max )
{
    //\todo static_assert<IsFloat<T>>
    return T( T(min) + T(max-min)*Random01<T>() );
}

template <typename T, typename T1, typename T2>
T RandomI( T1 min, T2 max )
{
    //\todo static_assert<IsInteger<T>>
    return T( T(min) + T( rand() % (max-min+1) ) );
}

/*
template<> inline int8 RandomI( int8 min, int8 max ) { return min + int8( rand() % (max-min+1) ); }
template<> inline uint8 RandomI( uint8 min, uint8 max ) { return min + uint8( rand() % (max-min+1) ); }
template<> inline int16 RandomI( int16 min, int16 max ) { return min + int16( rand() % (max-min+1) ); }
template<> inline uint16 RandomI( uint16 min, uint16 max ) { return min + uint16( rand() % (max-min+1) ); }
template<> inline int32 RandomI( int32 min, int32 max ) { return min + int32( rand() % (max-min+1) ); }
template<> inline uint32 RandomI( uint32 min, uint32 max ) { return min + uint32( rand() % (max-min+1) ); }
//\todo int64 WILL NOT GEN 64b values because rand() is 32b... see if there's any rand64...
template<> inline int64 RandomI( int64 min, int64 max ) { MAL_ASSERT(false); return min + int64( rand() % (max-min+1) ); }
template<> inline uint64 RandomI( uint64 min, uint64 max ) { MAL_ASSERT(false); return min + uint64( rand() % (max-min+1) ); }
*/

template <typename T, int N>
GVec<T,N> RandomV( const GVec<T,N> &vec_min, const GVec<T,N> &vec_max )
{
    GVec<T,N> res;
    for( int i=0; i<N; i++ ) res[i] = RandomF<T>(vec_min[i],vec_max[i]);
    return res;
}


//---- Radial random distributions

//! Call with RandomUnitVec<T,N>(min,max) for dimension N unitary vector
template <typename T, int N>
GVec<T,N> RandomUnitVec()
{
    GVec<T,N> dir;
    for( int i=0; i<N; i++ ) dir[i] = RandomF<T>(-1,1);
    if( dir.NormSq() == T(0) ) dir[0] = T(1);
    return dir.Normalized();
}

//! Call with RandomRadialVec<N>(min,max) for dimension N vector
template <int N, typename T>
GVec<T,N> RandomRadialVec( T radius_min, T radius_max )
{
    // Random radius
    return RandomF<T>(radius_min,radius_max) * RandomUnitVec<T,N>();
}


} // namespace mal

#endif //MAL_G_RANDOM_H
