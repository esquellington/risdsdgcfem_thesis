#ifndef GEO_BV_TEST_OVERLAP_H
#define GEO_BV_TEST_OVERLAP_H

#include <Geo/bv/BoundingVolume.h>

namespace geo {
namespace bv {

bool TestOverlap( const BoundingVolume2* p_bv1, const BoundingVolume2* p_bv2 );
inline bool TestOverlap( const BoundingVolume2& bv1, const BoundingVolume2& bv2 ) { return TestOverlap(&bv1,&bv2); }

bool TestOverlap( const BoundingVolume3* p_bv1, const BoundingVolume3* p_bv2 );
inline bool TestOverlap( const BoundingVolume3& bv1, const BoundingVolume3& bv2 ) { return TestOverlap(&bv1,&bv2); }

//\todo These methods exist as legacy, because they are used in dimension-independent IBroadPhase. MAY disappear if IBroadPhase is splito into IBroadPhase2 and IBroadPhase3 someday
bool TestOverlap( const IBoundingVolume* p_bv1, const IBoundingVolume* p_bv2 );
inline bool TestOverlap( const IBoundingVolume& bv1, const IBoundingVolume& bv2 ) { return TestOverlap(&bv1,&bv2); }

}} //namespace geo::bv

#endif // GEO_BV_TEST_OVERLAP_H
