#ifndef GEO_BV_TEST_RAY_H
#define GEO_BV_TEST_RAY_H

#include <Geo/bv/BoundingVolume.h>
#include <Geo/np/RayCast.h>

namespace geo {
namespace bv {

bool TestRay( const BoundingVolume2* p_bv, const np::Ray2& ray );
inline bool TestRay( const BoundingVolume2& bv, const np::Ray2& ray ) { return TestRay(&bv,ray); }

bool TestRay( const BoundingVolume3* p_bv, const np::Ray3& ray );
inline bool TestRay( const BoundingVolume3& bv, const np::Ray3& ray ) { return TestRay(&bv,ray); }

}} //namespace geo::bv

#endif // GEO_BV_TEST_RAY_H
