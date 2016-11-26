#ifndef GEO_MP_TEST_RAY_H
#define GEO_MP_TEST_RAY_H

#include <Geo/IObject.h>
#include <Geo/np/RayCast.h>

namespace geo {
namespace mp {

bool TestRay( const IObject2* p_obj, const np::Ray2& ray, np::RayHit2& rh, np::RayCache2* p_rc );
bool TestRay( const IObject3* p_obj, const np::Ray3& ray, np::RayHit3& rh, np::RayCache3* p_rc );

}} //namespace geo::mp

#endif // GEO_MP_TEST_RAY_H
