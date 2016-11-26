#ifndef GEO_NP_CLIP_H
#define GEO_NP_CLIP_H

#include <Geo/Config.h>
#include <Geo/np/Context.h>

/*! This file contains the prototypes for all default Clip_X_Y functionality.
*/

namespace geo {
namespace np {

/* Clip a Triangle3 against a Plane3
   - Computes a convex polygon resulting from the clipping, with the same orientation as the input triangle.
   - Returns the clipped polygon vertex count [3..4], or 0 if no intersection
   \pre vec_clipped_points MUST have, at least, 4 valid slots
*/
unsigned int Clip_Triangle3_Plane3( const Vec3& tri0, const Vec3& tri1, const Vec3& tri2,
                                    const Vec3& plane_normal, Real plane_coeff_d,
                                    Vec3* vec_polygon_clipped_pos,
                                    const Context* p_context = g_pDefaultContext );

/* Clip a Triangle3 against a Tetrahedron3:
   - Computes a convex polygon resulting from the clipping, with the same orientation as the input triangle.
   - Returns the clipped polygon vertex count [3..6], or 0 if no intersection
   \pre vec_clipped_points MUST have, at least, 6 valid slots
*/
unsigned int Clip_Triangle3_Tetrahedron3( const Vec3& tri0, const Vec3& tri1, const Vec3& tri2,
                                          const Vec3& tet0, const Vec3& tet1, const Vec3& tet2, const Vec3& tet3,
                                          Vec3* vec_polygon_clipped_pos,
                                          const Context* p_context = g_pDefaultContext );

}} //namespace geo::np

#endif // GEO_NP_CLIP_H
