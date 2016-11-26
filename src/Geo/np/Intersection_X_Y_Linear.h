#ifndef GEO_NP_INTERSECTION_X_Y_LINEAR_H
#define GEO_NP_INTERSECTION_X_Y_LINEAR_H

#include <Geo/Config.h>
#include "stats.h"

namespace geo {
namespace np {

/*\todo
  //2D
  Intersection_Line2_Line2
  Intersection_Line2_Segment2
  Intersection_Line2_Triangle2
  + Intersection_Segment2_Segment2
  Intersection_Segment2_Triangle2
  Intersection_Triangle2_Triangle2

  //3D
  Intersection_Line3_Line3
  Intersection_Line3_Segment3
  Intersection_Line3_Triangle3
  Intersection_Line3_Tetrahedron3
  Intersection_Segment3_Segment3
  Intersection_Segment3_Triangle3
  Intersection_Segment3_Tetrahedron3
  Intersection_Triangle3_Triangle3
  Intersection_Triangle3_Tetrahedron3
  Intersection_Tetrahedron3_Tetrahedron3
*/

//! Segment2 vs Segment 2D => (a1,b1) vs (a2,b2)
inline bool Intersection_Segment2_Segment2( const Vec2 &a1, const Vec2 &b1,
                                            const Vec2 &a2, const Vec2 &b2,
                                            //Real epsilon_length,
                                            Real &lambda1, Real &lambda2 )
{
    GEO_NP_STAT_INC(intersection.m_Segment2_Segment2);
    Vec2 d1( b1 - a1 );
    Vec2 d2( b2 - a2 );
    Real det( d1[1]*d2[0] - d1[0]*d2[1] );
    if( mal::Abs(det) > Real(0.00001f) ) //\todo epsilon parameter, or pass np::Context as in Overlap/Contact!!
    {
        Vec2 diff_ab( a2 - a1 );
        Real inv_det( mal::Rcp(det) );
        lambda1 = inv_det * (diff_ab[1]*d2[0] - diff_ab[0]*d2[1]);
        lambda2 = inv_det * (diff_ab[1]*d1[0] - diff_ab[0]*d1[1]);
        return lambda1 >= Real(0) && lambda1 <= Real(1)
            && lambda2 >= Real(0) && lambda2 <= Real(1);
    }
    else // \todo Parallel, may be partially coincident, we MUST handle it
    {
        /*\todo We'll need to discriminate between coincident cases
          CROSSING | INSIDE & NONCROSSING | OUTSIDE & NONCROSSING
          within a certain epsilon_length, because we're interested in
          shallow penetrations, where segments will tend to ALIGN
          while inside or outside (depending on contact response
          method)... this behaviour should be CONSISTENT, because
          incremental IM computation requires finding an outwards
          crossing for each inwards one.
        */
        return false;
    }
}

}} //namespace geo::np

#endif // GEO_NP_INTERSECTION_X_Y_LINEAR_H
