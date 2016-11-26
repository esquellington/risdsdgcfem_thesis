#ifndef GEO_NP_INTERSECTION_H
#define GEO_NP_INTERSECTION_H

#include <Geo/Config.h>
#include "Context.h"
#include "stats.h"

#ifdef __GEO_ENABLE_NP_SIMD
#  include "SIMD.h"
#endif

namespace geo {
namespace np {

/*\todo Intersection functionality should follow Overlap/Contact
  conventions and be split into an Intersection.h header and
  implementation-specific Intersection_X_Y_Method files

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

//! Segment2 vs Segment 2D => (p0,p1) vs (q0,q1)
inline bool Intersection_Segment2_Segment2( const Vec2& p0, const Vec2& p1,
                                            const Vec2& q0, const Vec2& q1,
                                            Real& lambda1, Real& lambda2,
                                            const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC(intersection.m_Segment2_Segment2);
    Vec2 d1( p1 - p0 );
    Vec2 d2( q1 - q0 );
    Real det( d1[1]*d2[0] - d1[0]*d2[1] );
    if( mal::Abs(det) > Real(1e-5) ) //\todo epsilon parameter, or use np::Context as in Overlap/Contact!!
    {
        Vec2 diff_ab( q0 - p0 );
        Real inv_det( mal::Rcp(det) );
        lambda1 = inv_det * (diff_ab[1]*d2[0] - diff_ab[0]*d2[1]);
        lambda2 = inv_det * (diff_ab[1]*d1[0] - diff_ab[0]*d1[1]);
        return lambda1 >= Real(0) && lambda1 <= Real(1)
            && lambda2 >= Real(0) && lambda2 <= Real(1);
    }
    else // \todo Parallel, may be partially coincident, we MUST handle it
    {
        /*\todo We'll need to discriminate between coincident cases
          CROSSING | INSIDE&  NONCROSSING | OUTSIDE&  NONCROSSING
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

/* Intersection T1 vs T2
   2 essentially different cases are possible V-F and E-E, where the V-F case is asymmetric V1-F2 != V2-F1
   - V-F: One vertex T1.V_i is behind the other triangle face T2.F, with 2 consecutive edges T1.E_i,T1.E_i+1 piercing T2.F
   - E-E: One eddge of each triangle T1.E_i,T2.E_j pierces the other triangle face T2.F,T1.F
   Thus, we always have 2 edges piercing 1 or 2 faces, the common data required to describe the intersection is:
   - V1-F2 | V2-F1 | E1-E2 case: 3 values == 2 bits
   - EdgeID1, EdgeID2: IT has value 0..2 ==> 2 bits
   ==> Both case and can be simplified using an EdgeID that contains the triangle 0/1 bit, thus if both IT2T_EdgeID have the same IT2T_TriangleID, we have a V-F case, otherwise an E-E case with E1-E2 implicit order
   - LambdaE1, LambdaE2
   - PointE1, PointE2: Points

   IT2T_EdgeID => 3b : 1b T1/T2 {0,1} | 2b EdgeInTriangle {0,1,2}

   V-F: The V is inbetween E0 and E1, and its VertexInTriangle index 0..2 corresponds to m_EdgeInT1 0..2
*/
struct intersection_triangle3_triangle3_result_type
{
    inline bool IsVF() const { return m_T0 == m_T1 && m_T0 == 0; }
    inline bool IsFV() const { return m_T0 == m_T1 && m_T0 == 1; }
    inline bool IsEE() const { return m_T0 != m_T1; } //implicitly T0=0 && T1=1
    inline uint32 GetT0() const { return m_T0; }
    inline uint32 GetT1() const { return m_T1; }
    inline uint32 GetEIT0() const { return m_EdgeInT0; }
    inline uint32 GetEIT1() const { return m_EdgeInT1; }
    inline void SetE( uint32 index, uint32 t, uint32 eit, Real lambda, const Vec3& point )
        {
            if( index == 0 ) { m_T0 = t; m_EdgeInT0 = eit; m_Lambda0 = lambda; m_Point0 = point; }
            else { m_T1 = t; m_EdgeInT1 = eit; m_Lambda1 = lambda; m_Point1 = point; }
        }
    /*
    inline void Flip()
        {
            std::swap(m_Point0,m_Point1);
            std::swap(m_Lambda0,m_Lambda1);
            uint8 tmp( m_T0 ); m_T0 = m_T1; m_T1 = tmp;
            tmp = m_EdgeInT0; m_EdgeInT0 = m_EdgeInT1; m_EdgeInT1 = tmp;
        }
    */

    Vec3 m_Point0, m_Point1;
    Real m_Lambda0, m_Lambda1;
    uint8 m_T0       :1;
    uint8 m_T1       :1;
    uint8 m_EdgeInT0 :2;
    uint8 m_EdgeInT1 :2;
};

bool Intersection_Triangle3_Triangle3( const Vec3& p0, const Vec3& p1, const Vec3& p2,
                                       const Vec3& q0, const Vec3& q1, const Vec3& q2,
                                       intersection_triangle3_triangle3_result_type& ittr,
                                       const Context* p_context = g_pDefaultContext );

#ifdef __GEO_ENABLE_NP_SIMD
bool Intersection_Triangle3_Triangle3_SIMD( const simd::V3f& p0, const simd::V3f& p1, const simd::V3f& p2,
                                            const simd::V3f& q0, const simd::V3f& q1, const simd::V3f& q2,
                                            intersection_triangle3_triangle3_result_type& ittr,
                                            const Context* p_context = g_pDefaultContext );
#endif

}} //namespace geo::np

#endif // GEO_NP_INTERSECTION_H
