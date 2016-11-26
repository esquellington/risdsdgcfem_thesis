#include "Intersection.h"
#include "RayCast.h"

#ifdef __GEO_ENABLE_NP_SIMD
#  include "SIMD.h"
#  define __ENABLE_SIMD_INTERSECTION_TRI_TRI
#endif

namespace geo {
namespace np {

/* Specialization of RayCast_Triangle3_DoubleSided for
   - Ray == Segment
   - Triangle normal passed explicitly
   - Only rh lambda and point are computed
*/
inline bool RayCast_Triangle3_DoubleSided_FAST_SEGMENT_EXPLICIT_EDGES( const Vec3& segment_p, const Vec3& segment_v, //segment_v == pq == -qp
                                                                       const Vec3& triangle_p0, const Vec3& triangle_v01, const Vec3& triangle_v02, const Vec3& triangle_normal_non_unitary,
                                                                       Vec3& point, Real& lambda,
                                                                       const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( raycast.m_Triangle3_DS );
    Vec3 ab( triangle_v01 );
    Vec3 ac( triangle_v02 );
    //segment_v == Vec3 pq( segment_q - segment_p ); //Vec3 qp( segment_p - segment_q );
    Vec3 n( triangle_normal_non_unitary );
    Real d( -mal::Dot( segment_v, n ) ); //== mal::Dot( qp, n )
    // Coplanar?
    //if( mal::Abs(d) < g_pDefaultContext->m_Epsilon_Dir ) return false;
    if( d <= Real(0) )
    {
        // Invert triangle orientation
        n = -n;
        ab = triangle_v02;
        ac = triangle_v01;
        d = -d;
    }
    Vec3 ap( segment_p - triangle_p0 );
    Real t( mal::Dot( ap, n ) );
    if( t < Real(0) ) return false;
    if( t > d ) return false;
    Vec3 e( mal::Cross( ap, segment_v ) ); //== mal::Cross( qp, ap ) );
    Real v( mal::Dot( ac, e ) );
    if( v < Real(0) || v > d ) return false;
    Real w( -mal::Dot( ab, e ) );
    if( w < Real(0) || v + w > d ) return false;
    // Hit found, fill data and return true
    lambda = t/d;
    point = segment_p + lambda * segment_v; //== p - lambda * qp
    return true;
}


#ifdef __ENABLE_SIMD_INTERSECTION_TRI_TRI
/* Specialization of RayCast_Triangle3_DoubleSided for
   - Ray == Segment
   - Triangle normal passed explicitly
   - Only rh lambda and point are computed
*/
inline bool RayCast_Triangle3_DoubleSided_FAST_SEGMENT_EXPLICIT_EDGES_SIMD( const simd::V3f& segment_p, const simd::V3f& segment_v, //segment_v == pq == -qp
                                                                            const simd::V3f& triangle_p0, const simd::V3f& triangle_v01, const simd::V3f& triangle_v02, const simd::V3f& triangle_normal_non_unitary,
                                                                            Vec3& point, Real& lambda,
                                                                            const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( raycast.m_Triangle3_DS );
    simd::V3f ab( triangle_v01 );
    simd::V3f ac( triangle_v02 );
    //segment_v == simd::V3f pq( segment_q - segment_p ); //simd::V3f qp( segment_p - segment_q );
    simd::V3f n( triangle_normal_non_unitary );
    Real d( -simd::Dot( segment_v, n ) ); //== simd::Dot( qp, n )
    // Coplanar?
    //if( simd::Abs(d) < g_pDefaultContext->m_Epsilon_Dir ) return false;
    if( d <= Real(0) )
    {
        // Invert triangle orientation
        n = -n;
        ab = triangle_v02;
        ac = triangle_v01;
        d = -d;
    }
    simd::V3f ap( segment_p - triangle_p0 );
    Real t( simd::Dot( ap, n ) );
    if( t < Real(0) ) return false;
    if( t > d ) return false;
    simd::V3f e( simd::Cross( ap, segment_v ) ); //== simd::Cross( qp, ap ) );
    Real v( simd::Dot( ac, e ) );
    if( v < Real(0) || v > d ) return false;
    Real w( -simd::Dot( ab, e ) );
    if( w < Real(0) || v + w > d ) return false;
    // Hit found, fill data and return true
    lambda = t/d;
    point = (segment_p + lambda * segment_v).AsVec3f(); //== p - lambda * qp
    return true;
}
#endif //__ENABLE_SIMD_INTERSECTION_TRI_TRI

/* Intersection T1 vs T2
   2 essentially different cases are possible V-F and E-E, where the V-F case is asymmetric V1-F2 != V2-F1
   - V-F: One vertex T1.V_i is behind the other triangle face T2.F, with 2 consecutive edges T1.E_i,T1.E_i+1 piercing T2.F
   - E-E: One eddge of each triangle T1.E_i,T2.E_j pierces the other triangle face T2.F,T1.F
   Thus, we always have 2 edges piercing 1 or 2 faces, the common data required to describe the intersection is:
   - V1-F2 | V2-F1 | E1-E2 case: 3 values == 2 bits
   - EdgeID1, EdgeID2: IT has value 0..2 ==> 2 bits
   ==> Both case and can be simplified using an EdgeID that contains
       the triangle 0/1 bit, thus if both IT2T_EdgeID have the same
       IT2T_TriangleID, we have a V-F case, otherwise an E-E case with
       E1-E2 implicit order
   - LambdaE1, LambdaE2
   - PointE1, PointE2: Points, can be derived from previous info, so
     it MAY NOT be part of
     intersection_triangle3_triangle3_result_type, but should be
     returned by Intersection_Triangle3_Triangle3() in any case.

   IT2T_EdgeID => 3b : 1b T1/T2 {0,1} | 2b EdgeInTriangle {0,1,2}

   V-F: The V is inbetween E0 and E1, and its VertexInTriangle index 0..2 corresponds to m_EdgeInT1 0..2

   \todo The cost of function is CRITICAL in
   Contact_DCRTS3_Plane3_BVH_TSS_Lazy_DCR(), and could be GREATLY
   REDUCED using early outs, and factoring/unrolling per-edge raycast
   computations, or using SIMD or other low-level optimizations.

   OPTIMIZATIONS:
   - DS double sided
   - C/I coherent/incoherent
   - FST: Fast segment raycast, implies DSI
   - EATP: Early out using triangle planes

   RESULTS for test_k3d.skr

scn_k3d_tetsolid( name="TetSolid.Box",
                  p=[0,2,0],
                  file = "s2s/sl/box3d.txt" );
scn_k3d_tetsolid( name = "TetSolid.Armadillo",
                  p  = [0.66,2,0],
                  file = "s2s/sl/armadillo_10K.txt" );
scn_k3d_tetsolid( name = "TetSolid.Dragon",
                  p  = [1.33,2,0],
                  file = "s2s/sl/dragon_10K.txt" );
scn_k3d_tetsolid( name = "TetSolid.Bunny",
                  p  = [2,2,0],
                  file = "s2s/sl/bunny3d.txt" );

   - Box-Armadillo (317 hit)
     - BF 270 DSC, 180 DSI, 120 FST, 27 EATP !!!!
     - BDT 170 DSC, 125 DSI, 90 FST, 35 EATP => 23 EATP if MaxLeafTests = 32
   - Armadillo-Dragon (370 hit)
     - BF 825 DSC, 570 DSI, 361 FST, 60 EATP
     - BDT 55 DSC, 40 DSI, 32 FST, 16 EATP ==> 12 if MaxLeafs = 8
   - Dragon-Bunny (603 hit)
     - BF 250, 170 DSC, 110 FST, 23 EATP
     - BDT 55 DSC, 40 DSI, 31 FST, 17 EATP
   DISCUSSION: (HD/MD/LD = high/mid/low def)
   - BDT slows down LD2HD case (Box is LARGE and overlaps a good part of the Armadillo)
     - The cost of BDOP unpacking, AABB(BDOP) conversion and AABB2AABB tests hits hard here
     => If we increase MaxLeafTests the cost DECREASES noticeably, below the BF cost
   - BDT improves dramatically the HD2HD case (Armadillo/Dragon)
     => Increasing MaxLeafTests does not alter performance noticeably
   - BDT improves HD2MD
     => Increasing MaxLeafTests there is further improvement

   \see RTCD 5.2.10
*/
bool Intersection_Triangle3_Triangle3( const Vec3 &p0, const Vec3 &p1, const Vec3 &p2,
                                       const Vec3 &q0, const Vec3 &q1, const Vec3 &q2,
                                       intersection_triangle3_triangle3_result_type& ittr,
                                       const Context* p_context )
{
#ifdef __ENABLE_SIMD_INTERSECTION_TRI_TRI
    return Intersection_Triangle3_Triangle3_SIMD( simd::V3f(p0), simd::V3f(p1), simd::V3f(p2),
                                                  simd::V3f(q0), simd::V3f(q1), simd::V3f(q2),
                                                  ittr, p_context );
#else
    GEO_NP_STAT_INC(intersection.m_Triangle3_Triangle3);

    // Compute triangle non-unit normals Early-out with triangle planes \todo MAYBE the 2xtriple-and can be unified/simplified somehow...
    const Vec3 q01(q1-q0);
    const Vec3 q02(q2-q0);
    const Vec3 nq( mal::Cross(q01,q02) );
    Real vec_dot_p_nq[3] = { mal::Dot(p0-q0,nq), mal::Dot(p1-q0,nq), mal::Dot(p2-q0,nq) };
    if( (vec_dot_p_nq[0] < 0 && vec_dot_p_nq[1] < 0 && vec_dot_p_nq[2] < 0)
        || (vec_dot_p_nq[0] > 0 && vec_dot_p_nq[1] > 0 && vec_dot_p_nq[2] > 0) ) return false;
    const Vec3 p01(p1-p0);
    const Vec3 p02(p2-p0);
    const Vec3 np( mal::Cross(p01,p02) );
    Real vec_dot_q_np[3] = { mal::Dot(q0-p0,np), mal::Dot(q1-p0,np), mal::Dot(q2-p0,np) };
    if( (vec_dot_q_np[0] < 0 && vec_dot_q_np[1] < 0 && vec_dot_q_np[2] < 0)
        || (vec_dot_q_np[0] > 0 && vec_dot_q_np[1] > 0 && vec_dot_q_np[2] > 0) ) return false;

    // Compute intersection dir d = np x nq and early-out projecting along it
    //=> This reduces RayCast calls to 1/2-1/3, but does NOT seem to speed up at all
    const Vec3 d( mal::Cross(np,nq) );
    Interval interval_p_d = Interval(mal::Dot(p0,d)).Merge(mal::Dot(p1,d)).Merge(mal::Dot(p2,d));
    Interval interval_q_d = Interval(mal::Dot(q0,d)).Merge(mal::Dot(q1,d)).Merge(mal::Dot(q2,d));
    if( !interval_p_d.TestOverlap(interval_q_d) ) return false;

    /* Most-parallel coord axis saves 6xDot => NOT faster than full-d
       version, but discards around 10% LESS and yields uglier code.
    int best_axis(0);
    Real abs_d_best(mal::Abs(d[0]));
    Real abs_d1(mal::Abs(d[1]));
    if( abs_d1 > abs_d_best ) { best_axis = 1; abs_d_best = abs_d1; }
    if( mal::Abs(d[2]) > abs_d_best ) best_axis = 2;
    Interval interval_p_d = Interval(p0[best_axis]).Merge(p1[best_axis]).Merge(p2[best_axis]);
    Interval interval_q_d = Interval(q0[best_axis]).Merge(q1[best_axis]).Merge(q2[best_axis]);
    if( !interval_p_d.TestOverlap(interval_q_d) ) return false;
    */

    // No early out left, perform 6 RC
    uint32 num_hits(0);
    Vec3 point;
    Real lambda;
    // Tri0 edges vs Tri1
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT_EXPLICIT_EDGES( p0, p01, q0, q01, q02, nq, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 0, 0, lambda, point ); }
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT_EXPLICIT_EDGES( p1, p2-p1, q0, q01, q02, nq, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 0, 1, lambda, point ); }
    if( num_hits == 2 ) return true;
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT_EXPLICIT_EDGES( p2, -p02, q0, q01, q02, nq, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 0, 2, lambda, point ); }
    if( num_hits == 2 ) return true;
    // Tri1 edges vs Tri0
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT_EXPLICIT_EDGES( q0, q01, p0, p01, p02, np, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 1, 0, lambda, point ); }
    if( num_hits == 2 ) return true;
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT_EXPLICIT_EDGES( q1, q2-q1, p0, p01, p02, np, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 1, 1, lambda, point ); }
    if( num_hits == 2 ) return true;
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT_EXPLICIT_EDGES( q2, -q02, p0, p01, p02, np, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 1, 2, lambda, point ); }
    if( num_hits == 2 ) return true;

    // GEO_LOG_ERROR_IF( num_hits != 2, "early outs should have discarded this intersection!" );

    GEO_LOG_ERROR_IF( num_hits != 0 && num_hits != 2, "geo::np::Intersection_Triangle3_Triangle3() num_hits = %u", num_hits );
    return num_hits > 0;
#endif //__ENABLE_SIMD_INTERSECTION_TRI_TRI
}

bool Intersection_Triangle3_Triangle3_OLD( const Vec3 &p0, const Vec3 &p1, const Vec3 &p2,
                                           const Vec3 &q0, const Vec3 &q1, const Vec3 &q2,
                                           intersection_triangle3_triangle3_result_type& ittr,
                                           const Context* p_context )
{
    GEO_NP_STAT_INC(intersection.m_Triangle3_Triangle3);

    // Compute triangle non-unit normals Early-out with triangle planes \todo MAYBE the 2xtriple-and can be unified/simplified somehow...
    const Vec3 nq( mal::Cross(q1-q0,q2-q0) );
    Real vec_dot_p_nq[3] = { mal::Dot(p0-q0,nq), mal::Dot(p1-q0,nq), mal::Dot(p2-q0,nq) };
    if( (vec_dot_p_nq[0] < 0 && vec_dot_p_nq[1] < 0 && vec_dot_p_nq[2] < 0)
        || (vec_dot_p_nq[0] > 0 && vec_dot_p_nq[1] > 0 && vec_dot_p_nq[2] > 0) ) return false;
    /*
    if( mal::Min( mal::Min(vec_dot_p_nq[0], vec_dot_p_nq[1]), vec_dot_p_nq[2])
        * mal::Max( mal::Max(vec_dot_p_nq[0], vec_dot_p_nq[1]), vec_dot_p_nq[2]) > 0 ) return false;
        */

    const Vec3 np( mal::Cross(p1-p0,p2-p0) );
    Real vec_dot_q_np[3] = { mal::Dot(q0-p0,np), mal::Dot(q1-p0,np), mal::Dot(q2-p0,np) };
    if( (vec_dot_q_np[0] < 0 && vec_dot_q_np[1] < 0 && vec_dot_q_np[2] < 0)
        || (vec_dot_q_np[0] > 0 && vec_dot_q_np[1] > 0 && vec_dot_q_np[2] > 0) ) return false;
    /*
    if( mal::Min( mal::Min(vec_dot_q_np[0], vec_dot_q_np[1]), vec_dot_q_np[2])
        * mal::Max( mal::Max(vec_dot_q_np[0], vec_dot_q_np[1]), vec_dot_q_np[2]) > 0 ) return false;
        */

    // No early out, perform 6 RC
    uint32 num_hits(0);
    Vec3 point;
    Real lambda;
    // Tri0 edges vs Tri1
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT( p0, p1, q0, q1, q2, nq, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 0, 0, lambda, point ); }
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT( p1, p2, q0, q1, q2, nq, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 0, 1, lambda, point ); }
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT( p2, p0, q0, q1, q2, nq, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 0, 2, lambda, point ); }
    // Tri1 edges vs Tri0
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT( q0, q1, p0, p1, p2, np, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 1, 0, lambda, point ); }
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT( q1, q2, p0, p1, p2, np, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 1, 1, lambda, point ); }
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT( q2, q0, p0, p1, p2, np, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 1, 2, lambda, point ); }

    /* Generic RayCast: 25-50% slower than FAST_SEGMENT
    RayHit3 rh;
    const Interval interval01(0,1);
    // Tri0 edges vs Tri1
    if( RayCast_Triangle3_DoubleSided( p0, p1-p0, interval01, q0, q1, q2, rh, p_context ) ) { ittr.SetE( num_hits++, 0, 0, rh.m_Interval.Min(), rh.m_Point ); }
    if( RayCast_Triangle3_DoubleSided( p1, p2-p1, interval01, q0, q1, q2, rh, p_context ) ) { ittr.SetE( num_hits++, 0, 1, rh.m_Interval.Min(), rh.m_Point ); }
    if( RayCast_Triangle3_DoubleSided( p2, p0-p2, interval01, q0, q1, q2, rh, p_context ) ) { ittr.SetE( num_hits++, 0, 2, rh.m_Interval.Min(), rh.m_Point ); }
    // Tri1 edges vs Tri0
    if( RayCast_Triangle3_DoubleSided( q0, q1-q0, interval01, p0, p1, p2, rh, p_context ) ) { ittr.SetE( num_hits++, 1, 0, rh.m_Interval.Min(), rh.m_Point ); }
    if( RayCast_Triangle3_DoubleSided( q1, q2-q1, interval01, p0, p1, p2, rh, p_context ) ) { ittr.SetE( num_hits++, 1, 1, rh.m_Interval.Min(), rh.m_Point ); }
    if( RayCast_Triangle3_DoubleSided( q2, q0-q2, interval01, p0, p1, p2, rh, p_context ) ) { ittr.SetE( num_hits++, 1, 2, rh.m_Interval.Min(), rh.m_Point ); }
    */

    /* COHERENT: 25-50% slower than non-coherent
    // Tri0 edges vs Tri1
    if( RayCast_Triangle3_DoubleSided_COHERENT( p0, p1-p0, interval01, q0, q1, q2, rh, p_context ) ) { ittr.SetE( num_hits++, 0, 0, rh.m_Interval.Min(), rh.m_Point ); }
    if( RayCast_Triangle3_DoubleSided_COHERENT( p1, p2-p1, interval01, q0, q1, q2, rh, p_context ) ) { ittr.SetE( num_hits++, 0, 1, rh.m_Interval.Min(), rh.m_Point ); }
    if( RayCast_Triangle3_DoubleSided_COHERENT( p2, p0-p2, interval01, q0, q1, q2, rh, p_context ) ) { ittr.SetE( num_hits++, 0, 2, rh.m_Interval.Min(), rh.m_Point ); }
    // Tri1 edges vs Tri0
    if( RayCast_Triangle3_DoubleSided_COHERENT( q0, q1-q0, interval01, p0, p1, p2, rh, p_context ) ) { ittr.SetE( num_hits++, 1, 0, rh.m_Interval.Min(), rh.m_Point ); }
    if( RayCast_Triangle3_DoubleSided_COHERENT( q1, q2-q1, interval01, p0, p1, p2, rh, p_context ) ) { ittr.SetE( num_hits++, 1, 1, rh.m_Interval.Min(), rh.m_Point ); }
    if( RayCast_Triangle3_DoubleSided_COHERENT( q2, q0-q2, interval01, p0, p1, p2, rh, p_context ) ) { ittr.SetE( num_hits++, 1, 2, rh.m_Interval.Min(), rh.m_Point ); }
    */

    GEO_LOG_ERROR_IF( num_hits != 0 && num_hits != 2, "geo::np::Intersection_Triangle3_Triangle3() num_hits = %u", num_hits );
    return num_hits > 0;
}

#ifdef __ENABLE_SIMD_INTERSECTION_TRI_TRI
bool Intersection_Triangle3_Triangle3_SIMD( const simd::V3f& p0, const simd::V3f& p1, const simd::V3f& p2,
                                            const simd::V3f& q0, const simd::V3f& q1, const simd::V3f& q2,
                                            intersection_triangle3_triangle3_result_type& ittr,
                                            const Context* p_context )
{
    GEO_NP_STAT_INC(intersection.m_Triangle3_Triangle3);
    // Compute triangle non-unit normals Early-out with triangle planes \todo MAYBE the 2xtriple-and can be unified/simplified somehow...
    const simd::V3f q01(q1 - q0);
    const simd::V3f q02(q2 - q0);
    const simd::V3f nq( simd::Cross(q01,q02) );
    Real vec_dot_p_nq[3] = { simd::Dot(p0-q0,nq),
                             simd::Dot(p1-q0,nq),
                             simd::Dot(p2-q0,nq) };
    if( (vec_dot_p_nq[0] < 0 && vec_dot_p_nq[1] < 0 && vec_dot_p_nq[2] < 0)
        || (vec_dot_p_nq[0] > 0 && vec_dot_p_nq[1] > 0 && vec_dot_p_nq[2] > 0) ) return false;
    const simd::V3f p01(p1-p0);
    const simd::V3f p02(p2-p0);
    const simd::V3f np( simd::Cross(p01,p02) );
    Real vec_dot_q_np[3] = { simd::Dot(q0-p0,np),
                             simd::Dot(q1-p0,np),
                             simd::Dot(q2-p0,np) };
    if( (vec_dot_q_np[0] < 0 && vec_dot_q_np[1] < 0 && vec_dot_q_np[2] < 0)
        || (vec_dot_q_np[0] > 0 && vec_dot_q_np[1] > 0 && vec_dot_q_np[2] > 0) ) return false;

    // Compute intersection dir d = np x nq and early-out projecting along it
    //=> This reduces RayCast calls to 1/2-1/3, but does NOT seem to speed up at all
    const simd::V3f d( simd::Cross(np,nq) );
    Interval interval_p_d = Interval(simd::Dot(p0,d)).Merge(simd::Dot(p1,d)).Merge(simd::Dot(p2,d));
    Interval interval_q_d = Interval(simd::Dot(q0,d)).Merge(simd::Dot(q1,d)).Merge(simd::Dot(q2,d));
    if( !interval_p_d.TestOverlap(interval_q_d) ) return false;

    /* Most-parallel coord axis saves 6xDot => NOT faster than full-d
       version, but discards around 10% LESS and yields uglier code.
    int best_axis(0);
    Real abs_d_best(mal::Abs(d[0]));
    Real abs_d1(mal::Abs(d[1]));
    if( abs_d1 > abs_d_best ) { best_axis = 1; abs_d_best = abs_d1; }
    if( mal::Abs(d[2]) > abs_d_best ) best_axis = 2;
    Interval interval_p_d = Interval(p0[best_axis]).Merge(p1[best_axis]).Merge(p2[best_axis]);
    Interval interval_q_d = Interval(q0[best_axis]).Merge(q1[best_axis]).Merge(q2[best_axis]);
    if( !interval_p_d.TestOverlap(interval_q_d) ) return false;
    */

    // No early out left, perform 6 RC
    uint32 num_hits(0);
    Vec3 point;
    Real lambda;
    // Tri0 edges vs Tri1
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT_EXPLICIT_EDGES_SIMD( p0, p01, q0, q01, q02, nq, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 0, 0, lambda, point ); }
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT_EXPLICIT_EDGES_SIMD( p1, p2-p1, q0, q01, q02, nq, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 0, 1, lambda, point ); }
    if( num_hits == 2 ) return true;
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT_EXPLICIT_EDGES_SIMD( p2, -p02, q0, q01, q02, nq, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 0, 2, lambda, point ); }
    if( num_hits == 2 ) return true;
    // Tri1 edges vs Tri0
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT_EXPLICIT_EDGES_SIMD( q0, q01, p0, p01, p02, np, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 1, 0, lambda, point ); }
    if( num_hits == 2 ) return true;
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT_EXPLICIT_EDGES_SIMD( q1, q2-q1, p0, p01, p02, np, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 1, 1, lambda, point ); }
    if( num_hits == 2 ) return true;
    if( RayCast_Triangle3_DoubleSided_FAST_SEGMENT_EXPLICIT_EDGES_SIMD( q2, -q02, p0, p01, p02, np, point, lambda, p_context ) ) { ittr.SetE( num_hits++, 1, 2, lambda, point ); }
    if( num_hits == 2 ) return true;

    // GEO_LOG_ERROR_IF( num_hits != 2, "early outs should have discarded this intersection!" );

    GEO_LOG_ERROR_IF( num_hits != 0 && num_hits != 2, "geo::np::Intersection_Triangle3_Triangle3() num_hits = %u", num_hits );
    return num_hits > 0;
}
#endif //__ENABLE_SIMD_INTERSECTION_TRI_TRI

}} //namespace geo::np
