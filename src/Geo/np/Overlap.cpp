#include "Overlap.h"

/*\todo Include implementations actually used
#include "Overlap_X_Y_Analytic.h"
#include "Overlap_X_Y_BruteForce.h"
*/
#include <Geo/shape/GPolygonalShape.h> //TEMP while proper Overlap_X_Y_Z is not available
#include <Geo/shape/TriSurfaceShape3.h> //TEMP while proper Overlap_X_Y_Z is not available
#include <Mal/GRandom.h>
#include <Geo/np/RayCast.h>
#include <Geo/np/stats.h>

namespace geo {
namespace np {

//----------------------------------------------------------------
// Point Vs X
//----------------------------------------------------------------
bool Overlap_Point2_Triangle2( const Vec2& point,
                               const Vec2& tri0, const Vec2& tri1, const Vec2& tri2,
                               OverlapCache2* p_oc,
                               const Context* p_context )
{
    GEO_NP_STAT_INC( overlap.m_Point2_Triangle2 );
    // Barycentric transform B using tri0 as the local origin
    Mat2x2 invB = mal::Inverse( mal::GMat2x2_From_Columns( tri1-tri0, tri2-tri0 ) );
    Vec2 q = invB * (point - tri0);
    // Compute and test barycentric coords in [0..1] and \sum b_i = 1
    Vec3 b = Vec3( 1-q[0]-q[1], q[0], q[1] );
    for( int i=0; i<3; i++ )
        if( b[i] < Real(0) || b[i] > Real(1) ) //\note this makes dragon10K crash //if( b[i] < Real(-0.001f) || b[i] > Real(1.001f) )
            return false;
    return mal::ApproxEq( mal::Sum(b), Real(1) );//\todo consider non-default epsilon...
}

bool Overlap_Point2_Triangle2( const Vec2& point,
                               const Mat3x3& tri_invBs,
                               OverlapCache2* p_oc,
                               const Context* p_context )
{
    GEO_NP_STAT_INC( overlap.m_Point2_Triangle2 );
    // Compute and test barycentric coords in [0..1] and \sum b_i = 1
    Vec3 b( tri_invBs * mal::Concat(1,point) );
    for( int i=0; i<3; i++ )
        if( b[i] < Real(0) || b[i] > Real(1) ) //\note this makes dragon10K crash //if( b[i] < Real(-0.001f) || b[i] > Real(1.001f) )
            return false;
    return mal::ApproxEq( mal::Sum(b), Real(1) );//\todo consider non-default epsilon...
}

/* Counts number of hits with an infinite ray, odd if interior, even if exterior.
   \todo Works but it's slow, should use some kind of BVH...
*/
bool Overlap_Point2_Polygonal2( const Vec2& point,
                                const PolygonalShape2* p_polygonal, const Transform2& polygonal_tr, const Real* p_polygonal_dof,
                                OverlapCache2* p_oc,
                                const Context* p_context )
{
    GEO_NP_STAT_INC( overlap.m_Point2_Polygonal2 );
    if( !p_polygonal->IsClosed() ) return false;
    const Vec2* vec_sdof( reinterpret_cast<const Vec2*>(p_polygonal_dof) );
    Vec2 point_local( polygonal_tr.Inverse()*point );
    unsigned int num_edges( p_polygonal->GetNumVertices() );
    unsigned int num_hits(0);
    bv::AABB2 aabb;
    p_polygonal->ComputeBVD( aabb, polygonal_tr, vec_sdof );
    Real ray_length( 3*mal::Norm( aabb.GetHalfSizes() ) );
    Ray2 ray( point_local, mal::RandomUnitVec<Real,2>(), Interval(0,ray_length), Real(0) );
    for( unsigned int it_edge=0; it_edge < num_edges; it_edge++ )
    {
        Vec2f p0( p_polygonal->V_Pos( it_edge, vec_sdof ) );
        int next_vid( (it_edge+1) % p_polygonal->GetNumVertices() );
        Vec2f p1( p_polygonal->V_Pos( next_vid, vec_sdof ) );
        RayHit2 tmp_rh;
        if( RayCast_Segment2_SingleSided( ray.m_Pos, ray.m_Dir, ray.m_Interval, p0, p1, tmp_rh ) )
            num_hits++;
    }
    return num_hits % 2 == 1;
}

/* Counts number of hits with an infinite ray, odd if interior, even if exterior.
   \todo Works but it's slow, should use some kind of BVH...
*/
bool Overlap_Point3_TriSurface3( const Vec3& point,
                                 const TriSurfaceShape3* p_surface, const Transform3& surface_tr, const Real* p_surface_dof,
                                 OverlapCache3* p_oc,
                                 const Context* p_context )
{
    GEO_NP_STAT_INC( overlap.m_Point3_TriSurface3 );
    if( !p_surface->IsClosed() ) return false;
    const Vec3* vec_sdof( reinterpret_cast<const Vec3*>(p_surface_dof) );
    Vec3 point_local( surface_tr.Inverse()*point );
    unsigned int num_hits(0);
    // Compute ray length = longest AABB diagonal
    Real ray_length( mal::Infinity<Real>() );
#ifdef __GEO_TRISS_ENABLE_BVH
    const BVH_TriSurfaceShape3* pBVH = p_surface->GetBVH(); //\todo ASSUME UP-TO-DATE, but IN GLOBAL COORDS?!
    if( pBVH ) // Use BVH if exists
        ray_length = 3*mal::Norm( pBVH->GetRoot().m_Geometry.m_BV.GetHalfSizes() );
    else //fall-back to AABB otherwise
#endif
    {
        bv::AABB3 aabb;
        p_surface->ComputeBVD( aabb, surface_tr, vec_sdof );
        ray_length = 2*mal::Norm( aabb.GetHalfSizes() );
    }
    Ray3 ray( point_local, mal::RandomUnitVec<Real,3>(), Interval(0,ray_length), Real(0) );
    for( unsigned int it_tri=0; it_tri < p_surface->GetNumT(); it_tri++ )
    {
        Vec3f p0( p_surface->V_Pos( p_surface->T_VID(it_tri,0), vec_sdof ) );
        Vec3f p1( p_surface->V_Pos( p_surface->T_VID(it_tri,1), vec_sdof ) );
        Vec3f p2( p_surface->V_Pos( p_surface->T_VID(it_tri,2), vec_sdof ) );
        RayHit3 tmp_rh;
        if( RayCast_Triangle3_DoubleSided( ray.m_Pos, ray.m_Dir, ray.m_Interval, p0, p1, p2, tmp_rh ) )
            num_hits++;
    }
    return num_hits % 2 == 1;
}

bool Overlap_Point3_Tetrahedron3( const Vec3& point,
                                  const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3,
                                  OverlapCache3* p_oc,
                                  const Context* p_context )
{
    GEO_NP_STAT_INC( overlap.m_Point3_Tetrahedron3 );
    // Barycentric transform B using p0 as the local origin
    Mat3x3 invB = mal::Inverse( mal::GMat3x3_From_Columns( p1-p0, p2-p0, p3-p0 ) );
    Vec3 q = invB * (point - p0);
    // Compute and test barycentric coords in [0..1] and \sum b_i = 1
    Vec4 b = Vec4( 1-q[0]-q[1]-q[2], q[0], q[1], q[2] );
    for( int i=0; i<4; i++ )
        if( b[i] < Real(0) || b[i] > Real(1) ) //\note this makes dragon10K crash //if( b[i] < Real(-0.001f) || b[i] > Real(1.001f) )
            return false;
    return mal::ApproxEq( mal::Sum(b), Real(1) );//\todo consider non-default epsilon...
}

//----------------------------------------------------------------
// Segment Vs X
//----------------------------------------------------------------
bool Overlap_Segment2_Triangle2( const Vec2& s0, const Vec2& s1,
                                 const Vec2& tri0, const Vec2& tri1, const Vec2& tri2,
                                 OverlapCache2* p_oc,
                                 const Context* p_context )
{
    GEO_NP_STAT_INC( overlap.m_Segment2_Triangle2 );
    //-- Test 2 segment points inside triangle \todo Could optimize as Overlap_Point2_Triangle2 recomputes invB on each call!
    if( Overlap_Point2_Triangle2( s0, tri0, tri1, tri2, p_oc, p_context ) ) return true;
    if( Overlap_Point2_Triangle2( s1, tri0, tri1, tri2, p_oc, p_context ) ) return true;
    //-- Test 3 triangle edges against segment \todo Could optimize, I guess
    RayHit2 rh;
    if( RayCast_Segment2_SingleSided( s0, s1-s0, Interval(0,1), tri0, tri1, rh ) ) return true;
    if( RayCast_Segment2_SingleSided( s0, s1-s0, Interval(0,1), tri1, tri2, rh ) ) return true;
    if( RayCast_Segment2_SingleSided( s0, s1-s0, Interval(0,1), tri2, tri0, rh ) ) return true;
    return false;
}
bool Overlap_Segment2_Triangle2( const Vec2& s0, const Vec2& s1,
                                 const Mat3x3& tri_invBs, //triangle inverse barycentric transform
                                 OverlapCache2* p_oc,
                                 const Context* p_context )
{
    GEO_NP_STAT_INC( overlap.m_Segment2_Triangle2 );
    //-- Test 2 segment points inside triangle
    // Compute and test barycentric coords in [0..1] and \sum b_i = 1
    Vec3 b0( tri_invBs * mal::Concat(1,s0) );
    bool bOverlap0(true);
    for( int i=0; i<3; i++ )
        if( b0[i] < Real(0) || b0[i] > Real(1) ) //\note this makes dragon10K crash //if( b[i] < Real(-0.001f) || b[i] > Real(1.001f) )
            bOverlap0 = false;
    if( bOverlap0 && mal::ApproxEq( mal::Sum(b0), Real(1) )  ) return true; //\todo consider non-default epsilon...
    Vec3 b1( tri_invBs * mal::Concat(1,s1) );
    bool bOverlap1(true);
    for( int i=0; i<3; i++ )
        if( b1[i] < Real(0) || b1[i] > Real(1) ) //\note this makes dragon10K crash //if( b[i] < Real(-0.001f) || b[i] > Real(1.001f) )
            bOverlap1 = false;
    if( bOverlap1 && mal::ApproxEq( mal::Sum(b1), Real(1) )  ) return true; //\todo consider non-default epsilon...

    //-- Test barycentric segment against 3 barycentric triangle edges
    Vec3 v( b1-b0 );
    for( int i=0; i<3; i++ )
        if( b0[i] > 0 && v[i] < 0 )
        {
            int j( (i+1)%3 );
            Real lambda( -b0[i] / v[i] );
            Real bj( b0[j] + lambda*v[j] );
            if( bj >= 0 && bj <= 1 ) //&& b[k] >= 0 && b[k] <= 1 ) //\note b[i]=0, b[j]=1-b[k], no need to check 0<=b[k]<=1
                return true;
            /* b[i]=0, b[j]=1-b[k], no need to compute nor check 0<=b[k]<=1
            int k( (k+2)%3 );
            Vec3 b( b0 + lambda*v );
            if( b[j] >= 0 && b[j] <= 1 ) //&& b[k] >= 0 && b[k] <= 1 )
                return true;
            */
        }
    return false;
}

/* Similar to GRayCast_CenteredAABB but using BSlabs
   Clip ray_interval vs bdop intervals (bslabs in bcoords) and, if non-empty, overlap found
*/
bool RayCast_Triangle2_BDOP_UNFINISHED( const Vec2& ray_pos, const Vec2& ray_dir, const Interval& interval,
                                        const Mat3x3& tri_invBs, const bv::BDOP3& bdop,
                                        RayHit2& rh, RayCache2* p_rc,
                                        const Context* p_context )
{
    Vec3 b0( tri_invBs * mal::Concat(1,ray_pos+interval.Min()*ray_dir) );
    Vec3 b1( tri_invBs * mal::Concat(1,ray_pos+interval.Max()*ray_dir) );
    Vec3 dir_b( b1-b0 );
    Interval interval_b(0,1);
    // int first_hit_axis(0);
    for( int i=0; i<3; i++ )
    {
        if( mal::Abs(dir_b[i]) < p_context->m_Epsilon_Dir )
        {
            if( b0[i] < bdop[i].Min() || b0[i] > bdop[i].Max() )
            {
                //Empty interval and return false
                rh.m_Interval = Interval::Empty();
                return false;
            }
            // Otherwise, current axis does NOT clip the interval, and
            // other axis must be checked as usual.
        }
        else
        {
            Real inv_dir_b_i( mal::Rcp(dir_b[i]) );
            Real lambda0( (bdop[i].Min()-b0[i]) * inv_dir_b_i );
            Real lambda1( (bdop[i].Max()-b0[i]) * inv_dir_b_i );
            if( lambda1 < lambda0 )
            {
                Real tmp( lambda0 );
                lambda0 = lambda1;
                lambda1 = tmp;
            }
            // Clip lambda-interval and update first-axis
            if( interval_b.Min() < lambda0 ) { interval_b.Min() = lambda0; /*first_hit_axis = i;*/ }
            if( interval_b.Max() > lambda1 ) interval_b.Max() = lambda1;
            if( interval_b.IsEmpty() ) return false;
        }
    }
    if( interval_b.IsEmpty() )
        return false;
    else
    {
        //\todo Compute all rh stuff in GLOBAL COORDS
        // rh.m_Interval = !!!interval_b; //\todo COULD USE RATIO interval_b RATIOS to rescale original interval min/max
        // rh.m_Point  = ray_pos + rh.m_Interval.Min() * ray_dir;
        //\todo Compute proper normal in global coords FROM first hit axis
        // rh.m_Normal = Vec2::Zero();
        // rh.m_Normal[ first_hit_axis ] = ( b0[first_hit_axis] < Real(0) ) ? Real(-1) : Real(1);
        rh.m_FeatureId = feature_id(); //\todo Consider reporting V/E/F sub-features
        return true;
    }
}

bool Overlap_Segment2_Triangle2_BDOP( const Vec2& s0, const Vec2& s1,
                                      const Mat3x3& tri_invBs,
                                      const bv::BDOP3& bdop,
                                      OverlapCache2* p_oc,
                                      const Context* p_context )
{
    GEO_NP_STAT_INC( overlap.m_Segment2_Triangle2 );
    RayHit2 rh;
    return RayCast_Triangle2_BDOP_UNFINISHED( s0, s1-s0, Interval(0,1), tri_invBs, bdop, rh, 0, p_context );
}

/* Analogous to Overlap_Segment2_Triangle2_BDOP(), but with specific
   code instead of using "generic"
   RayCast_Tetrahedron3_BDOP_UNFINISHED()
*/
bool Overlap_Segment3_Tetrahedron3_BDOP( const Vec3& s0, const Vec3& s1,
                                         const Mat4x4& tet_invBs, //tetrahedron inverse barycentric transform
                                         const bv::BDOP4& bdop,
                                         OverlapCache3* p_oc,
                                         const Context* p_context )
{
    Vec3 ray_pos(s0);
    Vec3 ray_dir(s1-s0);
    Interval interval(0,1);
    Vec4 b0( tet_invBs * mal::Concat(1,ray_pos+interval.Min()*ray_dir) );
    Vec4 b1( tet_invBs * mal::Concat(1,ray_pos+interval.Max()*ray_dir) );
    Vec4 dir_b( b1-b0 );
    Interval interval_b(0,1);
    for( int i=0; i<4; i++ )
    {
        if( mal::Abs(dir_b[i]) < p_context->m_Epsilon_Dir )
        {
            if( b0[i] < bdop[i].Min() || b0[i] > bdop[i].Max() ) return false;
            // Otherwise, current axis does NOT clip the interval, and
            // other axis must be checked as usual.
        }
        else
        {
            Real inv_dir_b_i( mal::Rcp(dir_b[i]) );
            Real lambda0( (bdop[i].Min()-b0[i]) * inv_dir_b_i );
            Real lambda1( (bdop[i].Max()-b0[i]) * inv_dir_b_i );
            if( lambda1 < lambda0 )
            {
                Real tmp( lambda0 );
                lambda0 = lambda1;
                lambda1 = tmp;
            }
            // Clip lambda-interval and update first-axis
            if( interval_b.Min() < lambda0 ) interval_b.Min() = lambda0;
            if( interval_b.Max() > lambda1 ) interval_b.Max() = lambda1;
            if( interval_b.IsEmpty() ) return false;
        }
    }
    return true; //if we arrived here, no B-axis discards, therefore overlap
}

//----------------------------------------------------------------
// Triangle Vs X
//----------------------------------------------------------------

/* Straighforward test, suboptimal, \todo see RTCD Section 5.2.10, pg 172 for alternatives
 */
bool Overlap_Triangle3_Triangle3( const Vec3& p0, const Vec3& p1, const Vec3& p2,
                                  const Vec3& q0, const Vec3& q1, const Vec3& q2,
                                  OverlapCache3* p_oc,
                                  const Context* p_context )
{
    GEO_NP_STAT_INC( overlap.m_Triangle3_Triangle3 );
    /* \todo EARLY OUT from Intersection_Triangle3_Triangle3 should
       improve aggregate cost of lots of Overlap_Triangle3_Triangle3()
       calls:

       //Compute triangle non-unit normals Early-out with triangle planes
       const Vec3 nq( mal::Cross(q1-q0,q2-q0) );
       Real vec_dot_p_nq[3] = { mal::Dot(p0-q0,nq), mal::Dot(p1-q0,nq), mal::Dot(p2-q0,nq) };
       if( vec_dot_p_nq[0] < 0 && vec_dot_p_nq[1] < 0 && vec_dot_p_nq[2] < 0 ) return false;
       else if( vec_dot_p_nq[0] > 0 && vec_dot_p_nq[1] > 0 && vec_dot_p_nq[2] > 0 ) return false;
       const Vec3 np( mal::Cross(p1-p0,p2-p0) );
       Real vec_dot_q_np[3] = { mal::Dot(q0-p0,np), mal::Dot(q1-p0,np), mal::Dot(q2-p0,np) };
       if( vec_dot_q_np[0] < 0 && vec_dot_q_np[1] < 0 && vec_dot_q_np[2] < 0 ) return false;
       else if( vec_dot_q_np[0] > 0 && vec_dot_q_np[1] > 0 && vec_dot_q_np[2] > 0 ) return false;
    */

    RayHit3 rh;
    // p-edges vs q-tri
    if( RayCast_Triangle3_DoubleSided_COHERENT( p0, p1-p0, Interval(0,1), q0, q1, q2, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( p1, p2-p1, Interval(0,1), q0, q1, q2, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( p2, p0-p2, Interval(0,1), q0, q1, q2, rh ) ) return true;
    // q-edges vs p-tri
    if( RayCast_Triangle3_DoubleSided_COHERENT( q0, q1-q0, Interval(0,1), p0, p1, p2, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( q1, q2-q1, Interval(0,1), p0, p1, p2, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( q2, q0-q2, Interval(0,1), p0, p1, p2, rh ) ) return true;
    return false;
}

/* \todo OPTIMIZATION: Consider SAT test with 4+1 F axis and 3x6 EE
   axis = 23 axis... not trivial but probably faster than 6+12
   double-sided raycasts and 3 point-in-tetrahedron tests.
   \todo CONSIDER using RayCast_Triangle3_DoubleSided_COHERENT() from Clip.cpp...
*/
bool Overlap_Triangle3_Tetrahedron3( const Vec3& tri0, const Vec3& tri1, const Vec3& tri2,
                                     const Vec3& tet0, const Vec3& tet1, const Vec3& tet2, const Vec3& tet3,
                                     OverlapCache3* p_oc,
                                     const Context* p_context )
{
    GEO_NP_STAT_INC( overlap.m_Triangle3_Tetrahedron3 );

    // // Early-out with triangle plane
    // const Vec3 ntri( mal::Cross(tri1-tri0,tri2-tri0) );
    // Real vec_dot_tet_ntri[4] = { mal::Dot(tet0-tri0,ntri), mal::Dot(tet1-tri0,ntri), mal::Dot(tet2-tri0,ntri), mal::Dot(tet3-tri0,ntri) };
    // if( vec_dot_tet_ntri[0] < 0 && vec_dot_tet_ntri[1] < 0 && vec_dot_tet_ntri[2] < 0 && vec_dot_tet_ntri[3] < 0 ) return false;
    // else if( vec_dot_tet_ntri[0] > 0 && vec_dot_tet_ntri[1] > 0 && vec_dot_tet_ntri[2] > 0 && vec_dot_tet_ntri[3] > 0 ) return false;

#define __ENABLE_COHERENT
#ifdef __ENABLE_COHERENT //ensures same result as performing Tri2Tri for all 4 Tet faces, required by Fit_TetSolidShape3_To_TriSurfaceShape3()
    //-- Test 3 triangle points inside tetrahedron \todo Could optimize as Overlap_Point3_Tetrahedron3 recomputes invB on each call!
    if( Overlap_Point3_Tetrahedron3( tri0, tet0, tet1, tet2, tet3 ) ) return true;
    if( Overlap_Point3_Tetrahedron3( tri1, tet0, tet1, tet2, tet3 ) ) return true;
    if( Overlap_Point3_Tetrahedron3( tri2, tet0, tet1, tet2, tet3 ) ) return true;
    //-- Test tri against 4 tet triangles
    if( Overlap_Triangle3_Triangle3(tri0,tri1,tri2,tet0,tet1,tet2) ) return true;
    if( Overlap_Triangle3_Triangle3(tri0,tri1,tri2,tet0,tet2,tet3) ) return true;
    if( Overlap_Triangle3_Triangle3(tri0,tri1,tri2,tet0,tet3,tet1) ) return true;
    if( Overlap_Triangle3_Triangle3(tri0,tri1,tri2,tet1,tet2,tet3) ) return true;
    return false;
#else
    // Original code... should be correct, but caused divergence between Tri2Tet and Tri2Tri on specific Tet faces...
    //-- Test 3 triangle points inside tetrahedron \todo Could optimize as Overlap_Point3_Tetrahedron3 recomputes invB on each call!
    if( Overlap_Point3_Tetrahedron3( tri0, tet0, tet1, tet2, tet3 ) ) return true;
    if( Overlap_Point3_Tetrahedron3( tri1, tet0, tet1, tet2, tet3 ) ) return true;
    if( Overlap_Point3_Tetrahedron3( tri2, tet0, tet1, tet2, tet3 ) ) return true;
    RayHit3 rh;
    //-- Test 6 tetrahedron edges against the triangle \todo Could optimize, I guess...
    if( RayCast_Triangle3_DoubleSided_COHERENT( tet0, tet1-tet0, Interval(0,1), tri0, tri1, tri2, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( tet0, tet2-tet0, Interval(0,1), tri0, tri1, tri2, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( tet0, tet3-tet0, Interval(0,1), tri0, tri1, tri2, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( tet1, tet2-tet1, Interval(0,1), tri0, tri1, tri2, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( tet2, tet3-tet2, Interval(0,1), tri0, tri1, tri2, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( tet3, tet1-tet3, Interval(0,1), tri0, tri1, tri2, rh ) ) return true;
    //-- Test 3 triangle edges against 4 tetrahedron faces \todo Could optimize, I guess...
    // tri.e0
    if( RayCast_Triangle3_DoubleSided_COHERENT( tri0, tri1-tri0, Interval(0,1), tet0, tet1, tet2, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( tri0, tri1-tri0, Interval(0,1), tet0, tet2, tet3, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( tri0, tri1-tri0, Interval(0,1), tet0, tet3, tet1, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( tri0, tri1-tri0, Interval(0,1), tet1, tet2, tet3, rh ) ) return true;
    // tri.e1
    if( RayCast_Triangle3_DoubleSided_COHERENT( tri1, tri2-tri1, Interval(0,1), tet0, tet1, tet2, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( tri1, tri2-tri1, Interval(0,1), tet0, tet2, tet3, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( tri1, tri2-tri1, Interval(0,1), tet0, tet3, tet1, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( tri1, tri2-tri1, Interval(0,1), tet1, tet2, tet3, rh ) ) return true;
    // tri.e2
    if( RayCast_Triangle3_DoubleSided_COHERENT( tri2, tri0-tri2, Interval(0,1), tet0, tet1, tet2, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( tri2, tri0-tri2, Interval(0,1), tet0, tet2, tet3, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( tri2, tri0-tri2, Interval(0,1), tet0, tet3, tet1, rh ) ) return true;
    if( RayCast_Triangle3_DoubleSided_COHERENT( tri2, tri0-tri2, Interval(0,1), tet1, tet2, tet3, rh ) ) return true;
#endif
    return false;
}


bool Overlap_Triangle3_TriSurface3( const Vec3& tri0, const Vec3& tri1, const Vec3& tri2,
                                    const TriSurfaceShape3* p_surface, const Transform3& surface_tr, const Real* p_surface_dof,
                                    OverlapCache3* p_oc,
                                    const Context* p_context )
{
    GEO_NP_STAT_INC( overlap.m_Triangle3_TriSurface3 );
    Transform3 inv_tr( surface_tr.Inverse() );
    Vec3 vec_tri_pos_local[3] = { inv_tr*tri0, inv_tr*tri1, inv_tr*tri2 };
    const Vec3* vec_sdof( reinterpret_cast<const Vec3*>(p_surface_dof) );
#ifdef __GEO_TRISS_ENABLE_BVH
    const BVH_TriSurfaceShape3* pBVH = p_surface->GetBVH(); //\todo ASSUME UP-TO-DATE, but IN GLOBAL COORDS?!
    if( pBVH )
    {
        //\todo NOTE, triBV built in GLOBAL coords because BVH is expected there by now... this is a *MESS* I'll regret someday, I'm sure...
        auto triBV = BVH_TriSurfaceShape3::bv_type(tri0).Merge(tri1).Merge(tri2);
        triBV.Extend(0.01f);
        std::vector< BVH_TriSurfaceShape3::entry_index_type > vecOverlaps; //\todo alloc in scope allocator!
        if( pBVH->Test( triBV, vecOverlaps ) )
        {
            for( unsigned int it_o=0; it_o<vecOverlaps.size(); it_o++ )
            {
                uint32 tid = vecOverlaps[it_o];
                if( Overlap_Triangle3_Triangle3( vec_tri_pos_local[0], vec_tri_pos_local[1], vec_tri_pos_local[2],
                                                 p_surface->V_Pos( p_surface->T_VID(tid,0), vec_sdof ),
                                                 p_surface->V_Pos( p_surface->T_VID(tid,1), vec_sdof ),
                                                 p_surface->V_Pos( p_surface->T_VID(tid,2), vec_sdof ),
                                                 0, p_context ) )
                    return true;
            }
        }
    }
    else
    {
        for( unsigned int it_tri=0; it_tri < p_surface->GetNumT(); it_tri++ )
            if( Overlap_Triangle3_Triangle3( vec_tri_pos_local[0], vec_tri_pos_local[1], vec_tri_pos_local[2],
                                             p_surface->V_Pos( p_surface->T_VID(it_tri,0), vec_sdof ),
                                             p_surface->V_Pos( p_surface->T_VID(it_tri,1), vec_sdof ),
                                             p_surface->V_Pos( p_surface->T_VID(it_tri,2), vec_sdof ),
                                             0, p_context ) )
            return true;
    }
#else
    for( unsigned int it_tri=0; it_tri < p_surface->GetNumT(); it_tri++ )
        if( Overlap_Triangle3_Triangle3( vec_tri_pos_local[0], vec_tri_pos_local[1], vec_tri_pos_local[2],
                                         p_surface->V_Pos( p_surface->T_VID(it_tri,0), vec_sdof ),
                                         p_surface->V_Pos( p_surface->T_VID(it_tri,1), vec_sdof ),
                                         p_surface->V_Pos( p_surface->T_VID(it_tri,2), vec_sdof ),
                                         0, p_context ) )
            return true;
#endif
    return false;
}

}} //namespace geo::np
