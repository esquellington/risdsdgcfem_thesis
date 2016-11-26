#include "Contact_X_DCR3_BVH.h"
#include <Geo/shape/TetSolidShape3.h>
#include "Intersection.h"
#include "RayCast.h"
#include "Overlap.h"
#include "Clip.h"
#include <Geo/bv/GBDOP.h>
#include "../util/GSimpleSpatialHash.h"
#include <Mal/GRandom.h>

#ifdef __GEO_ENABLE_NP_SCRATCHPAD
#  include <util/ScopedAllocator.h>
#endif

#include <boost/bind.hpp>
#include <unordered_set>
#include <unordered_map>

#define __ENABLE_DCR3_BDT_EXTRA_AXIS
#ifdef __ENABLE_DCR3_BDT_EXTRA_AXIS
#  define __ENABLE_DCR3_BDT_EXTRA_AXIS_ELEMENT
//#  define __ENABLE_DCR3_BDT_EXTRA_AXIS_PATCH //\todo Doing it per-patch is slightly slower
#  define __ENABLE_DCR3_BDT_EXTRA_AXIS_DOP3K14 //Use all 4 extra DOP3_K14 axis (cube diags) AFTER cube faces (AABB3)
#endif

#define __ENABLE_DCR3_ELEMENT_CACHE //\todo Saves around 1ms out of 20ms for Armadillo100K vs Dragon100K and yields much cleaner code, enable by default once sufficiently tested
#ifdef __ENABLE_DCR3_ELEMENT_CACHE
#  define __ENABLE_DCR3_ELEMENT_CACHE
#  ifdef __GEO_ENABLE_NP_SIMD
#    define __ENABLE_SIMD_PRETRANSFORM
#    ifdef __ENABLE_SIMD_PRETRANSFORM
#      include "SIMD.h"
//#      define __ENABLE_SIMD_PRETRANSFORM_DOT4
#      define __ENABLE_SIMD_PRETRANSFORM_DOT3
#    endif
#  endif
#endif

namespace geo {
namespace bv {

template <unsigned D>
inline GAABB<D> GAABB_From( const GSphere<D>& sphere ) { return GAABB<D>().SetPosHalfSizes( sphere.GetPos(), GAABB<D>::vec_type(sphere.GetRadius()) ); }
template <unsigned D>
inline GAABB<D> GAABB_From( const GAABB<D>& aabb ) { return aabb; }
template <unsigned D, unsigned K>
inline GAABB<D> GAABB_From( const GKDOP<D,K>& kdop )
{
    GAABB<D> aabb;
    for( unsigned int i=0; i<D; i++ )
    {
        aabb.m_Pos[i] = kdop.GetInterval(i).Mid();
        aabb.m_HalfSizes[i] = Real(0.5)*kdop.GetInterval(i).Length();
    }
    return aabb;
}

template <unsigned D>
inline GAABB<D> GAABB_Intersect( const GAABB<D>& aabb1, const GAABB<D>& aabb2 )
{
    GAABB<D> aabb;
    for( unsigned int i=0; i<D; i++ )
    {
        Interval interval = Interval( aabb1.m_Pos[i]-aabb1.m_HalfSizes[i], aabb1.m_Pos[i]+aabb1.m_HalfSizes[i] )
                            .Intersect( Interval( aabb2.m_Pos[i]-aabb2.m_HalfSizes[i], aabb2.m_Pos[i]+aabb2.m_HalfSizes[i] ) );
        aabb.m_Pos[i] = interval.Mid();
        aabb.m_HalfSizes[i] = Real(0.5)*interval.Length();
    }
    return aabb;
}

}}

namespace geo {
namespace np {

//! Transform a plane equation \todo Generalize to GTransformPlane<D> 2D/3D and move elsewhere
inline void TransformPlane( const Vec3& plane_n, Real plane_d,
                            const Transform3& tr,
                            Vec3& transformed_n, Real& transformed_d )
{
    Vec3 plane_p( -plane_d*plane_n );
    transformed_n = tr.m_Rot * plane_n;
    Vec3 transformed_p( tr * plane_p );
    transformed_d = -mal::Dot(transformed_n,transformed_p);
}

//--------------------------------------------------------------------------------------------------------------------------------
// TetSolidShape3/DCR Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

/*! TetSolidShape3/DCR Vs Plane */
bool Contact_DCRTS3_Plane3_BVH_TSS_Lazy_DCR( const TetSolidShape3* p_solid, const Transform3& solid_tr, const Real* vec_dof,
                                             const Vec3& plane_normal, Real plane_coeff_d,
                                             ContactData3& cd, ContactCache3* p_cc,
                                             const Context* p_context )
{
    bool bViz(p_context->m_DCR2DCR_Viz_Enabled);

    // GEO_LOG("Contact_DCRTS3_Plane3_BVH_TSS_Lazy_DCR()");
    cd.Begin();

    //TEMP: REQUIRE DCR
    const DCR_TetSolidShape3* pDCR( p_solid->GetDCR() );
    if( 0 == pDCR )
    {
        GEO_LOG_WARNING("Contact_DCRTS3_Plane3_BVH_TSS_Lazy_DCR() with no DCR, ignoring");
        return false;
    }

    //TEMP: REQUIRE BVH
    BVH_TetSolidShape3* pBVH( p_solid->GetBVH() );
    if( 0 == pBVH )
    {
        GEO_LOG_WARNING("Contact_DCRTS3_Plane3_BVH_TSS_Lazy_DCR() with no BVH, ignoring");
        return false;
    }

    //\todo THIS IS dangerous... Should have a generic SRV type with safe casts...
    const Vec3* vec_sdof( reinterpret_cast<const Vec3*>(vec_dof) );
    const Vec3* default_sdof( p_solid->GetVecDefaultSDOF() );

    /*\todo IF p_solid/DOF/Transform->m_TimeStamp != pBVH->GetTimeStamp()
      Perform lazy Refit() or ASSUME up-to-date?
    */

    /* Refit BVH
       TEMP: ASSUME UP TO DATE!
       \todo See Contact_DCRMS2_DCRMS2_BVH_MSS_Lazy_DCR() for detailed discussion
    switch( p_context->m_BVH_Method )
    {
    case Context::eBVHM_BV_E:
        pBVH->Refit( boost::bind<void>( &GEBV_TetSolidShape3_E<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
                                        p_solid, solid_tr, vec_sdof, _1, _2) );
        break;
    case Context::eBVHM_BV_BSlab:
        pBVH->Refit( boost::bind<void>( &GEBV_TetSolidShape3_BSlab<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
                                        p_solid, pDCR, solid_tr, vec_sdof, _1, _2) );
        break;
    case Context::eBVHM_BV_BDOP:
        pBVH->Refit( boost::bind<void>( &GEBV_TetSolidShape3_BDOP<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
                                        p_solid, pDCR, solid_tr, vec_sdof, _1, _2) );
        break;
    default: break;
    }
    */

    /* Test overlap
       \todo We should compute the actual set of DCR.E that overlap
       the halfspace. It's NOT a good a idea to use the BVH here, as
       the planeBV will be often infinite (if unbounded and not
       axis-aligned). By now, we just consider that all DCR.E overlap
       the plane.
    */
    std::vector< BVH_TetSolidShape3::entry_index_type > vecOverlaps;
    for( unsigned int it_e=0; it_e<pDCR->m_NumElements; it_e++ )
        vecOverlaps.push_back( it_e );

    // Transform plane to mesh-local coords
    Vec3 plane_normal_1;
    Real plane_coeff_d_1;
    Transform3 inv_solid_tr( mal::Inverse(solid_tr) );
    TransformPlane( plane_normal, plane_coeff_d, inv_solid_tr, plane_normal_1, plane_coeff_d_1 );

    // Process overlaps
    for( auto eid : vecOverlaps )
    {
        const DCR_TetSolidShape3::Element& ed( pDCR->m_vecE[eid] );
        // Get element nodes in mesh-local coords
        uint32 vec_nid[4] = { p_solid->T_VID(eid,0), p_solid->T_VID(eid,1), p_solid->T_VID(eid,2), p_solid->T_VID(eid,3) };
        // Test overlap DCR.E vs Plane in mesh-local coords
        // \todo THIS WOULD BE MUCH TIGHTER if we used BDOP/BSlab
        // vertices, and more DCR.E could be classified as interior
        // (much cheaper to handle than piercing DCR.E)
        Real vec_dist[4] = { mal::Dot( plane_normal_1, vec_sdof[vec_nid[0]] ) + plane_coeff_d_1,
                             mal::Dot( plane_normal_1, vec_sdof[vec_nid[1]] ) + plane_coeff_d_1,
                             mal::Dot( plane_normal_1, vec_sdof[vec_nid[2]] ) + plane_coeff_d_1,
                             mal::Dot( plane_normal_1, vec_sdof[vec_nid[3]] ) + plane_coeff_d_1 };
        if( vec_dist[0] <= 0 || vec_dist[1] <= 0 || vec_dist[2] <= 0 || vec_dist[3] <= 0 ) // interior | piercing
        {
            // Get barycentric transform
            Mat4x4 Bm( 1, 1, 1, 1,
                       default_sdof[vec_nid[0]].x(), default_sdof[vec_nid[1]].x(), default_sdof[vec_nid[2]].x(), default_sdof[vec_nid[3]].x(),
                       default_sdof[vec_nid[0]].y(), default_sdof[vec_nid[1]].y(), default_sdof[vec_nid[2]].y(), default_sdof[vec_nid[3]].y(),
                       default_sdof[vec_nid[0]].z(), default_sdof[vec_nid[1]].z(), default_sdof[vec_nid[2]].z(), default_sdof[vec_nid[3]].z() );
            Mat4x4 invBm( mal::Inverse( Bm ) ); //\todo could be precomputed!
            if( vec_dist[0] <= 0 && vec_dist[1] <= 0 && vec_dist[2] <= 0 && vec_dist[3] <= 0 ) //interior
            {
                /* Compute Single CP for each overlapping DCR.P
                   with barycentrically averaged pos/depth/normal from
                   DCR.T
                */
                for( unsigned int it_patch=0; it_patch < ed.m_NumPatches; it_patch++ )
                {
                    //\todo THIS IS VERY SLOW, AVOID per-triangle loop
                    //using precomputed per-patch data, do NOT compute
                    //O(T) total_area, not even O(1) total_area_sq, as
                    //it's NOT USED by now
                    Vec4 avg_p_b(0);
                    Vec3 avg_p(0);
                    Real total_area(0);
                    const geo::DCR_TetSolidShape3::Patch& pd( pDCR->m_vecP[ed.m_FirstPID + it_patch] );
                    for( unsigned int it_tid=0; it_tid < pd.m_NumTriangles; it_tid++ )
                    {
                        unsigned int tid( pd.m_FirstTID + it_tid );
                        uint32 vid0( pDCR->m_vecT[tid].GetVID(0) );
                        uint32 vid1( pDCR->m_vecT[tid].GetVID(1) );
                        uint32 vid2( pDCR->m_vecT[tid].GetVID(2) );
                        // Compute barycentric coords for p,q
                        Vec4 p_b( invBm * mal::Concat(1,pDCR->m_vecV[vid0]) ); //\todo could be precomputed!
                        Vec4 q_b( invBm * mal::Concat(1,pDCR->m_vecV[vid1]) ); //\todo could be precomputed!
                        Vec4 r_b( invBm * mal::Concat(1,pDCR->m_vecV[vid2]) ); //\todo could be precomputed!
                        // Compute DCR.V global coords from node positions \todo IMPORTANT: These are NOT global, we have NOT YET applied solid_tr
                        Vec3 p( p_b[0] * vec_sdof[vec_nid[0]] + p_b[1] * vec_sdof[vec_nid[1]] + p_b[2] * vec_sdof[vec_nid[2]] + p_b[3] * vec_sdof[vec_nid[3]] );
                        Vec3 q( q_b[0] * vec_sdof[vec_nid[0]] + q_b[1] * vec_sdof[vec_nid[1]] + q_b[2] * vec_sdof[vec_nid[2]] + q_b[3] * vec_sdof[vec_nid[3]] );
                        Vec3 r( r_b[0] * vec_sdof[vec_nid[0]] + r_b[1] * vec_sdof[vec_nid[1]] + r_b[2] * vec_sdof[vec_nid[2]] + r_b[3] * vec_sdof[vec_nid[3]] );
                        // Area-weighted average pos (\todo consider using precomputed pose areas)
                        Real area( Real(0.5) * mal::Norm(mal::Cross( q - p, r - p ) ) );
                        avg_p_b += area * Real(0.333333) * (p_b + q_b + r_b);
                        avg_p += area * Real(0.333333) * (p + q + r);
                        total_area += area;
                        //DEBUG
                        if( bViz )
                        {
                            // cd.m_VD.m_vecDCR3_Vs_Plane_vecInternalT.push_back( util::make_triad(solid_tr*p,solid_tr*q,solid_tr*r) );
                        }
                    }
                    if( total_area > 0 )
                    {
                        Real rcp_total_area( mal::Rcp(total_area) );
                        Vec3 avg_p_0( solid_tr * (avg_p * rcp_total_area) ); //globalize using solid_tr
                        avg_p_b = avg_p_b * rcp_total_area;
                        const Real cBaryEps(0.01f);
                        GEO_LOG_ERROR_IF( avg_p_b[0] < -cBaryEps || avg_p_b[0] > 1+cBaryEps
                                          || avg_p_b[1] < -cBaryEps || avg_p_b[1] > 1+cBaryEps
                                          || avg_p_b[2] < -cBaryEps || avg_p_b[2] > 1+cBaryEps
                                          || avg_p_b[3] < -cBaryEps || avg_p_b[3] > 1+cBaryEps,
                                          "Interior E[%u] AvgPos_b = (%f,%f,%f,%f)",
                                          eid, avg_p_b[0], avg_p_b[1], avg_p_b[2], avg_p_b[3] );
                        Real avg_dist( mal::Dot(plane_normal,avg_p_0) + plane_coeff_d );
                        cd.AddCP( avg_p_0, //on mesh
                                  avg_p_0 - avg_dist*plane_normal, //on plane
                                  plane_normal, -avg_dist,
                                  total_area );
                        cd.AddPOF( PointOnFeature( feature_id(eFT_Tetrahedron,eid),avg_p_b),
                                   PointOnFeature() );
                        GEO_NP_STAT_INC( mig2015.m_Num_CP );
                    }
                }
            }
            else //piercing
            {
                //\todo CONSIDER transforming plane into element-local coords (barycentric transform and plane are affine, SHOULD WORK)
                // Transform all DCR.E.V
                // \todo THIS MAY BE a uselss optimization, as we'll need the barycentric coords anyway (for piercing triangles, at least) and the global pos is easily computed from there...
                util::ScopedAllocator scoped_allocator( p_context->m_ScratchPad, "Contact_DCRTS3_Plane3_BVH_TSS_Lazy_DCR" );
                Vec3* vec_v1_0 = scoped_allocator.NewArrayPOD<Vec3>(ed.m_NumVertices);
                {
                    Mat3x3 Ds( mal::GMat3x3_From_Columns(vec_sdof[vec_nid[1]]-vec_sdof[vec_nid[0]],
                                                         vec_sdof[vec_nid[2]]-vec_sdof[vec_nid[0]],
                                                         vec_sdof[vec_nid[3]]-vec_sdof[vec_nid[0]]) );
                    Mat3x3 Dm( mal::GMat3x3_From_Columns(default_sdof[vec_nid[1]]-default_sdof[vec_nid[0]],
                                                         default_sdof[vec_nid[2]]-default_sdof[vec_nid[0]],
                                                         default_sdof[vec_nid[3]]-default_sdof[vec_nid[0]]) ); //\todo Could be precomputed in DCR::ED
                    Mat3x3 Ds_invDm( Ds * mal::Inverse(Dm) );
                    Transform3 tr_Ds_invDm( solid_tr.m_Pos, solid_tr.m_Rot * Ds_invDm );
                    Vec3 p0_0( solid_tr.m_Rot * vec_sdof[vec_nid[0]] );
                    // Transform ed1 m_NumVertices vertices from m_First and store them in 0-based indices inside vec_v1_0
                    for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
                        vec_v1_0[it_vie] = tr_Ds_invDm * (pDCR->m_vecV[ed.m_FirstVID + it_vie] - default_sdof[vec_nid[0]]) + p0_0;
                    GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed.m_NumVertices );
                }
                // Get inv barycentric transform (local to solid_tr) \todo TRY TO use Ds instead... 4x4 inverse is EXPENSIVE!
                Mat4x4 invBs( mal::Inverse( Mat4x4( 1, 1, 1, 1,
                                                    vec_sdof[vec_nid[0]].x(), vec_sdof[vec_nid[1]].x(), vec_sdof[vec_nid[2]].x(), vec_sdof[vec_nid[3]].x(),
                                                    vec_sdof[vec_nid[0]].y(), vec_sdof[vec_nid[1]].y(), vec_sdof[vec_nid[2]].y(), vec_sdof[vec_nid[3]].y(),
                                                    vec_sdof[vec_nid[0]].z(), vec_sdof[vec_nid[1]].z(), vec_sdof[vec_nid[2]].z(), vec_sdof[vec_nid[3]].z() ) ) );

                /* Compute Single CP for each piercing DCR.P with
                   barycentrically averaged pos/depth/normal from
                   DCR.T
                */
                for( unsigned int it_patch=0; it_patch < ed.m_NumPatches; it_patch++ )
                {
                    const geo::DCR_TetSolidShape3::Patch& pd( pDCR->m_vecP[ed.m_FirstPID + it_patch] );
                    Vec4 avg_p_b(0);
                    Vec3 avg_p_0(0);
                    Real total_area(0);
                    for( unsigned int it_tid=0; it_tid < pd.m_NumTriangles; it_tid++ )
                    {
                        unsigned int tid( pd.m_FirstTID + it_tid );
                        uint32 vid0( pDCR->m_vecT[tid].GetVID(0) );
                        uint32 vid1( pDCR->m_vecT[tid].GetVID(1) );
                        uint32 vid2( pDCR->m_vecT[tid].GetVID(2) );
                        Vec3 p_0( vec_v1_0[ vid0 - ed.m_FirstVID ] );
                        Vec3 q_0( vec_v1_0[ vid1 - ed.m_FirstVID ] );
                        Vec3 r_0( vec_v1_0[ vid2 - ed.m_FirstVID ] );
                        // Clip triangle vs plane and integrate over overlapping surface
                        Vec3 vec_clipped_points[4];
                        unsigned int num_clipped = Clip_Triangle3_Plane3( p_0, q_0, r_0, plane_normal, plane_coeff_d, vec_clipped_points, p_context );
                        if( num_clipped > 0 )
                        {
                            Vec3 cp0_0( vec_clipped_points[0] );
                            Real acc_area(0);
                            Vec3 acc_bary_0(0);
                            for( unsigned int i=1; i<num_clipped-1; i++ )
                            {
                                Vec3 cp1_0( vec_clipped_points[i] );
                                Vec3 cp2_0( vec_clipped_points[i+1] );
                                Real area( Real(0.5) * mal::Norm( mal::Cross( cp1_0 - cp0_0, cp2_0 - cp0_0 ) ) ); //scalar addition because poly IS guaranteed to be convex, otherwise should add vector-areas
                                acc_area += area;
                                Vec3 bary_0( Real(0.333333) * (cp0_0 + cp1_0 + cp2_0) );
                                acc_bary_0 += area * bary_0;
                            }
                            avg_p_0 += acc_bary_0;
                            total_area += acc_area;

                            /* \todo This was DISCONTINUOUS wrt non-piercing DCR.T barycenters and DCR.P CP "jumped" when num_clipped changes between 3 and 4
                            Vec3 acc_cp(0);
                            for( unsigned int i=0; i<num_clipped; i++ ) acc_cp += vec_clipped_points[i];
                            acc_cp *= acc_area / num_clipped;
                            avg_p_0 += acc_cp;
                            avg_p_b += invDs * mal::Concat(1,acc_cp);
                            */
                            // GEO_LOG("Added num_clipped=%u total_area=%f",num_clipped,total_area);

                            //TEMP DEBUG Clip_Triangle3_Plane3()
                            if( bViz )
                            {
                                if( num_clipped > 3 )
                                {
                                    for( unsigned int i=0; i<num_clipped; i++ )
                                    {
                                        /*
                                          Vec3 p(vec_clipped_points[i]);
                                          Vec3 q(vec_clipped_points[(i+1)%num_clipped]);
                                          cd.AddCP( p, q, mal::SafeNormalized(q-p), mal::Norm(q-p), total_area );
                                          cd.AddPOF( PointOnFeature(), PointOnFeature() );
                                        */
                                        cd.m_VD.m_vecDCR3_Vs_Primitive_vecClippedT_Segments.push_back( std::make_pair(vec_clipped_points[i],vec_clipped_points[(i+1)%num_clipped]) );
                                    }
                                }
                            }
                        }
                    }
                    if( total_area > 0 )
                    {
                        Real rcp_total_area( mal::Rcp(total_area) );
                        avg_p_0 = avg_p_0 * rcp_total_area;
                        avg_p_b = invBs * mal::Concat(1,inv_solid_tr * avg_p_0);
                        const Real cBaryEps(0.01f);
                        GEO_LOG_ERROR_IF( avg_p_b[0] < -cBaryEps || avg_p_b[0] > 1+cBaryEps
                                          || avg_p_b[1] < -cBaryEps || avg_p_b[1] > 1+cBaryEps
                                          || avg_p_b[2] < -cBaryEps || avg_p_b[2] > 1+cBaryEps
                                          || avg_p_b[3] < -cBaryEps || avg_p_b[3] > 1+cBaryEps,
                                          "Pierced E[%u] AvgPos_b = (%f,%f,%f,%f)",
                                          eid, avg_p_b[0], avg_p_b[1], avg_p_b[2], avg_p_b[3] );
                        Real avg_dist( mal::Dot(plane_normal,avg_p_0) + plane_coeff_d );
                        cd.AddCP( avg_p_0, //on mesh
                                  avg_p_0 - avg_dist*plane_normal, //on plane
                                  plane_normal, -avg_dist,
                                  total_area );
                        cd.AddPOF( PointOnFeature( feature_id(eFT_Tetrahedron,eid),avg_p_b),
                                   PointOnFeature() );
                        GEO_NP_STAT_INC( mig2015.m_Num_CP );
                    }
                }
            }
        }
        // else exterior
    }
    if( cd.Size() > 0 ) cd.SetNumDisjointManifolds( 1 ); //\todo Count it properly, for a deformable object it may be > 1
    return cd.End();
}

//\todo IMPORTANT: This CANNOT be the SAME class as 2D version because they get confused, put into anonymous namespace or use specific crossing_segment/triangle_pair_type name
struct crossing_triangle_pair_type
{
    struct endpoint_id_type
    {
        uint32 m_Bits;
        inline endpoint_id_type() : m_Bits(0xFFFFFFFF) {}
        inline endpoint_id_type( uint32 ctpid, int epictp ) : m_Bits(2*ctpid + epictp) { GEO_ASSERT(epictp == 0 || epictp == 1); }
        inline bool IsValid() const { return m_Bits != 0xFFFFFFFF; }
        inline uint32 CTPID() const { return m_Bits/2; }
        inline int EPICTP() const { return m_Bits%2; }
    };
    intersection_triangle3_triangle3_result_type m_ITTR;
    uint32 m_EID[2]; //EID on each object
    uint32 m_PID[2]; //PID on each object
    uint32 m_TID[2]; //TID on each object
    // Triangle points (CCW)
    Vec3 m_vecPos_0[2][3]; //[object 0..1][vertex 0..3]
    // IC data
    uint32 m_NumConnected[2]; //!< CTP segment endpoint connection count (ideally 1)
    endpoint_id_type m_NeighbourEPID[2];
    uint32 m_ParentId; //!< MF-Set parent for IC computation
    uint32 m_CurveId; //!< MF-Set parent for IC computation
    inline crossing_triangle_pair_type( const intersection_triangle3_triangle3_result_type& ittr,
                                        uint32 eid1, uint32 pid1, uint32 tid1, const Vec3& a1_0, const Vec3& b1_0, const Vec3& c1_0,
                                        uint32 eid2, uint32 pid2, uint32 tid2, const Vec3& a2_0, const Vec3& b2_0, const Vec3& c2_0 )
    : m_ITTR(ittr)
    , m_ParentId(0xFFFFFFFF), m_CurveId(0xFFFFFFFF)
        {
            m_EID[0] = eid1; m_EID[1] = eid2;
            m_PID[0] = pid1; m_PID[1] = pid2;
            m_TID[0] = tid1; m_TID[1] = tid2;
            m_vecPos_0[0][0] = a1_0;
            m_vecPos_0[0][1] = b1_0;
            m_vecPos_0[0][2] = c1_0;
            m_vecPos_0[1][0] = a2_0;
            m_vecPos_0[1][1] = b2_0;
            m_vecPos_0[1][2] = c2_0;
            m_NumConnected[0] = 0;
            m_NumConnected[1] = 0;
            m_NeighbourEPID[0] = endpoint_id_type();
            m_NeighbourEPID[1] = endpoint_id_type();
        }
};

inline uint32 CTP_ComputeRootId( std::vector<crossing_triangle_pair_type>& vec_ctp, uint32 ctpid )
{
    if( vec_ctp[ctpid].m_ParentId != ctpid ) vec_ctp[ctpid].m_ParentId = CTP_ComputeRootId(vec_ctp,vec_ctp[ctpid].m_ParentId);
    return vec_ctp[ctpid].m_ParentId;
}

/* Check adjacency of ctp/octp around edge E0 on ctp that generates it_point_in_ctp_segment.
   Adjacency requires:
   - Finding the "other" triangle T0 across E0
   - That the triangle T1 on the "other object" is the same for both ctp/octp
     - This helps discard matching among the 4 E-E ctp that result
       from 2 pairs of triangles intersecting their their shared
       edges.  The 4 ctp may be VERY close/coincident, but the
       ctp/octp match must bet topologically compatible.
   \todo I'm sure this code can be greatly simplified, it's too "symmetric" not to
*/
inline bool CTP_IsAdjacent( const crossing_triangle_pair_type& ctp, uint32 it_point_in_ctp_segment,
                            const crossing_triangle_pair_type& octp,
                            const DCR_TetSolidShape3* pDCR1, const DCR_TetSolidShape3* pDCR2 )
{
    bool bAdjacent(false);
    if( ctp.m_ITTR.IsVF() )
    {
        if( it_point_in_ctp_segment == 0 )
            bAdjacent = pDCR1->m_vecT[ctp.m_TID[0]].GetNTID(ctp.m_ITTR.m_EdgeInT0) == octp.m_TID[0] && ctp.m_TID[1] == octp.m_TID[1];
        else
            bAdjacent = pDCR1->m_vecT[ctp.m_TID[0]].GetNTID(ctp.m_ITTR.m_EdgeInT1) == octp.m_TID[0] && ctp.m_TID[1] == octp.m_TID[1];
    }
    else if( ctp.m_ITTR.IsFV() )
    {
        if( it_point_in_ctp_segment == 0 )
            bAdjacent = pDCR2->m_vecT[ctp.m_TID[1]].GetNTID(ctp.m_ITTR.m_EdgeInT0) == octp.m_TID[1] && ctp.m_TID[0] == octp.m_TID[0];
        else
            bAdjacent = pDCR2->m_vecT[ctp.m_TID[1]].GetNTID(ctp.m_ITTR.m_EdgeInT1) == octp.m_TID[1] && ctp.m_TID[0] == octp.m_TID[0];
    }
    else if( ctp.m_ITTR.IsEE() )
    {
        // GEO_ASSERT( ctp.m_ITTR.m_T0 == 0 && ctp.m_ITTR.m_T1 == 1 );
        if( it_point_in_ctp_segment == 0 )
            bAdjacent = pDCR1->m_vecT[ctp.m_TID[0]].GetNTID(ctp.m_ITTR.m_EdgeInT0) == octp.m_TID[0] && ctp.m_TID[1] == octp.m_TID[1];
        else
            bAdjacent = pDCR2->m_vecT[ctp.m_TID[1]].GetNTID(ctp.m_ITTR.m_EdgeInT1) == octp.m_TID[1] && ctp.m_TID[0] == octp.m_TID[0];
    }
    return bAdjacent;
}

inline bool CTP_AreAdjacent( const crossing_triangle_pair_type& ctp, uint32 it_point_in_ctp_segment,
                             const crossing_triangle_pair_type& octp, uint32 it_point_in_octp_segment,
                             const DCR_TetSolidShape3* pDCR1, const DCR_TetSolidShape3* pDCR2 )
{
    return CTP_IsAdjacent( ctp, it_point_in_ctp_segment, octp, pDCR1, pDCR2 )
        && CTP_IsAdjacent( octp, it_point_in_octp_segment, ctp, pDCR1, pDCR2 ); //\todo THIS discards adjacencies and therefore IS REQUIRED...
}

struct dcr3_intersection_curve_type
{
    std::vector< uint32 > m_vecCTPID;
    uint32 m_IBID[2]; //IB_1, IB_2
    inline dcr3_intersection_curve_type() { m_IBID[0] = 0xFFFFFFFF; m_IBID[1] = 0xFFFFFFFF; }
};

class dcr3_intersection_boundary_type
{
public:
    // DCR.Patch aggregate data
    class Patch
    {
    public:
        explicit inline Patch( uint32 pid, uint32 parent_ibid ) : m_PID(pid), m_ParentIBID(parent_ibid), m_NumTriangles(0), m_AvgPos_0(0), m_VectorArea_0(0), m_AreaSq(0), m_Area(0) {}
    public:
        uint32 m_PID; uint32 m_ParentIBID; uint32 m_NumTriangles; Vec3 m_AvgPos_0; Vec3 m_VectorArea_0; Real m_AreaSq; Real m_Area;
    };

public:
    explicit inline dcr3_intersection_boundary_type( uint32 icid )
    : m_ParentId(0xFFFFFFFF), m_NumTriangles(0), m_AvgPos_0(0), m_VectorArea_0(0), m_AreaSq(0), m_Area(0) { m_vecICID.push_back(icid); }

    void Merge( dcr3_intersection_boundary_type& other )
        {
            // Add IC
            for( auto icid : other.m_vecICID ) m_vecICID.push_back(icid);
            // Add or Merge crossing patches
            for( auto opatch : other.m_mapCrossingP )
            {
                auto it_patch( m_mapCrossingP.find(opatch.second.m_PID) );
                if( it_patch == m_mapCrossingP.end() ) //Add new patch with updated parent
                {
                    opatch.second.m_ParentIBID = m_ParentId;
                    m_mapCrossingP.insert( opatch );
                }
                else //Aggregate patch
                {
                    it_patch->second.m_NumTriangles += opatch.second.m_NumTriangles;
                    it_patch->second.m_AvgPos_0 += opatch.second.m_AvgPos_0;
                    it_patch->second.m_VectorArea_0 += opatch.second.m_VectorArea_0;
                    it_patch->second.m_AreaSq += opatch.second.m_AreaSq;
                    it_patch->second.m_Area += opatch.second.m_Area;
                }
            }
            // Add patches setting parent to root
            for( const auto& patch : other.m_vecP )
            {
                m_vecP.push_back(patch);
                m_vecP.back().m_ParentIBID = m_ParentId;
            }
            // Merge aggregate data
            m_NumTriangles += other.m_NumTriangles;
            m_AvgPos_0 += other.m_AvgPos_0;
            m_VectorArea_0 += other.m_VectorArea_0;
            m_AreaSq += other.m_AreaSq;
            m_Area += other.m_Area;
            // Clear other
            other.m_vecICID.clear();
            other.m_mapCrossingP.clear();
            other.m_vecP.clear();
        }

public:
    uint32 m_ParentId;
    std::vector<uint32> m_vecICID; //list of IC that bound the IB
    std::unordered_map<uint32,Patch> m_mapCrossingP; //crossing patches, may contribute to other IB
    // std::unordered_map<uint32,Patch> m_mapInternalP; //fully interior patches \todo NO, this map is GLOBAL to detect IB merging
    std::vector<Patch> m_vecP; //This will be filled at the end (*AFTER* IB merging) and contain exclusive internal patches first and potentially non-exclusive crossing patches
    // IB aggregate data
    uint32 m_NumTriangles;
    Vec3 m_AvgPos_0;
    Vec3 m_VectorArea_0;
    Real m_AreaSq;
    Real m_Area;
    // std::vector< util::triad<Vec3,Vec3,Vec3> > m_vecT;
};

inline uint32 IB_ComputeRootId( std::vector<dcr3_intersection_boundary_type>& vec_ib, uint32 ibid )
{
    if( vec_ib[ibid].m_ParentId != ibid ) vec_ib[ibid].m_ParentId = IB_ComputeRootId( vec_ib, vec_ib[ibid].m_ParentId );
    return vec_ib[ibid].m_ParentId;
}

/* Per-element cache
   - Stores per-element which used multiple times during a Contact() test and is constant and too expensive to recompute
   - Will be filled on demand
   \todo CONSIDER scoped allocation on init (vecECD)
   IMPORTANT: Magnitures are stored RELATIVE to TetSolidShape3 transform, NOT GLOBAL!
*/
class dcr3_element_cache_type
{
public:
    class ECD //ElementCachedData
    {
    public:
        inline ECD() : m_IsValid(false), m_FirstTransformedVID(0xFFFFFFFF), m_Stats_NumTransform(0) {}
        inline bool IsValid() const { return m_IsValid; }
        inline const Mat4x4& Bs() const { GEO_ASSERT(m_IsValid); return m_Bs; }
        inline const Mat4x4& invBs() const { GEO_ASSERT(m_IsValid); return m_invBs; }
        inline const Mat4x4& invBm() const { GEO_ASSERT(m_IsValid); return m_invBm; }
        inline const Mat4x4& Bs_invBm() const { GEO_ASSERT(m_IsValid); return m_Bs_invBm; }
    private:
        friend class dcr3_element_cache_type;
        inline void Update( uint32 eid, const TetSolidShape3* p_solid, const DCR_TetSolidShape3* p_dcr, const Vec3* vec_sdof )
            {
                uint32 vec_nid[4] = { p_solid->T_VID(eid,0), p_solid->T_VID(eid,1), p_solid->T_VID(eid,2), p_solid->T_VID(eid,3) };
                m_Bs = Mat4x4( 1, 1, 1, 1,
                               vec_sdof[vec_nid[0]].x(), vec_sdof[vec_nid[1]].x(), vec_sdof[vec_nid[2]].x(), vec_sdof[vec_nid[3]].x(),
                               vec_sdof[vec_nid[0]].y(), vec_sdof[vec_nid[1]].y(), vec_sdof[vec_nid[2]].y(), vec_sdof[vec_nid[3]].y(),
                               vec_sdof[vec_nid[0]].z(), vec_sdof[vec_nid[1]].z(), vec_sdof[vec_nid[2]].z(), vec_sdof[vec_nid[3]].z() );
                m_invBs = mal::Inverse(m_Bs); //\todo MAY be singular...
                const Vec3* default_sdof( p_solid->GetVecDefaultSDOF() );
                m_invBm = mal::Inverse( Mat4x4( 1, 1, 1, 1,
                                                default_sdof[vec_nid[0]].x(), default_sdof[vec_nid[1]].x(), default_sdof[vec_nid[2]].x(), default_sdof[vec_nid[3]].x(),
                                                default_sdof[vec_nid[0]].y(), default_sdof[vec_nid[1]].y(), default_sdof[vec_nid[2]].y(), default_sdof[vec_nid[3]].y(),
                                                default_sdof[vec_nid[0]].z(), default_sdof[vec_nid[1]].z(), default_sdof[vec_nid[2]].z(), default_sdof[vec_nid[3]].z() ) ); //\todo COULD BE PRECOMPUTED in DCR.E
                m_Bs_invBm = m_Bs*m_invBm;
                m_IsValid = true;
            }
    private:
        bool m_IsValid;
        Mat4x4 m_Bs;
        Mat4x4 m_invBs;
        Mat4x4 m_invBm;
        Mat4x4 m_Bs_invBm;
    public:
        uint32 m_FirstTransformedVID;
        mutable uint32 m_Stats_NumTransform;
    };
public:
    inline dcr3_element_cache_type( const TetSolidShape3* p_solid, const DCR_TetSolidShape3* p_dcr, const Transform3& tr, const Vec3* vec_sdof )
    : m_pSolid(p_solid), m_pDCR(p_dcr), m_Tr(tr), m_vecSDOF(vec_sdof)
        {
            m_vecECD.resize( p_solid->GetNumT() );
#ifdef __ENABLE_SIMD_PRETRANSFORM //TEMP: Testing pre-reserve N/10
            m_vecTransformedV_SIMD.reserve( m_pDCR->m_NumTriangles / 10 );
#else
            m_vecTransformedV.reserve( m_pDCR->m_NumTriangles / 10 );
#endif
        }
    const ECD& operator[]( uint32 eid ) const
        {
            ECD& ecd(m_vecECD[eid]);
            if( !ecd.IsValid() ) ecd.Update( eid, m_pSolid, m_pDCR, m_vecSDOF );
            return ecd;
        }
    const Vec3* GetTransformedV( uint32 eid, const Transform3& tr ) const
        {
            ECD& ecd(m_vecECD[eid]);
            if( !ecd.IsValid() ) ecd.Update( eid, m_pSolid, m_pDCR, m_vecSDOF );
            if( ecd.m_FirstTransformedVID == 0xFFFFFFFF )
            {
                // ALLOC eid.m_NumVertices
                const DCR_TetSolidShape3::Element& ed( m_pDCR->m_vecE[eid] );
                ecd.m_FirstTransformedVID = m_vecTransformedV.size();
                m_vecTransformedV.resize( ecd.m_FirstTransformedVID + ed.m_NumVertices );
                // TRANSFORM E.V using tr and saved bary transforms
                Vec3* vec_v( &m_vecTransformedV[ecd.m_FirstTransformedVID] );
                if( true )//\todo tr != Transform3::Identity() ) THIS would accelerate element_cache2 case, where tr == Identity(), but only marginally, so there's no point
                {
                    Mat3x3 tr_F( tr.Rot() * mal::GRange<1,1,3,3>(ecd.Bs_invBm()) ); //Rot * F
                    Vec3 tr_T( tr * mal::GRange<1,3>( mal::GColumn<0>(ecd.Bs_invBm()) ) ); //Tr * T
                    // Transform ed1 m_NumVertices vertices from m_First and store them in 0-based indices inside vec_v1_2
                    for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
                        vec_v[it_vie] = tr_F * m_pDCR->m_vecV[ed.m_FirstVID + it_vie] + tr_T;
                }
                else
                {
                    Mat3x3 tr_F( mal::GRange<1,1,3,3>(ecd.Bs_invBm()) ); //Rot * F
                    Vec3 tr_T( mal::GRange<1,3>( mal::GColumn<0>(ecd.Bs_invBm()) ) ); //Tr * T
                    // Transform ed1 m_NumVertices vertices from m_First and store them in 0-based indices inside vec_v1_2
                    for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
                        vec_v[it_vie] = tr_F * m_pDCR->m_vecV[ed.m_FirstVID + it_vie] + tr_T;
                }
                ecd.m_Stats_NumTransform++;
            }
            return &m_vecTransformedV[ecd.m_FirstTransformedVID];
        };
#ifdef __ENABLE_SIMD_PRETRANSFORM
    const simd::V3f* GetTransformedV_SIMD( uint32 eid, const Transform3& tr ) const
        {
            ECD& ecd(m_vecECD[eid]);
            if( !ecd.IsValid() ) ecd.Update( eid, m_pSolid, m_pDCR, m_vecSDOF );
            if( ecd.m_FirstTransformedVID == 0xFFFFFFFF )
            {
                // ALLOC eid.m_NumVertices
                const DCR_TetSolidShape3::Element& ed( m_pDCR->m_vecE[eid] );
                ecd.m_FirstTransformedVID = m_vecTransformedV_SIMD.size();
                m_vecTransformedV_SIMD.resize( ecd.m_FirstTransformedVID + ed.m_NumVertices );
                // TRANSFORM E.V using tr and saved bary transforms
                simd::V3f* vec_v( &m_vecTransformedV_SIMD[ecd.m_FirstTransformedVID] );
                const Mat3x3 tr_F( tr.Rot() * mal::GRange<1,1,3,3>(ecd.Bs_invBm()) ); //Rot * F
                const Vec3 tr_T( tr * mal::GRange<1,3>( mal::GColumn<0>(ecd.Bs_invBm()) ) ); //Tr * T
                // Transform ed1 m_NumVertices vertices from m_First and store them in 0-based indices inside vec_v1_2
#ifdef __ENABLE_SIMD_PRETRANSFORM_DOT4
                /* SIMD transform Dot4() IMPORTANT: This changes floating point results and DOES NOT produce the SAME simulation as other branches NoSIMD, SIMD_NoDot and SIMD_Dot3
                   However, results seem "correct", and performance is slightly better in multi_A1K, must test further
                 */
                const simd::V3f row0( tr_F(0,0), tr_F(0,1), tr_F(0,2), tr_T.x() );
                const simd::V3f row1( tr_F(1,0), tr_F(1,1), tr_F(1,2), tr_T.y() );
                const simd::V3f row2( tr_F(2,0), tr_F(2,1), tr_F(2,2), tr_T.z() );
                for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
                {
                    const simd::V3f v( m_pDCR->m_vecV[ed.m_FirstVID + it_vie], 1 );
                    vec_v[it_vie] = simd::V3f( simd::Dot4(row0,v), simd::Dot4(row1,v), simd::Dot4(row2,v) );
                }
                /* This version avoids 3 _mm_cvtss_f32 using 3 sums of
                   dp instead without ever converting from _m128 to
                   float, but produces different results, not sure if it's correct.

                for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
                {
                    const simd::V3f v( m_pDCR->m_vecV[ed.m_FirstVID + it_vie], 1 );
                    vec_v[it_vie]._data = _mm_add_ps( _mm_dp_ps(row0._data,v._data,0xF1),
                                                      _mm_add_ps( _mm_dp_ps(row1._data,v._data,0xF2),
                                                                  _mm_dp_ps(row2._data,v._data,0xF4) ) );
                }
                */
#elif defined(__ENABLE_SIMD_PRETRANSFORM_DOT3)
                // SIMD transform Dot3()
                const simd::V3f row0( tr_F(0,0), tr_F(0,1), tr_F(0,2) );
                const simd::V3f row1( tr_F(1,0), tr_F(1,1), tr_F(1,2) );
                const simd::V3f row2( tr_F(2,0), tr_F(2,1), tr_F(2,2) );
                for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
                {
                    const simd::V3f v( m_pDCR->m_vecV[ed.m_FirstVID + it_vie] );
                    vec_v[it_vie] = simd::V3f( simd::Dot(row0,v), simd::Dot(row1,v), simd::Dot(row2,v) );
                    vec_v[it_vie] += simd::V3f(tr_T);
                }
#else
                // SIMD NoDot
                for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
                    vec_v[it_vie] = simd::V3f( tr_F * m_pDCR->m_vecV[ed.m_FirstVID + it_vie] + tr_T );
#endif
                ecd.m_Stats_NumTransform++;
            }
            return &m_vecTransformedV_SIMD[ecd.m_FirstTransformedVID];
        };
#endif

private:
    const TetSolidShape3* m_pSolid;
    const DCR_TetSolidShape3* m_pDCR;
    const Transform3& m_Tr;
    const Vec3* m_vecSDOF;
    mutable std::vector<ECD> m_vecECD; //\todo SHOULD be scope-allocated
    mutable std::vector<Vec3> m_vecTransformedV; //\todo SHOULD be scope-allocated or .reserve() a reasonable fraction (eg num_vertices/10) at startup
#ifdef __ENABLE_SIMD_PRETRANSFORM
    mutable std::vector<simd::V3f> m_vecTransformedV_SIMD; //\todo SHOULD be scope-allocated or .reserve() a reasonable fraction (eg num_vertices/10) at startup
#endif
};


//! Compact half_edge_id representation as <triangle_id:22b, edge_in_triangle_index:2b>
class dcr3_half_edge_id_type
{
public:
    inline dcr3_half_edge_id_type() : m_Bits(0xFFFFFFFF), m_PID(0xFFFFFFFF) {}
    inline dcr3_half_edge_id_type( uint32 pid, uint32 tid, uint32 eit )
    : m_Bits( (tid << 2) + eit ), m_PID(pid) { GEO_ASSERT( tid < 0x3FFFFFFF && eit < 4 ); }
    inline bool operator==( const dcr3_half_edge_id_type& other ) const { return m_Bits == other.m_Bits; }
    inline uint32 PID() const { return m_PID; }
    inline uint32 TID() const { return m_Bits >> 2; }
    inline uint32 EIT() const { return m_Bits & 0x3; }
    inline operator uint32() const { return m_Bits; } //\todo Ugly way to use hash<uint32> in unordered_map/set
private:
    uint32 m_Bits;
    uint32 m_PID; //useful to avoid re-searching it
};

/* Flood IB (on a given side object 0/1)

   IMPORTANT: This comment block is partially deprecated by the actual
   implementation... TRUST THE CODE.

   CTP:
   - VF: Flood only on V1-side through V1-edges, no flood on F2-side
   - FV: Symmetric
   - EE: Flood on both sides through E1 and E2 edges

   Alg:
    0) Get an unconsumed IC and create its IB
    1) FloodIC()
       Iterate over all IC.CTP and:
       1.1) For each piercing edge E, if the inwards vertex V is INTERIOR, add CCW interior triangles {T} adjacent to E->V
       1.2) For each pierced triangle T, gather all cutting segments S (from both VF and EE CTP)
       1.3) For each pierced triangle T,
            - process all (normal, are, etc...) using cutting segments S
              and mark as flooded reagarding IC
            - ONLY process the area of T internal to IC
    2) Flood through interior triangles T starting from E->V->{T} seeds
       2.1) Process T and mark T flooded
       2.2) For all T-adjacent triangles
            - If not flooded
              - If pierced
                - If CTP are all in the same IC, T was already handled in 2.1, ignore
                - If T has CTP in a different IC2, IC and IC2 become "connected"
                  - Perform 1) and 2) on IC2, as IC and IC2 bound the
                    same IB and they should be flooded together, or
                    DEFER TO step 4) where they will be merged after
                    being flooded individually
              - Else
                - Add to flood stack
    3) Global processing of a fully-internal patch P requires:
       - Detecting T1->T2 flood-into-different-patch
         - TID are consecutive in a DCR.P
         - VID are consecutive in a DCR.P
       - Testing IsCrossing(P)
         - vecIsCrossingP
       - NOT PROCESSING any T individually if it's in a fully-internal patch
         - Detecting T is in a different patch
           - Ignore them in step 2.2
       - Flood to neighbour patches P1->P2
         - At the patch level if P2 is fully-internal
           - Requires P1->P2 topology
         - Flood through P1->P2 shared edges otherwise
           - Requires iteration over P1->P2 shared edges
    4) Merge IC/IB with pre-existing connected IC2/IB2
       - AFTER flooding a given IC and generating its IB, connections
         with other IC2/IB2 may have been found.
         - Merge them (should be FAST)
       - IMPORTANT: This requires setFloodedT and setFloodedP to be
         GLOBAL, as an IC2 must not flood into an internal T/P already
         associated to a previous IC1

    IMPORTANT:
    - Vertices are SPLIT at DCR.E boundaries... this means that
      symmetric edges accross neighbour DCR.P DO NOT refer the SAME
      vertex pair.
    - Triangles and Edges may be pierced by SEVERAL CTP that belong to
      DIFFERENT IC, thus, any "check if pierced" must be done by
      testing the SSH for adjacent matching CTP, not only the CTP on a
      given IC.
    - Edge-flooding needs to account for the lambda param along the
      edge, due to potential multiple in/out intervals along it.
    - Piercing edges E may be pierced multiple times, though, any
      in/out analysis must consider multiple in/out intervals along E
    - In FloodIC() Piercing edges E will only be considered for a
      given CTP if INWARDS. Otherwise, the symmetric edge will exist
      and will be processed instead.
    - Pierced triangles T may have SEVERAL DISJOINT areas interior to
      the either the same IC or a DIFFERENT one.

    Implementation Phases:
    +A) Basic flood on DCR.T, no DCR.P flood, no IB.P, no aggregation, no IC connection
       => Draw flooded DCR.T
    B) Process crossing DCR.T
       => Draw crossing DCR.T
    D) DCR.T flood with IB.P creation and aggregate data
       => Draw aggregate data
    E) DCR.P flood with O(1) aggregation
       => Enable/Disable at runtime to check identical results but faster
    C) Handle IC connection
       => Draw IC connection on an IB

   SIMPLIFICATION:
   - If we accept flooding all IC at the same time, we can AVOID
     having to flood from fully-internal DCR.P to its boundary
     DCR.P.T.
     - Start flood from any IC1, flood inwards through piercing DCR.T
       and, if we arrive at an internal DCR.P, start flooding DCR.P
       while internal, NEVER transitioning back to DCR.T from a DCR.P
       (but continue with already open DCR.T anyway)
     - This flood will finish when all IC1 piercing DCR.T and all
       accessible DCR.P have been visited. The IC1 will have an
       associated IB1.
       - If the IB1 is simple, we've finished, all DCR.T have been
         flooded either individually or included in internal DCR.P
       - OTHERWISE, we have NOT ARRIVED at an IB1-connected IC2 != IC1.
       - Now, we can start flooding any other unconsumed IC_j in the
         same way, stopping when we find any fully internal DCR.P that
         was ALREADY included in a previously floodedd IC_i. IC_i and
         IC_j become connected, and their "finished" IB_i and
         "ongoing" IB_j can be merged into the SAME IB...
         - EP do NOT NEED to be merged, only concatenated.
   => This AVOIDS having to store the
      DCR.P.NeighbourPatches[npip].vecBoundaryTID[] array, saving
      memory and reducing complexity.

   IMPORTANT: This comment block is partially deprecated by the actual
   implementation... TRUST THE CODE.

   \todo Subst std::vector stack with std::array
*/
bool FloodIB( std::vector<dcr3_intersection_boundary_type>& vecIB,
              const std::vector<crossing_triangle_pair_type>& vecCTP,
              std::vector<dcr3_intersection_curve_type>& vecIC,
              const TetSolidShape3* p_solid, const DCR_TetSolidShape3* p_dcr, const Transform3& tr, const Vec3* vec_sdof,
              int side,
              bool b_inwards, //otherwise outwards
              const dcr3_element_cache_type& element_cache,
              //Results
              ContactData3& cd, //only required to store debug/viz
              const Context* p_context = g_pDefaultContext )
{
    bool bLog(p_context->m_DCR2DCR_Log_Enabled);
    bool bViz(p_context->m_DCR2DCR_Viz_Enabled);
    int oside( 1-side );

    // Perform FloodT/P if enabled
    if( p_context->m_DCR2DCR_IB_Method != Context::eDCR2DCR_IB_None )
    {
        GEO_NP_BEGIN_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_FloodIB_IT_ms );

        // Flag all crossing DCR.P
        std::vector<bool> vecIsCrossingP(p_dcr->m_NumPatches,false);
        for( const auto& ctp : vecCTP )
            vecIsCrossingP[ ctp.m_PID[side] ] = true;

        // Gather set of crossing triangles in all IC
        std::unordered_set< uint32 > setCrossingTID;
        for( const auto& ic : vecIC )
            for( auto ctpid : ic.m_vecCTPID )
                setCrossingTID.insert( vecCTP[ctpid].m_TID[side] );

        /* Flood each IC from its inwards edges, navigating only to fully
           internal T/P.

           \note An internal/non-crossing T or P CANNOT be flooded from
           different IC if they are not connected by a single IB, thus,
           internal T and P are unambiguous witnesses for IC/IB
           merging. Crossing T are complementary IC/IB merge witnesses
           that will be handled after the flood.
        */
        std::unordered_map<uint32,dcr3_intersection_boundary_type::Patch> mapFloodedInternalP; //flooded internal Patches
        std::unordered_map<uint32,uint32> mapFloodedInternalTID; //flooded internal Triangles with their parent IBID
        uint32 num_flood_iter(0);
        uint32 num_flood_triangles(0);
        uint32 num_flood_patches(0);
        for( unsigned int it_ic=0; it_ic<vecIC.size(); it_ic++ )
        {
            uint32 icid(it_ic);
            uint32 ibid(it_ic);
            const dcr3_intersection_curve_type& ic( vecIC[icid] );
            //---- 0) Init IB
            /* Create an empty IB for each IC, that may be MF-Set merged
               afterwards. Each IC will be flooded individually, contributing
               to its exclusive IB. The relationship IC<->IB will be 1:1
               during per-IC floods.
               After all IC/IB are flooded, the IB may be compacted into their
               root-IB.
            */
            vecIB.push_back( dcr3_intersection_boundary_type(icid) );
            dcr3_intersection_boundary_type& ib( vecIB.back() );
            ib.m_ParentId = it_ic;
            if( bViz )
            {
                cd.m_VD.m_vecDCR2DCR_IB_vecInternalT_Side_x_IB[side].push_back( std::vector< util::triad<Vec3,Vec3,Vec3> >() ); //add [ibid] entry
            }

            //---- 1) Process the whole IC
            /* We store all DCR half-edges that appear in the IC and are
               inwards on this side in a map together with their interior
               interval, which will be clipped successively if other CTP
               reference it. After all CTP have been processed, all
               inwards HE that have a non-empty interval that includes the
               endpoint (which is inside DCR[side]) will be used as FloodT
               seeds.
            */
            std::unordered_map< dcr3_half_edge_id_type, Interval, std::hash<uint32> > mapHE; //<HEID,Interval>
            for( auto ctpid : ic.m_vecCTPID )
            {
                const crossing_triangle_pair_type& ctp( vecCTP[ctpid] );
                //1.1) For each piercing edge E, if the inwards vertex V is INTERIOR, add CCW interior triangles {T} adjacent to E->V
                /* TODO: Sometimes the the flood covers BOTH SIDES from an
                   IC, it may be caused by an incorrectly classified IN/OUT
                   direction or from a flood spill due to unclosed IC, make
                   this robust or the approach is USELESS.
                */
                const intersection_triangle3_triangle3_result_type& ittr( ctp.m_ITTR );
                if( (ittr.IsVF() && side == 0) || (ctp.m_ITTR.IsFV() && side == 1) )
                {
                    Interval i0, i1;
                    Vec3 other_normal_0( mal::Cross( ctp.m_vecPos_0[oside][1]-ctp.m_vecPos_0[oside][0], ctp.m_vecPos_0[oside][2]-ctp.m_vecPos_0[oside][0] ) );
                    if( mal::Dot( ctp.m_vecPos_0[side][(ittr.m_EdgeInT0+1)%3] - ctp.m_vecPos_0[side][ittr.m_EdgeInT0], other_normal_0 ) < 0 ) //E0 inwards
                    {
                        i0 = Interval( ittr.m_Lambda0, 1 );
                        i1 = Interval( 0, ittr.m_Lambda1 );
                    }
                    else if( mal::Dot( ctp.m_vecPos_0[side][(ittr.m_EdgeInT1+1)%3] - ctp.m_vecPos_0[side][ittr.m_EdgeInT1], other_normal_0 ) < 0 ) //E0 outwards ==> E1 must be inwards \todo unless numeric problems strike again
                    {
                        i1 = Interval( ittr.m_Lambda1, 1 );
                        i0 = Interval( 0, ittr.m_Lambda0 );
                    }
                    else
                    {
                        //\todo else, no inwards edge, DON't add the HE!?
                        if( bLog ) GEO_LOG("!!!!!!!!!!!!!!!!!!! Coincident edge !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    }
                    // add or clip he0,he1 according to inwards/outwards
                    dcr3_half_edge_id_type heid0( ctp.m_PID[side], ctp.m_TID[side], ittr.m_EdgeInT0 );
                    dcr3_half_edge_id_type heid1( ctp.m_PID[side], ctp.m_TID[side], ittr.m_EdgeInT1 );
                    auto it_he0( mapHE.find( heid0 ) );
                    auto it_he1( mapHE.find( heid1 ) );
                    if( it_he0 == mapHE.end() ) mapHE[heid0] = i0;
                    else it_he0->second.Intersect( i0 );
                    if( it_he1 == mapHE.end() ) mapHE[heid1] = i1;
                    else it_he1->second.Intersect( i1 );
                }
                else if( ctp.m_ITTR.IsEE() )
                {
                    Vec3 other_normal_0( mal::Cross( ctp.m_vecPos_0[oside][1]-ctp.m_vecPos_0[oside][0], ctp.m_vecPos_0[oside][2]-ctp.m_vecPos_0[oside][0] ) );
                    if( side == 0 )
                    {
                        // add or clip inwards/outwards he0
                        Interval i0;
                        if( mal::Dot( ctp.m_vecPos_0[side][(ittr.m_EdgeInT0+1)%3] - ctp.m_vecPos_0[side][ittr.m_EdgeInT0], other_normal_0 ) < 0 )
                            i0 = Interval( ittr.m_Lambda0, 1 );
                        else
                            i0 = Interval( 0, ittr.m_Lambda0 );
                        dcr3_half_edge_id_type heid0( ctp.m_PID[side], ctp.m_TID[side], ittr.m_EdgeInT0 );
                        auto it_he0( mapHE.find( heid0 ) );
                        if( it_he0 == mapHE.end() ) mapHE[heid0] = i0;
                        else it_he0->second.Intersect( i0 );
                    }
                    else //side == 1
                    {
                        // add or clip inwards/outwards he1
                        Interval i1;
                        if( mal::Dot( ctp.m_vecPos_0[side][(ittr.m_EdgeInT1+1)%3] - ctp.m_vecPos_0[side][ittr.m_EdgeInT1], other_normal_0 ) < 0 )
                            i1 = Interval( ittr.m_Lambda1, 1 );
                        else
                            i1 = Interval( 0, ittr.m_Lambda1 );
                        // add or clip inwards/outwards he1
                        dcr3_half_edge_id_type heid1( ctp.m_PID[side], ctp.m_TID[side], ittr.m_EdgeInT1 );
                        auto it_he1( mapHE.find( heid1 ) );
                        if( it_he1 == mapHE.end() ) mapHE[heid1] = i1;
                        else it_he1->second.Intersect( i1 );
                    }
                }
            }

            //---- 2) Gather flood seeds
            // Use INWARDS piercing edges to find adjacent INTERNAL triangles and add them as flood seeds
            struct flood_triangle_type { uint32 m_PID; uint32 m_TID; inline flood_triangle_type( uint32 pid, uint32 tid ) : m_PID(pid), m_TID(tid) {} };
            std::vector< flood_triangle_type > stackT; //triangle flood stack    \todo Subst std::vector stack with std::array
            for( const auto& it_he : mapHE )
            {
                //find other CTP on the same edge, check if any has lambda that avoids E->V reaching V
                if( !it_he.second.IsEmpty() && it_he.second.Max() == 1 ) //ONLY non-empty he intervals that end at 1 have the end vertex INTERIOR
                {
                    dcr3_half_edge_id_type heid( it_he.first );
                    // Add all T adjacent to HE end vertex but NOT he.left/he.right or any crossing T
                    uint32 tid( heid.TID() );
                    uint32 pid( heid.PID() );
                    uint32 eit( heid.EIT() );
                    // from tid, advance 2 triangles CCW (first shares the piercing edge and thus is crossing, second may be internal)
                    uint32 ntid( p_dcr->m_vecT[tid].GetNTID(eit) );
                    uint32 neint( p_dcr->Find_NEINT(tid,ntid) );
                    uint32 npid( p_dcr->Find_PID_From_TID_Hint(ntid,pid) );
                    uint32 nntid( p_dcr->m_vecT[ntid].GetNTID( (neint+2) % 3 ) );
                    uint32 nneint( p_dcr->Find_NEINT(ntid,nntid) );
                    uint32 nnpid( p_dcr->Find_PID_From_TID_Hint(nntid,npid) );
                    /*BUG: with armadillo10K vs box the following loop
                     * never ends some times, with tid=12?49 pid=505, seems
                     * to step into an inconsistent loop of faces around a
                     * vertex that DOES NEVER get back to the original
                     * tid, BUT yet cycles (nntid is periodic, but NEVER
                     * goes back to tid)... this is probably an error in
                     * armadillo topology that WE SHOULD DETECT at
                     * build-time and fix if possible...
                     */
                    while( nntid != tid )
                    {
                        /*IMPORTANT: Ideally, we would ONLY consider
                          Merging HERE, NOT DURING Flood.

                          \todo WHILE we COULD handle Patch-level stuff
                          here, it's a LOT SIMPLER if we just stackT and
                          the flood converts suitable stackT entries into
                          stackP entries AT A SINGLE POINT

                          IMPORTANT: If I'm lucky, THIS IS THE
                          ONLY place where an IC/IB may become
                          merged with a previous one, because
                          there CANNOT be fully internal T/P that
                          are not accessible from the crossing T
                          of ANY connected IC.

                          IMPORTANT!!!! If THIS is the only merge-point,
                          and we always merge to a previous IB, we could
                          AVOID MERGE alltogether by NOT CREATING the IB
                          until we've confirmed there's no merge and, if
                          merge, we just retrieve that very same ib/IBID
                          and expand it with this IC-flood.

                          => The caveat are potential merges through
                          CROSSING triangles, at the end... if we manage
                          to solve them without whole-IB merging, we can
                          avoid merge overhead everywhere.
                        */
                        if( !vecIsCrossingP[nnpid] )
                        {
                            auto it_nnp( mapFloodedInternalP.find(nnpid) );
                            if( it_nnp == mapFloodedInternalP.end() ) // Unflooded
                                stackT.push_back( flood_triangle_type(nnpid,nntid) );
                            else if( ib.m_ParentId != it_nnp->second.m_ParentIBID ) //flooded and different parent IB
                            {
                                // Merge Root(ib.m_ParentId) and Root(it_nnp->second.m_ParentIBID)
                                uint32 ib_root_id( IB_ComputeRootId( vecIB, ib.m_ParentId ) );
                                uint32 nnp_root_id( IB_ComputeRootId( vecIB, it_nnp->second.m_ParentIBID ) );
                                if( ib_root_id != nnp_root_id )
                                {
                                    vecIB[ib_root_id].m_ParentId = nnp_root_id;
                                    ib.m_ParentId = nnp_root_id; //path-compresion
                                    if( bLog ) { GEO_LOG( "MF-Set connecting IB[%u] with root %u to root IB[%u] through NNP[%u] during FloodIC", ibid, ib_root_id, nnp_root_id, nnpid ); }
                                }
                            }
                        }
                        else if( setCrossingTID.find(nntid) == setCrossingTID.end() )
                        {
                            auto it_nnt( mapFloodedInternalTID.find(nntid) );
                            if( it_nnt == mapFloodedInternalTID.end() ) // Unflooded
                                stackT.push_back( flood_triangle_type(nnpid,nntid) );
                            else if( ib.m_ParentId != it_nnt->second ) //flooded and different parent IB
                            {
                                // Merge Root(ib.m_ParentId) and Root(it_nnt->second)
                                uint32 ib_root_id( IB_ComputeRootId( vecIB, ib.m_ParentId ) );
                                uint32 nnt_root_id( IB_ComputeRootId( vecIB, it_nnt->second ) );
                                if( ib_root_id != nnt_root_id )
                                {
                                    vecIB[ib_root_id].m_ParentId = nnt_root_id;
                                    ib.m_ParentId = nnt_root_id; //path-compresion
                                    if( bLog ) { GEO_LOG( "MF-Set connecting IB[%u] with root %u to root IB[%u] through NNT[%u] during FloodIC", ibid, ib_root_id, nnt_root_id, nntid ); }
                                }
                            }
                        }

                        // Advance
                        ntid = nntid;
                        neint = nneint;
                        npid = nnpid;
                        nntid = p_dcr->m_vecT[ntid].GetNTID( (neint+2) % 3 );
                        nneint = p_dcr->Find_NEINT(ntid,nntid);
                        nnpid = p_dcr->Find_PID_From_TID_Hint(nntid,npid);
                    }

                    //DEBUG seed HE
                    if( bViz )
                    {
                        uint32 eid( p_dcr->m_vecP[heid.PID()].m_EID );
                        Vec3 edge_p0( tr * mal::GRange<1,3>( element_cache[eid].Bs_invBm() * mal::Concat(1,p_dcr->m_vecV[p_dcr->m_vecT[heid.TID()].GetVID(heid.EIT()) ] ) ) );
                        Vec3 edge_p1( tr * mal::GRange<1,3>( element_cache[eid].Bs_invBm() * mal::Concat(1,p_dcr->m_vecV[p_dcr->m_vecT[heid.TID()].GetVID((heid.EIT()+1)%3) ] ) ) );
                        cd.m_VD.m_vecDCR2DCR_IB_vecSeedHE_Side[side].push_back( std::make_pair( (1-it_he.second.Min()) * edge_p0 + it_he.second.Min() * edge_p1,
                                                                                                (1-it_he.second.Max()) * edge_p0 + it_he.second.Max() * edge_p1 ) );
                    }
                }
            }

            //---- 3) Flood internal Triangles
            std::vector<uint32> stackP; //patch flood stack, initially empty as stackT only contains crossing triangles    \todo Subst std::vector stack with std::array
            while( !stackT.empty() || !stackP.empty() )
            {
                num_flood_iter++;
                // 3.1) FloodT
                if( !stackT.empty() )
                {
                    // pop
                    flood_triangle_type ft( stackT.back() );
                    stackT.pop_back();
                    uint32 tid( ft.m_TID );
                    uint32 pid( ft.m_PID );
                    uint32 eid( p_dcr->m_vecP[pid].m_EID );
                    /*\todo Phase E) If we want to simplify neighbour
                      evaluation here, we should detect that we're trying
                      to enter a non-crossing P and, if not-yet Flooded,
                      zoom out and push it onto the stackP.

                      It MAY NOT BE OPTIMAL, but the code is a LOT SIMPLER
                      if THIS IS THE ONLY PLACE where we zoom out from
                      stackT to stackP. This way we DO NOT NEED TO
                      consider patch-level flooding at IC seeding or at
                      neighbour expansion.
                    */
                    // 3.1.1) Analyze P and T status (ALWAYS consider P-level status first)
                    if( p_context->m_DCR2DCR_IB_Method == Context::eDCR2DCR_IB_FloodP
                        && !vecIsCrossingP[pid] )
                    {
                        auto it_np( mapFloodedInternalP.find(pid) );
                        if( it_np == mapFloodedInternalP.end() )
                            stackP.push_back( pid );
                        else if( ib.m_ParentId != it_np->second.m_ParentIBID ) //flooded and different parent IB
                        {
                            // Merge Root(ib.m_ParentId) and Root(it_np->second.m_ParentIBID)
                            uint32 ib_root_id( IB_ComputeRootId( vecIB, ib.m_ParentId ) );
                            uint32 np_root_id( IB_ComputeRootId( vecIB, it_np->second.m_ParentIBID ) );
                            if( ib_root_id != np_root_id )
                            {
                                vecIB[ib_root_id].m_ParentId = np_root_id;
                                ib.m_ParentId = np_root_id; //path-compresion
                                if( bLog ) { GEO_LOG( "MF-Set connecting IB[%u] with root %u to root IB[%u] through NP[%u] during FloodT", ibid, ib_root_id, np_root_id, pid ); }
                            }
                        }
                    }
                    else if( mapFloodedInternalTID.find( tid ) == mapFloodedInternalTID.end() )
                    {
                        mapFloodedInternalTID.insert( std::make_pair(tid,ibid) );
                        num_flood_triangles++;

                        const DCR_TetSolidShape3::Triangle& td( p_dcr->m_vecT[tid] );
                        //\todo CONSIDER using element_cache transformed V, not sure...
                        Vec3 tri0_0( tr * mal::GRange<1,3>( element_cache[eid].Bs_invBm() * mal::Concat(1,p_dcr->m_vecV[td.GetVID(0)] ) ) );
                        Vec3 tri1_0( tr * mal::GRange<1,3>( element_cache[eid].Bs_invBm() * mal::Concat(1,p_dcr->m_vecV[td.GetVID(1)] ) ) );
                        Vec3 tri2_0( tr * mal::GRange<1,3>( element_cache[eid].Bs_invBm() * mal::Concat(1,p_dcr->m_vecV[td.GetVID(2)] ) ) );
                        Vec3 v_area_0( mal::Cross(tri1_0-tri0_0,tri2_0-tri0_0) );
                        // Aggregate whole-T stuff to parent P in parent IB
                        {
                            // If parent patch is crossing, add in current IB, otherwise, add in its ROOT IB
                            dcr3_intersection_boundary_type::Patch* pPatch(0);
                            if( vecIsCrossingP[pid] )
                            {
                                auto it_patch( ib.m_mapCrossingP.find(pid) );
                                if( it_patch == ib.m_mapCrossingP.end() )
                                    pPatch = &ib.m_mapCrossingP.insert( std::make_pair( pid, dcr3_intersection_boundary_type::Patch(pid,ibid) ) ).first->second;
                                else
                                    pPatch = &it_patch->second;
                            }
                            else //\todo if eDCR2DCR_IB_FloodP, this branch is never entered, handled though stackP
                            {
                                auto it_patch( mapFloodedInternalP.find(pid) );
                                if( it_patch == mapFloodedInternalP.end() )
                                    pPatch = &mapFloodedInternalP.insert( std::make_pair( pid, dcr3_intersection_boundary_type::Patch(pid,ibid) ) ).first->second;
                                else
                                    pPatch = &it_patch->second;
                            }
                            pPatch->m_NumTriangles++;
                            // pPatch->m_AvgPos_0 += Real(0.333333)*(tri0_0+tri1_0+tri2_0); //\todo this may be wrong, we need something that can be aggregated for crossing T and transformed for internal P
                            Real area( Real(0.5)*mal::Norm(v_area_0) );
                            pPatch->m_AvgPos_0 += area*Real(0.333333)*(tri0_0+tri1_0+tri2_0); //\todo this may be wrong, we need something that can be aggregated for crossing T and transformed for internal P
                            pPatch->m_VectorArea_0 += v_area_0;
                            pPatch->m_AreaSq += mal::Sq(area);
                            pPatch->m_Area += area;
                        }
                        // Flood to non-crossing, non-flooded neighbours
                        for( int it_ntit=0; it_ntit<3; it_ntit++ )
                        {
                            uint32 ntid( p_dcr->m_vecT[tid].GetNTID(it_ntit) );
                            GEO_ASSERT( ntid != DCR_TetSolidShape3::Triangle::cInvalidTID );
                            /* Phase E) Here we ALLOW flooding into a
                               DIFFERENT npid that MAY be non-crossing and
                               already flooded. We'll detect both cases at
                               3.1.1), NOT HERE, to avoid code duplication.

                               IMPORTANT: If FloodP, here ntid MAY BE in a
                               flooded internal npid that is in a PREVIOUS
                               IB_0, and therefore IB and IB_0 NEED to be
                               merge. This situation may arise when the
                               currently flooded IB is only adjacent to
                               IB_0 through a closed band of internal P
                               that effectively "isolates" the current IC
                               from the previous IB at the T-level. This
                               is handled at 3.1.1)
                            */
                            if( setCrossingTID.find(ntid) == setCrossingTID.end()
                                &&
                                mapFloodedInternalTID.find(ntid) == mapFloodedInternalTID.end() )
                            {
                                uint32 npid( p_dcr->Find_PID_From_TID_Hint(ntid,pid) );
                                stackT.push_back( flood_triangle_type(npid,ntid) );
                            }
                        }
                        //DEBUG internal triangles
                        if( bViz )
                        {
                            // cd.m_VD.m_vecDCR2DCR_IB_vecInternalT_Side[side].push_back( util::make_triad(tri0,tri1,tri2) ); //superseded by m_vecDCR2DCR_IB_vecInternalT_Side_x_IB
                            cd.m_VD.m_vecDCR2DCR_IB_vecInternalT_Side_x_IB[side][ibid].push_back( util::make_triad(tri0_0,tri1_0,tri2_0) );
                        }
                    }
                }

                // 3.2) FloodP
                if( !stackP.empty() )
                {
                    // pop
                    uint32 pid( stackP.back() );
                    stackP.pop_back();
                    uint32 eid( p_dcr->m_vecP[pid].m_EID );
                    // 3.2.1) Analyze P status
                    if( mapFloodedInternalP.find( pid ) == mapFloodedInternalP.end() )
                    {
                        dcr3_intersection_boundary_type::Patch& patch( mapFloodedInternalP.insert( std::make_pair( pid, dcr3_intersection_boundary_type::Patch(pid,ibid) ) ).first->second );
                        num_flood_patches++;
                        //\todo PROCESS P as a whole by affine transform
                        GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_FloodIB_IT_ms );
                        const DCR_TetSolidShape3::Patch& pd( p_dcr->m_vecP[pid] );
                        patch.m_NumTriangles += pd.m_NumTriangles;
                        for( unsigned int it_tip=0; it_tip<pd.m_NumTriangles; it_tip++ )
                        {
                            uint32 tid( pd.m_FirstTID + it_tip );
                            const DCR_TetSolidShape3::Triangle& td( p_dcr->m_vecT[tid] );
                            Vec3 tri0_0( tr * mal::GRange<1,3>( element_cache[eid].Bs_invBm() * mal::Concat(1,p_dcr->m_vecV[td.GetVID(0)] ) ) );
                            Vec3 tri1_0( tr * mal::GRange<1,3>( element_cache[eid].Bs_invBm() * mal::Concat(1,p_dcr->m_vecV[td.GetVID(1)] ) ) );
                            Vec3 tri2_0( tr * mal::GRange<1,3>( element_cache[eid].Bs_invBm() * mal::Concat(1,p_dcr->m_vecV[td.GetVID(2)] ) ) );
                            Vec3 v_area_0( mal::Cross(tri1_0-tri0_0,tri2_0-tri0_0) );
                            // patch.m_AvgPos_0 += Real(0.333333)*(tri0_0+tri1_0+tri2_0); //\todo this may be wrong, we need something that can be aggregated for crossing T and transformed for internal P
                            Real area( Real(0.5)*mal::Norm(v_area_0) );
                            patch.m_AvgPos_0 += area*Real(0.333333)*(tri0_0+tri1_0+tri2_0);
                            patch.m_VectorArea_0 += v_area_0;
                            patch.m_AreaSq += mal::Sq(area);
                            patch.m_Area += area;
                            //DEBUG internal triangles
                            if( bViz )
                            {
                                // cd.m_VD.m_vecDCR2DCR_IB_vecInternalT_Side[side].push_back( util::make_triad(tri0,tri1,tri2) ); //superseded by m_vecDCR2DCR_IB_vecInternalT_Side_x_IB
                                cd.m_VD.m_vecDCR2DCR_IB_vecInternalT_Side_x_IB[side][ibid].push_back( util::make_triad(tri0_0,tri1_0,tri2_0) );
                            }
                        }
                        GEO_NP_BEGIN_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_FloodIB_IT_ms );
                        // Flood to non-crossing, non-flooded neighbours
                        for( unsigned int it_npip=0; it_npip<pd.m_NumNeighbours; it_npip++ )
                        {
                            uint32 npid( p_dcr->m_vecNPID[pd.m_FirstIndexNPID + it_npip] );
                            if( !vecIsCrossingP[npid]
                                &&
                                mapFloodedInternalP.find( npid ) == mapFloodedInternalP.end() )
                            {
                                stackP.push_back( npid );
                            }
                            /*else NPID is either crossing, which will be
                              flooded at the T-level from an inwards edge
                              at some IC, or already flooded, which needs
                              no further processing.
                            */
                        }
                    }
                }
            }
        }
        if( bLog )
        {
            GEO_LOG( "Flood[%d] #tri = %u, #patch = %u, #iter/#tri = %f, #iter/#patch = %f",
                     side, num_flood_triangles, num_flood_patches,
                     (num_flood_triangles>0) ? float(num_flood_iter)/num_flood_triangles : 0,
                     (num_flood_patches>0) ? float(num_flood_iter)/num_flood_patches : 0 );
        }

        // Aggregate InternalP to their root IB \todo MUST do it BEFORE merging IB, because afterwards the original IBID become invalid
        for( auto patch : mapFloodedInternalP )
        {
            uint32 root_ibid( IB_ComputeRootId(vecIB,patch.second.m_ParentIBID) );
            vecIB[root_ibid].m_vecP.push_back( patch.second );
        }
        GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_FloodIB_IT_ms );
        GEO_NP_STAT_ADD( mig2015.m_Num_IB_FT, num_flood_triangles );
        GEO_NP_STAT_ADD( mig2015.m_Num_IB_FP, num_flood_patches );
    }
    else // If no FloodT/P, create an empty IB for each IC, required by CT treatment
    {
        for( unsigned int it_ic=0; it_ic<vecIC.size(); it_ic++ )
        {
            uint32 icid(it_ic);
            uint32 ibid(it_ic);
            //---- 0) Init IB
            /* Create an empty IB for each IC, that may be MF-Set merged
               afterwards. Each IC will be flooded individually, contributing
               to its exclusive IB. The relationship IC<->IB will be 1:1
               during per-IC floods.
               After all IC/IB are flooded, the IB may be compacted into their
               root-IB.
            */
            vecIB.push_back( dcr3_intersection_boundary_type(icid) );
            dcr3_intersection_boundary_type& ib( vecIB.back() );
            ib.m_ParentId = it_ic;
        }
    }

    // Process Clipped Triangles (CT)
    if( p_context->m_DCR2DCR_IB_EnableCT )
    {
        GEO_NP_BEGIN_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_FloodIB_CT_ms );

        /*CT treatment:
          - For each CT
          - Gather all CTP
          - Build all CCW closed polygons formed by adjacent CTP-segment
          edges and original triangle edges
          - For each closed polygon
          - Integrate stuff and accumulate into its parent crossingP
          - If formed by CTP-segment edges from DIFFERENT IC, merge
          these IC into a single IB

          - A given T may have CTP from different IC that MAY or MAY NOT bound the same IB.
          - A CT may have several disjoint regions interior to several different IB.
          - Each region of a crossing T must only contribute to its related IB, that may be bound by different IC.

          IMPORTANT: A shared CT MAY BE THE ONLY way in which
          Different IC that are directly connected (eg: imagine a
          tesellated cylinder cut perpendicularly by a box so flat
          that no cyl-triangle is completely inside. There will be a
          "strip" of triangles around the cyl that are clipped by the
          box top and down faces, with no connection through fully
          internal triangles. the ONLY way to detect that upper and
          lower IC bound a single IB on the cyl is realizing that
          there are CT in the cyl that have closed polygons bound by
          CTP from upper and lower IC)

          IMPORTANT: Ideally, given a consistent orientation (eg: CCW==inside)
          => It should be possible to use the add/substract
          integration trick of the VectorArea/Area on the raw list of
          oriented segments, WITHOUT building the actual
          open/closed/nested CTP-polygonals at all.

          - CAN we compute the "barycenter" using a similar trick?...
          => Maybe signed-area-weighted barycenter = sum of T.area*T.barycenter / TotalArea...
          => More like a CoM... can we affine-transform this?
          - If not, maybe use squared area...
          => ACTUALLY, for clipped-T we can use normal area
          safely, and square it afterwards if req. ONLY
          !crossingP may need the area^2 trick.
        */
        // For each pierced triangle T, gather all cutting segments S (from both VF and EE CTP)
        std::unordered_map< uint32, std::vector<uint32> > mapCT; //<TID,vecClippingCTPID>
        for( unsigned int it_ctp=0; it_ctp<vecCTP.size(); it_ctp++ )
            mapCT[ vecCTP[it_ctp].m_TID[side] ].push_back(it_ctp);

        /* There are 2 global cases:
           a) There are at least 2 CTP that clip T edges and define
              in/out regions limited by open polygonals that start/end at
              T edges.
              - T edges can contribute to several disjoint edge
                intervals that are internal to the DCR and "close" the
                open polygonals to defining internal areas.
                => In order to find this internal intervals, we sort all
                   CTP clips along the triangle perimeter s=[0..3), find
                   the smallest clip s, determine if it's inwards or
                   outwards and perform a whole cycle over all clips, in
                   CCW order, which are guaranteed to be alternating
                   inwards/outwards, so we create internal sub-segments
                   for (inwards,outwards) sub-intervals, considering any T
                   vertices inbetween.
              - If consecutive CTP that define an (inwards,outwards)
                sub-interval belong to different IC, these must be
                merged.

           b) There is no CTP that clips a any T edge, so T must
              completely enclose one or several IC:
              b.1) A negative hole, if all 3 T vertices and edges are inside DCR
              b.2) A positive hole, if all 3 T vertices and edges are outside DCR
              B) Any succession of nested, non-crossing, positive and negative holes.
              => We DO NOT WANT to test if T.V or T.E are
                 inside/outside DCR (requires 1 global raycast),
                 instead, we can discover it indirectly thanks to the
                 orientation convention: CCW on side[0]
                 - Compute VectorArea for all IC and, if negative,
                   the outermost "hole" must be negative, and
                   therefore the triangle was "positive" but we did
                   not take its area contribution into account, so we
                   add it (positive) to the accumulated vector area
                   to obtain a positive (along surface normal) total
                   as expected.
           a) + b) Can appear simultaneously on the same T (nested
                   holes and T-edge-clipping CTP), but cover disjoint
                   regions and therefore can be handled independently.
        */
        unsigned int num_ct(0);
        unsigned int num_ct_boundary_segments(0);
        for( auto& ct : mapCT )
        {
            uint32 tid( ct.first );
            uint32 pid( vecCTP[ct.second.front()].m_PID[side] );
            uint32 eid( p_dcr->m_vecP[pid].m_EID );
            const DCR_TetSolidShape3::Triangle& td( p_dcr->m_vecT[tid] );
            //\todo CONSIDER using element_cache transformed V, as they MUST be in the cache for CT
            Vec3 vec_pos_0[3] = { tr * mal::GRange<1,3>( element_cache[eid].Bs_invBm() * mal::Concat(1,p_dcr->m_vecV[td.GetVID(0)] ) ),
                                  tr * mal::GRange<1,3>( element_cache[eid].Bs_invBm() * mal::Concat(1,p_dcr->m_vecV[td.GetVID(1)] ) ),
                                  tr * mal::GRange<1,3>( element_cache[eid].Bs_invBm() * mal::Concat(1,p_dcr->m_vecV[td.GetVID(2)] ) ) };
            //DEBUG
            if( bViz )
            {
                cd.m_VD.m_vecDCR2DCR_IB_vecClippedT_Side[side].push_back( util::make_triad(vec_pos_0[0],vec_pos_0[1],vec_pos_0[2]) );
            }

            // Gather all CTP-clipped edges on T and all (b)-edges
            struct clipped_edge_type
            {
                uint32 m_CTPID; bool m_bInwards; int m_EIT; Real m_Lambda; //Sort by [eit+lambda) = [0,3)
                inline clipped_edge_type( uint32 ctpid, bool b_inwards, int eit, Real lambda )
                : m_CTPID(ctpid), m_bInwards(b_inwards), m_EIT(eit), m_Lambda(lambda) {}
            };
            std::vector< clipped_edge_type > vecCE;
            std::vector< util::triad<Vec3,Vec3,uint32> > vecBoundarySegmentsCCW; //<p0,p1,ctp.m_CurveId>
            for( auto ctpid : ct.second )
            {
                const crossing_triangle_pair_type& ctp( vecCTP[ctpid] );
                const intersection_triangle3_triangle3_result_type& ittr( ctp.m_ITTR );
                if( (ittr.IsVF() && side == 0) || (ctp.m_ITTR.IsFV() && side == 1) )
                {
                    Vec3 other_normal_0( mal::Cross( ctp.m_vecPos_0[oside][1]-ctp.m_vecPos_0[oside][0], ctp.m_vecPos_0[oside][2]-ctp.m_vecPos_0[oside][0] ) );
                    bool bInwardsE0( mal::Dot( ctp.m_vecPos_0[side][(ittr.m_EdgeInT0+1)%3] - ctp.m_vecPos_0[side][ittr.m_EdgeInT0], other_normal_0 ) < 0 );
                    bool bInwardsE1( !bInwardsE0 ); //\todo Unless numerical things go wrong...
                    vecCE.push_back( clipped_edge_type( ctpid, bInwardsE0, ittr.m_EdgeInT0, ittr.m_Lambda0 ) );
                    vecCE.push_back( clipped_edge_type( ctpid, bInwardsE1, ittr.m_EdgeInT1, ittr.m_Lambda1 ) );
                }
                else if( ctp.m_ITTR.IsEE() )
                {
                    Vec3 other_normal_0( mal::Cross( ctp.m_vecPos_0[oside][1]-ctp.m_vecPos_0[oside][0], ctp.m_vecPos_0[oside][2]-ctp.m_vecPos_0[oside][0] ) );
                    if( side == 0 )
                    {
                        bool bInwardsE0( mal::Dot( ctp.m_vecPos_0[side][(ittr.m_EdgeInT0+1)%3] - ctp.m_vecPos_0[side][ittr.m_EdgeInT0], other_normal_0 ) < 0 );
                        vecCE.push_back( clipped_edge_type( ctpid, bInwardsE0, ittr.m_EdgeInT0, ittr.m_Lambda0 ) );
                    }
                    else //side == 1
                    {
                        bool bInwardsE1( mal::Dot( ctp.m_vecPos_0[side][(ittr.m_EdgeInT1+1)%3] - ctp.m_vecPos_0[side][ittr.m_EdgeInT1], other_normal_0 ) < 0 );
                        vecCE.push_back( clipped_edge_type( ctpid, bInwardsE1, ittr.m_EdgeInT1, ittr.m_Lambda1 ) );
                    }
                }
                // else: VF side1 or FV side0, no T-edge involved

                /* Add CTP in CCW orientation according to side
                   \note This is the only place where ittr orientation
                   needs to be taken into account, all posterior
                   code is side-independent
                */
                if( side == 0 ) vecBoundarySegmentsCCW.emplace_back( ittr.m_Point0, ittr.m_Point1, ctp.m_CurveId );
                else vecBoundarySegmentsCCW.emplace_back( ittr.m_Point1, ittr.m_Point0, ctp.m_CurveId );
            }
            uint32 num_boundary_segments_from_CTP( vecBoundarySegmentsCCW.size() );

            // if( bLog ) { GEO_LOG("CT[%u] has #CTP=%u, #CE=%u", tid, (uint32)ct.second.size(), (uint32)vecCE.size() ); }
            GEO_LOG_ERROR_IF( vecCE.size() % 2 != 0, "CT[%u] has #CE=%u != 2*k, cancelling FloodIB", tid, (uint32)vecCE.size() );
            if( vecCE.size() % 2 != 0 )
            {
                GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_FloodIB_CT_ms );
                return false;
            }

            // Sort clipped edges CCW along the CT perimeter
            std::sort( vecCE.begin(), vecCE.end(),
                       []( const clipped_edge_type& a, const clipped_edge_type& b )
                       { return a.m_EIT < b.m_EIT || (a.m_EIT == b.m_EIT && a.m_Lambda < b.m_Lambda); } );

            /* Iterate CCW over vecCE starting at the first inwards ce
               and gather internal sub-segments, adding them to
               vecBoundarySegmentsCCW and connecting to any
               potentially adjacent IC.
            */
            uint32 first_inwards_ceid(0);
            while( first_inwards_ceid < vecCE.size() && !vecCE[first_inwards_ceid].m_bInwards ) first_inwards_ceid++;
            if( first_inwards_ceid != vecCE.size() )
            {
                uint32 it_ce0( first_inwards_ceid );
                uint32 it_ce1( (first_inwards_ceid+1)%vecCE.size() );
                while( it_ce1 != first_inwards_ceid )
                {
                    clipped_edge_type ce0( vecCE[it_ce0] );
                    clipped_edge_type ce1( vecCE[it_ce1] );
                    Vec3 p0( (1-ce0.m_Lambda)*vec_pos_0[ce0.m_EIT] + ce0.m_Lambda*vec_pos_0[(ce0.m_EIT+1)%3] );
                    Vec3 p1( (1-ce1.m_Lambda)*vec_pos_0[ce1.m_EIT] + ce1.m_Lambda*vec_pos_0[(ce1.m_EIT+1)%3] );
                    // Only for internal sub-intervals
                    if( ce0.m_bInwards )
                    {
                        GEO_LOG_ERROR_IF( ce1.m_bInwards, "CT[%u] has consecutive inwards CE, cancelling FloodIB", tid );
                        if( ce1.m_bInwards )
                        {
                            GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_FloodIB_CT_ms );
                            return false;
                        }

                        // We'll associate the CE to ce0's IC/IB, but later we may merge it into a different root IC/IB
                        uint32 ibid0( vecCTP[ce0.m_CTPID].m_CurveId );

                        // Connect 3 potential T.V inbetween ce0 and ce1
                        if( ce0.m_EIT != ce1.m_EIT  ) //different edge
                        {
                            vecBoundarySegmentsCCW.emplace_back( p0, vec_pos_0[(ce0.m_EIT+1)%3], ibid0 );
                            p0 = vec_pos_0[(ce0.m_EIT+1)%3];
                            if( (ce0.m_EIT+1)%3 != ce1.m_EIT  )
                            {
                                vecBoundarySegmentsCCW.emplace_back( p0, vec_pos_0[(ce0.m_EIT+2)%3], ibid0 );
                                p0 = vec_pos_0[(ce0.m_EIT+2)%3];
                            }
                        }
                        else if( ce0.m_Lambda > ce1.m_Lambda ) //same edge, previous lambda ==> include all 3 vertices inbetween
                        {
                            vecBoundarySegmentsCCW.emplace_back( p0, vec_pos_0[(ce0.m_EIT+1)%3], ibid0 );
                            vecBoundarySegmentsCCW.emplace_back( vec_pos_0[(ce0.m_EIT+1)%3], vec_pos_0[(ce0.m_EIT+2)%3], ibid0 );
                            vecBoundarySegmentsCCW.emplace_back( vec_pos_0[(ce0.m_EIT+2)%3], vec_pos_0[ce0.m_EIT], ibid0 );
                            p0 = vec_pos_0[ce0.m_EIT];
                        }
                        // Connect p0->p1, where p0 may have jumped to vertices between ce0 and ce1
                        vecBoundarySegmentsCCW.emplace_back( p0, p1, ibid0 );

                        // Merge IB
                        uint32 ibid1( vecCTP[ce1.m_CTPID].m_CurveId );
                        if( ibid0 != ibid1 )
                        {
                            if( vecIB[ibid0].m_ParentId != vecIB[ibid1].m_ParentId )
                            {
                                uint32 ib_root_id0( IB_ComputeRootId( vecIB, ibid0 ) );
                                uint32 ib_root_id1( IB_ComputeRootId( vecIB, ibid1 ) );
                                if( ib_root_id0 != ib_root_id1 )
                                {
                                    vecIB[ib_root_id0].m_ParentId = ib_root_id1;
                                    vecIB[ibid0].m_ParentId = ib_root_id1; //path-compresion
                                    vecIB[ibid1].m_ParentId = ib_root_id1;
                                    if( bLog ) { GEO_LOG( "MF-Set connecting IB[%u] with root %u to root IB[%u] through CTP pair (%u,%u)",
                                                           ibid0, ib_root_id0, ib_root_id1, ce0.m_CTPID, ce1.m_CTPID ); }
                                }
                            }
                        }
                    }
                    it_ce0 = (it_ce0+1) % vecCE.size();
                    it_ce1 = (it_ce1+1) % vecCE.size();
                }
            }
            // else //this can happen for hole-only CT, it's NOT a problem at all
            // {
            //     if( bLog ) { GEO_LOG( "CT[%u] has no inwards CE", tid ); }
            // }

            Vec3 triangle_v_area_0( mal::Cross( vec_pos_0[1] - vec_pos_0[0], vec_pos_0[2] - vec_pos_0[0] ) );
            Vec3 triangle_normal_0( mal::SafeNormalized(triangle_v_area_0) );

            // Integrate over all regions at once by accumulating the effect of each CT-region CCW boundary segment on the IC/IB that originated it
#define __ENABLE_ACC_CT_ON_IB //This is the PROPER way, it reduces
                              //DRAMATICALLY glitches in CT-regions,
                              //but there are STILL a few empty IB
                              //with empty crossingP generated, we're
                              //just ignoring them after the
                              //FloodIB(), but should not be gen in
                              //the first place...
#ifdef __ENABLE_ACC_CT_ON_IB
            struct ib_act_type { Vec3 m_VectorArea_0; Vec3 m_Centroid_0; inline ib_act_type() : m_VectorArea_0(0), m_Centroid_0(0) {} };
            std::vector<ib_act_type> vecIB_AccCT( vecIB.size() );
#else
            Vec3 ct_v_area_0(0);
            Vec3 ct_centroid_0(0); //Computed wrt triangle vertex 0 \see https://en.wikipedia.org/wiki/Centroid
#endif
            for( auto segment : vecBoundarySegmentsCCW )
            {
                Vec3 va( mal::Cross( segment.first-vec_pos_0[0], segment.second-vec_pos_0[0] ) ); //segment vector-area contribution
                Real a( Real(0.5)*mal::Dot(va,triangle_normal_0) );
                Vec3 c( a*(segment.first-vec_pos_0[0]+segment.second-vec_pos_0[0]) ); //segment centroid contribution
#ifdef __ENABLE_ACC_CT_ON_IB
                // Accumulate on root IB
                uint32 ib_root_id( IB_ComputeRootId( vecIB, segment.third ) );
                vecIB_AccCT[ib_root_id].m_VectorArea_0 += va;
                vecIB_AccCT[ib_root_id].m_Centroid_0 += c;
#else
                ct_v_area_0 += va;
                ct_centroid_0 += c;
#endif
            }

#ifdef __ENABLE_ACC_CT_ON_IB
            for( unsigned int it_ib=0; it_ib<vecIB_AccCT.size(); it_ib++ )
            {
                ib_act_type& act( vecIB_AccCT[it_ib] );
                Real act_area( Real(0.5)*mal::Dot(act.m_VectorArea_0,triangle_normal_0) );
                if( act_area > 0 )
                {
                    // Finish IB-esclusive CT-region centroid computation
                    act.m_Centroid_0 *= mal::Rcp(3*act_area);
                    act.m_Centroid_0 += vec_pos_0[0];

                    // Accumulate contribution on IB.m_mapCrossingP, create if required
                    uint32 ib_root_id( IB_ComputeRootId( vecIB, it_ib ) );
                    dcr3_intersection_boundary_type& ib( vecIB[ib_root_id] );
                    dcr3_intersection_boundary_type::Patch* pPatch(0);
                    auto it_patch( ib.m_mapCrossingP.find(pid) );
                    if( it_patch == ib.m_mapCrossingP.end() )
                        pPatch = &ib.m_mapCrossingP.insert( std::make_pair( pid, dcr3_intersection_boundary_type::Patch(pid,ib_root_id) ) ).first->second;
                    else
                        pPatch = &it_patch->second;
                    pPatch->m_NumTriangles++;
                    pPatch->m_AvgPos_0 += act_area*act.m_Centroid_0;
                    pPatch->m_VectorArea_0 += act.m_VectorArea_0;
                    pPatch->m_AreaSq += mal::Sq(act_area);
                    pPatch->m_Area += act_area;
                    //DEBUG: Add only the segments for the current ACT, connect with act centroid to show proper association
                    if( bViz )
                    {
                        cd.m_VD.m_vecDCR2DCR_IB_vecCT_Side[side].emplace_back( act.m_Centroid_0, act.m_VectorArea_0 );
                        for( auto segment : vecBoundarySegmentsCCW )
                            if( segment.third == it_ib )
                                cd.m_VD.m_vecDCR2DCR_IB_vecCE_Side[side].emplace_back( segment.first, segment.second );
                        // This "radial" segments are useful but cause confusion with actual triangle edges, consider specific vector and color, if Viz is required
                        // for( auto segment : vecBoundarySegmentsCCW )
                        //     if( segment.third == it_ib )
                        //         cd.m_VD.m_vecDCR2DCR_IB_vecCE_Side[side].emplace_back( act.m_Centroid_0, segment.first );
                    }
                }
            }
#else
            Real ct_area( Real(0.5)*mal::Dot(ct_v_area_0,triangle_normal_0) );

            /* TODO:
               If area is negative we may be in case b.2), where the CT
               does not contribute any edge to vecBoundarySegmentsCCW and
               the existing ones define a positive negative hole, so that
               the whole triangle area needs to be added.
            */
            // if( area < 0 && num_boundary_segments_from_CTP == vecBoundarySegmentsCCW.size() )
            // {
            //     v_area_0 += triangle_v_area_0;
            //     //\todo FIX centroid, should add solid triangle contribution, including area...
            //     // Real a( Real(0.5)*mal::Dot(va,triangle_normal_0) );
            //     // centroid_0 += a * (vec_pos_0[0]-vec_pos_0[0] + vec_pos_0[1]-vec_pos_0[]
            //     area = Real(0.5)*mal::Dot(v_area_0,triangle_normal_0);
            // }
            GEO_LOG_ERROR_IF( ct_area < 0 && num_boundary_segments_from_CTP != vecBoundarySegmentsCCW.size(),
                              "CT[%u] has negative area %f but not due hole-only CT...",
                              tid, ct_area );
            if( ct_area > 1e-9 )
            {
                ct_centroid_0 *= mal::Rcp(3*ct_area) ;
                ct_centroid_0 += vec_pos_0[0];

                /*\todo NO! THIS IS WRONG... we MUST NOT add all CT
                 segment stuff to THE SAME IB/IC , because they can
                 come from DIFFERENT IC... each CT clipped interior
                 region CAN BE in a different IC or connect several.
                 => In effect, we should consider that the
                 accumulation patch can be different for each
                 CT-segment, however, we can cluster them by ICID, and
                 use the ICID->IB.  IF a CT merges two IC/IB, a later
                 merge pass will handle them, no need to do it here.
                 */

                // Accumulate CT magnitudes on parent IB crossing Patch, creating it if not present
                uint32 ibid( vecCTP[ct.second.front()].m_CurveId );
                uint32 ib_root_id( IB_ComputeRootId( vecIB, ibid ) );
                dcr3_intersection_boundary_type& ib( vecIB[ib_root_id] );
                dcr3_intersection_boundary_type::Patch* pPatch(0);
                auto it_patch( ib.m_mapCrossingP.find(pid) );
                if( it_patch == ib.m_mapCrossingP.end() )
                    pPatch = &ib.m_mapCrossingP.insert( std::make_pair( pid, dcr3_intersection_boundary_type::Patch(pid,ib_root_id) ) ).first->second;
                else
                    pPatch = &it_patch->second;
                pPatch->m_NumTriangles++;
                pPatch->m_AvgPos_0 += ct_area*ct_centroid_0;
                pPatch->m_VectorArea_0 += ct_v_area_0;
                pPatch->m_AreaSq += mal::Sq(ct_area);
                pPatch->m_Area += ct_area;

                //DEBUG
                if( bViz )
                {
                    cd.m_VD.m_vecDCR2DCR_IB_vecCT_Side[side].emplace_back( ct_centroid_0, ct_v_area_0 );
                    for( auto segment : vecBoundarySegmentsCCW )
                        cd.m_VD.m_vecDCR2DCR_IB_vecCE_Side[side].emplace_back( ct_centroid_0, segment.first );
                    for( auto segment : vecBoundarySegmentsCCW )
                        cd.m_VD.m_vecDCR2DCR_IB_vecCE_Side[side].emplace_back( segment.first, segment.second );
                }
            }
            else
            {
                GEO_LOG_ERROR( "CT[%u] has small or negative area %f, ignored", tid, ct_area );
            }
#endif //__ENABLE_ACC_CT_ON_IB
            num_ct++;
            num_ct_boundary_segments += vecBoundarySegmentsCCW.size();
        }
        if( bLog ) { GEO_LOG( "#CT = %u, #BoundarySegments = %u", num_ct, num_ct_boundary_segments ); }

        GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_FloodIB_CT_ms );
        GEO_NP_STAT_ADD( mig2015.m_Num_IB_CT, num_ct );
    }

    //---- Merge IB into their MF-Set roots
    // IB contains non-trivial members, so be careful with efficiency here...
    //IMPORTANT: IF merge only happens through CT when processing each IC, WE COULD AVOID merging alltogether..., by NOT CREATING the IB untill non-merge has been confirmed!
    for( unsigned int it_ib=0; it_ib<vecIB.size(); it_ib++ )
    {
        //Merge to parent, which is GUARANTEED to have a smaller index (but needs to be found using path compression)
        uint32 root_ibid( IB_ComputeRootId(vecIB,it_ib) );
        if( it_ib != root_ibid )
        {
            if( bLog ) { GEO_LOG( "Merging IB[%u] into RootIB[%u]", it_ib, root_ibid ); }
            vecIB[root_ibid].Merge( vecIB[it_ib] );
        }
    }
    //Swap-out non-root from vecIB
    unsigned int num_root_ib(0);
    for( unsigned int it_ib=0; it_ib<vecIB.size(); it_ib++ )
    {
        if( it_ib != vecIB[it_ib].m_ParentId )
        {
            //TEMP We DO NOT compact the IB to avoid the HIGH cost of swapping them to throw them away afterwards
            // std::swap(vecIB[it_ib],vecIB.back());
            // vecIB.pop_back();
        }
        else
            num_root_ib++;
    }

    /* Gather root IB \todo THIS SEEMS TO LEAVE all patches
       EMPTY... probably I need a copy constructor or assignment op or
       some shit...

    std::vector<dcr3_intersection_boundary_type> vecRootIB;
    for( unsigned int it_ib=0; it_ib<vecIB.size(); it_ib++ )
        if( it_ib != vecIB[it_ib].m_ParentId )
            vecRootIB.push_back( vecIB[it_ib] );
    // Save root IB IMPORTANT: m_ParentId are out of date now
    std::swap(vecIB,vecRootIB);
    vecRootIB.clear();
    */

    // Compute aggregate stuff for unique IB
    if( bLog ) { GEO_LOG( "#RootIB[%d] = %u", side, num_root_ib ); }
    for( unsigned int it_ib=0; it_ib<vecIB.size(); it_ib++ )
    {
        if( it_ib == vecIB[it_ib].m_ParentId ) //is root
        {
            dcr3_intersection_boundary_type& ib( vecIB[it_ib] );
            uint32 num_crossing_p( ib.m_mapCrossingP.size() );
            uint32 num_internal_p( ib.m_vecP.size() );

            // Append crossing P to IB.P
            for( const auto& patch : ib.m_mapCrossingP )
                ib.m_vecP.push_back( patch.second );
            ib.m_mapCrossingP.clear();

            // Accumulate per-patch data into unique IB
            for( auto& patch : ib.m_vecP )
            {
                ib.m_NumTriangles += patch.m_NumTriangles;
                ib.m_AvgPos_0 += patch.m_AvgPos_0;
                ib.m_VectorArea_0 += patch.m_VectorArea_0;
                ib.m_AreaSq += patch.m_AreaSq;
                ib.m_Area += patch.m_Area;
                if( patch.m_Area > 0 )
                    patch.m_AvgPos_0 = mal::Rcp(patch.m_Area) * patch.m_AvgPos_0;
                GEO_LOG_ERROR_IF( patch.m_NumTriangles == 0, "RootIB[%d][%u].P[%u] with 0 triangles!", side, it_ib, patch.m_PID );
            }
            if( ib.m_Area > 0 )
                ib.m_AvgPos_0 = mal::Rcp(ib.m_Area) * ib.m_AvgPos_0;

            /* If no FloodT/P, compute global IB vector area using CCW
               boundary segments only

               IMPORTANT: This seems to be correct and is A LOT
               CHEAPER than per-patch integration, thus consider using
               it exclusively. Flood will still be necessary to
               discover internal T/P and compute their area for
               penalty pressure distribution.
             */
            if( p_context->m_DCR2DCR_IB_Method == Context::eDCR2DCR_IB_None )
            {
                ib.m_VectorArea_0 = Vec3::Zero();
                for( auto icid : ib.m_vecICID )
                {
                    for( auto ctpid : vecIC[icid].m_vecCTPID )
                    {
                        const crossing_triangle_pair_type& ctp( vecCTP[ctpid] );
                        ib.m_VectorArea_0 += Real(0.5)*mal::Cross( ctp.m_ITTR.m_Point0-ib.m_AvgPos_0, ctp.m_ITTR.m_Point1-ib.m_AvgPos_0 );
                    }
                }
                if( side == 1 ) ib.m_VectorArea_0 = -ib.m_VectorArea_0;
            }

            if( bLog ) { GEO_LOG( "RootIB[%d][%u] with #IC = %u, #P = %u+%u, #T = %u",
                                   side, it_ib, (uint32)ib.m_vecICID.size(), num_internal_p, num_crossing_p, ib.m_NumTriangles ); }
            GEO_LOG_ERROR_IF( ib.m_NumTriangles == 0, "RootIB[%d][%u] with 0 triangles, may be due to CT mis-association to a single IC/IB!", side, it_ib );
        }
    }
    //DEBUG
    {
        for( const auto& ib : vecIB )
        {
            if( bViz )
            {
                // if( ib.m_NumTriangles > 0 ) TEMP: Draw them to debug visually by now
                cd.m_VD.m_vecDCR2DCR_IB_PosAndNormal_Side_x_IB[side].push_back( std::make_pair(ib.m_AvgPos_0,mal::SafeNormalized(ib.m_VectorArea_0)) );
                cd.m_VD.m_vecDCR2DCR_IB_vecP_Side_x_IB[side].push_back( std::vector< std::pair<Vec3,Vec3> >() );
                for( const auto& patch : ib.m_vecP )
                    // if( patch.m_NumTriangles > 0 ) TEMP: Draw them to debug visually by now
                        cd.m_VD.m_vecDCR2DCR_IB_vecP_Side_x_IB[side].back().push_back( std::make_pair(patch.m_AvgPos_0,mal::SafeNormalized(patch.m_VectorArea_0)) );
            }
        }
    }
    return true;
}

// TODO: Add RaycastDCR() that uses BVH.RayCast() internally, try to make it orthogonal to RC_Method
bool RaycastDCR( const TetSolidShape3* p_solid,
                 const DCR_TetSolidShape3* p_dcr, const BVH_TetSolidShape3* p_bvh,
                 const Transform3& tr, const Transform3& inv_tr,
                 const dcr3_element_cache_type& element_cache,
                 const Vec3& p, const Vec3& n, Real max_length,
                 RayHit3& rh,
                 const Context* p_context )
{
    RayHit3 tmp_rh;
    rh.m_Interval.Set(max_length);
    std::vector< BVH_TetSolidShape3::entry_index_type > vecOverlaps;
    Vec3 p_local( inv_tr*p );
    Vec3 n_local( inv_tr.Rot()*n );

    GEO_NP_BEGIN_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_RaycastBVH_ms );
    if( p_bvh->RayCast( p, n, Interval(0,max_length), vecOverlaps ) )
    {
        GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_RaycastBVH_ms );

        // GEO_LOG( "RaycastDCR overlaps BVH" );
        for( auto eid : vecOverlaps )
        {
            if( eid < p_dcr->m_NumElements && p_dcr->m_vecE[eid].m_NumPatches > 0 ) //ignore non-DCR elements
            {
                // GEO_LOG( "RaycastDCR overlaps E[%u].BV", eid );
                const dcr3_element_cache_type::ECD& ecd( element_cache[eid] );
                if( Overlap_Segment3_Tetrahedron3_BDOP( p_local, p_local+max_length*n_local, ecd.invBs(), p_dcr->m_vecE[eid].m_BDOP ) )
                {
                    // GEO_LOG( "RaycastDCR overlaps E[%u]", eid );
                    GEO_NP_BEGIN_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_RaycastDCR_ms );
                    if( RayCast_DCR_E_DoubleSided( p_local, n_local, Interval(0,rh.m_Interval.Min()), //clipped ray to optimize multiple rh candidate search
                                                   p_dcr, eid,
                                                   ecd.Bs(), ecd.invBs(), ecd.Bs_invBm(),
                                                   tmp_rh, 0, p_context ) ) //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min()
                    {
                        rh = tmp_rh; //\note This already contains rh data in global coords as well as feature_id=EID and barycentric coords wrt DCR.E
                        rh.m_Point = tr*rh.m_Point;
                        rh.m_Normal = tr.Rot()*rh.m_Normal;
                    }
                    GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_RaycastDCR_ms );
                }
            }
        }
    }
    else
    {
        GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_RaycastBVH_ms );
    }
    return rh.m_Interval.Min() < max_length;
}

/* AddCFP_Intersection_DCR3_E_DCR3_E_BruteForce()

   Add crossing_triangle_pair_type to given vec_cfp, found by
   bruteforce DCR.E vs DCR.E all Patch.Triangle pairwise tests.
*/
unsigned int AddCFP_Intersection_DCR3_E_DCR3_E_BruteForce( const TetSolidShape3* p_solid1, const Vec3* vec_sdof1, const DCR_TetSolidShape3* p_dcr1, uint32 eid1,
                                                           const TetSolidShape3* p_solid2, const Vec3* vec_sdof2, const DCR_TetSolidShape3* p_dcr2, uint32 eid2,
                                                           const Transform3& tr1_2, const Transform3& tr2_0,
                                                           const dcr3_element_cache_type& element_cache1, const dcr3_element_cache_type& element_cache2,
                                                           std::vector< crossing_triangle_pair_type >& vec_cfp,
                                                           const Context* p_context = g_pDefaultContext )
{
    const DCR_TetSolidShape3::Element& ed1( p_dcr1->m_vecE[eid1] );
    const DCR_TetSolidShape3::Element& ed2( p_dcr2->m_vecE[eid2] );

    /* This can happen if non-layer[0] elements are present in
     * the DCR and contain no geometry patches
     */
    if( ed1.m_NumPatches == 0 || ed2.m_NumPatches == 0 ) return 0;

    // Transform ALL vertices e1 \todo and cache
    util::ScopedAllocator scoped_allocator_e1( p_context->m_ScratchPad, "Contact_DCRTS3_DCRTS3_BruteForce_MSS_Lazy_DCR_e1" );
    Vec3* vec_v1_2 = scoped_allocator_e1.NewArrayPOD<Vec3>(ed1.m_NumVertices);
    {
#ifdef __ENABLE_DCR3_ELEMENT_CACHE
        const Mat4x4& Bs_invBm( element_cache1[eid1].Bs_invBm() );
        // Compute tr1_2_F and tr1_2_T so that v1_2 = tr1_2_F*r1 + tr1_2_T
        Mat3x3 tr1_2_F( tr1_2.Rot() * mal::GRange<1,1,3,3>(Bs_invBm) );
        Vec3 tr1_2_T( tr1_2 * mal::GRange<1,3>( mal::GColumn<0>(Bs_invBm) ) );
        // Transform ed1 m_NumVertices vertices from m_First and store them in 0-based indices inside vec_v1_2
        for( unsigned int it_vie=0; it_vie<ed1.m_NumVertices; it_vie++ )
            vec_v1_2[it_vie] = tr1_2_F * p_dcr1->m_vecV[ed1.m_FirstVID + it_vie] + tr1_2_T;
            //\todo SLOWER version: vec_v1_2[it_vie] = tr1_2 * mal::GRange<1,3>( Bs_invBm * mal::Concat(1,p_dcr1->m_vecV[ed1.m_FirstVID + it_vie]) );
#else
        uint32 vec_nid[4] = { p_solid1->T_VID(eid1,0), p_solid1->T_VID(eid1,1), p_solid1->T_VID(eid1,2), p_solid1->T_VID(eid1,3) };
        Mat3x3 B( mal::GMat3x3_From_Columns(vec_sdof1[vec_nid[1]]-vec_sdof1[vec_nid[0]],
                                            vec_sdof1[vec_nid[2]]-vec_sdof1[vec_nid[0]],
                                            vec_sdof1[vec_nid[3]]-vec_sdof1[vec_nid[0]]) );
        const Vec3* default_sdof1( p_solid1->GetVecDefaultSDOF() );
        Mat3x3 Bm( mal::GMat3x3_From_Columns(default_sdof1[vec_nid[1]]-default_sdof1[vec_nid[0]],
                                             default_sdof1[vec_nid[2]]-default_sdof1[vec_nid[0]],
                                             default_sdof1[vec_nid[3]]-default_sdof1[vec_nid[0]]) ); //\todo Could be precomputed in DCR::ED
        Mat3x3 B_invBm( B * mal::Inverse(Bm) );
        Transform3 tr1_2_B_invBm( tr1_2.m_Pos, tr1_2.m_Rot * B_invBm );
        Vec3 p0_2( tr1_2.m_Rot * vec_sdof1[vec_nid[0]] );
        // Transform ed1 m_NumVertices vertices from m_First and store them in 0-based indices inside vec_v1_2
        for( unsigned int it_vie=0; it_vie<ed1.m_NumVertices; it_vie++ )
            vec_v1_2[it_vie] = tr1_2_B_invBm * (p_dcr1->m_vecV[ed1.m_FirstVID + it_vie] - default_sdof1[vec_nid[0]]) + p0_2;
#endif
        GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed1.m_NumVertices );
    }
    // Transform ALL vertices e2 \todo and cache
    util::ScopedAllocator scoped_allocator_e2( p_context->m_ScratchPad, "Contact_DCRTS3_DCRTS3_BruteForce_MSS_Lazy_DCR_e2" );
    Vec3* vec_v2_2 = scoped_allocator_e2.NewArrayPOD<Vec3>(ed2.m_NumVertices);
    {
#ifdef __ENABLE_DCR3_ELEMENT_CACHE
        const Mat4x4& Bs_invBm( element_cache2[eid2].Bs_invBm() );
        // Compute F and T so that v2_2 = F*r2 + T
        Mat3x3 F( mal::GRange<1,1,3,3>(Bs_invBm) );
        Vec3 T( mal::GRange<1,3>( mal::GColumn<0>(Bs_invBm) ) );
        // Transform ed1 m_NumVertices vertices from m_First and store them in 0-based indices inside vec_v1_2
        for( unsigned int it_vie=0; it_vie<ed2.m_NumVertices; it_vie++ )
            vec_v2_2[it_vie] = F * p_dcr2->m_vecV[ed2.m_FirstVID + it_vie] + T;
            //\todo SLOWER version: vec_v2_2[it_vie] = mal::GRange<1,3>( Bs_invBm * mal::Concat(1,p_dcr2->m_vecV[ed2.m_FirstVID + it_vie]) );
#else
        uint32 vec_nid[4] = { p_solid2->T_VID(eid2,0), p_solid2->T_VID(eid2,1), p_solid2->T_VID(eid2,2), p_solid2->T_VID(eid2,3) };
        Mat3x3 B( mal::GMat3x3_From_Columns(vec_sdof2[vec_nid[1]]-vec_sdof2[vec_nid[0]],
                                            vec_sdof2[vec_nid[2]]-vec_sdof2[vec_nid[0]],
                                            vec_sdof2[vec_nid[3]]-vec_sdof2[vec_nid[0]]) );
        const Vec3* default_sdof2( p_solid2->GetVecDefaultSDOF() );
        Mat3x3 Bm( mal::GMat3x3_From_Columns(default_sdof2[vec_nid[1]]-default_sdof2[vec_nid[0]],
                                             default_sdof2[vec_nid[2]]-default_sdof2[vec_nid[0]],
                                             default_sdof2[vec_nid[3]]-default_sdof2[vec_nid[0]]) ); //\todo Could be precomputed in DCR::ED
        Mat3x3 B_invBm( B * mal::Inverse(Bm) );
        Vec3 p0_2( vec_sdof2[vec_nid[0]] );
        for( unsigned int it_vie=0; it_vie<ed2.m_NumVertices; it_vie++ )
            vec_v2_2[it_vie] = B_invBm * (p_dcr2->m_vecV[ed2.m_FirstVID + it_vie] - default_sdof2[vec_nid[0]]) + p0_2;
#endif
        GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed2.m_NumVertices );
    }
    //-- Test all patch pairs using cached transformed vertices
    // For each patch1 in e1
    uint32 num_crossings(0);
    for( unsigned int it_patch1=0; it_patch1 < ed1.m_NumPatches; it_patch1++ )
    {
        uint32 pid1( ed1.m_FirstPID + it_patch1 );
        const geo::DCR_TetSolidShape3::Patch& pd1( p_dcr1->m_vecP[pid1] );
        // For each patch2 in e2
        for( unsigned int it_patch2=0; it_patch2 < ed2.m_NumPatches; it_patch2++ )
        {
            uint32 pid2( ed2.m_FirstPID + it_patch2 );
            const geo::DCR_TetSolidShape3::Patch& pd2( p_dcr2->m_vecP[pid2] );
            // For each triangle (a1,b1,c1)
            for( unsigned int it_tid1=0; it_tid1 < pd1.m_NumTriangles; it_tid1++ )
            {
                unsigned int tid1( pd1.m_FirstTID + it_tid1 );
                Vec3 a1_2( vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(0) - ed1.m_FirstVID ] );
                Vec3 b1_2( vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(1) - ed1.m_FirstVID ] );
                Vec3 c1_2( vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(2) - ed1.m_FirstVID ] );
                //\todo IMPORTANT! Early-out if tid1 does NOT overlap bv2!! \todo use bslab2 here, it's FAST to test and provides the smallest volume among all BV types
                // for each triangle (a2,b2,c2)
                for( unsigned int it_tid2=0; it_tid2 < pd2.m_NumTriangles; it_tid2++ )
                {
                    unsigned int tid2( pd2.m_FirstTID + it_tid2 );
                    Vec3 a2_2( vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(0) - ed2.m_FirstVID ] );
                    Vec3 b2_2( vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(1) - ed2.m_FirstVID ] );
                    Vec3 c2_2( vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(2) - ed2.m_FirstVID ] );
                    intersection_triangle3_triangle3_result_type ittr;
                    if( Intersection_Triangle3_Triangle3( a1_2, b1_2, c1_2, a2_2, b2_2, c2_2, ittr ) )
                    {
                        // save Crossing Feature Pair
                        ittr.m_Point0 = tr2_0 * ittr.m_Point0;
                        ittr.m_Point1 = tr2_0 * ittr.m_Point1;
                        vec_cfp.push_back( crossing_triangle_pair_type( ittr,
                                                                        eid1, pid1, tid1,
                                                                        tr2_0*a1_2, tr2_0*b1_2, tr2_0*c1_2,
                                                                        eid2, pid2, tid2,
                                                                        tr2_0*a2_2, tr2_0*b2_2, tr2_0*c2_2 ) );
                        num_crossings++;
                    }
                }
            }
        }
    }
    return num_crossings;
}


/* Intersect DCR3 Element vs DCR3 Element using b-DOP-Tree
  \todo Subst std::vector stack with std::array
*/
unsigned int AddCFP_Intersection_DCR3_E_DCR3_E_BDT( const TetSolidShape3* p_solid1, const Vec3* vec_sdof1, const DCR_TetSolidShape3* p_dcr1, uint32 eid1,
                                                    const TetSolidShape3* p_solid2, const Vec3* vec_sdof2, const DCR_TetSolidShape3* p_dcr2, uint32 eid2,
                                                    const Transform3& tr1_2, const Transform3& tr2_0,
                                                    const dcr3_element_cache_type& element_cache1, const dcr3_element_cache_type& element_cache2,
                                                    std::vector< crossing_triangle_pair_type >& vec_cfp,
                                                    const Context* p_context = g_pDefaultContext )
{
    uint32 num_crossings(0);
    const DCR_TetSolidShape3::Element& ed1( p_dcr1->m_vecE[eid1] );
    const DCR_TetSolidShape3::Element& ed2( p_dcr2->m_vecE[eid2] );

    /* This can happen if non-layer[0] elements are present in
     * the DCR and contain no geometry patches
     */
    if( ed1.m_NumPatches == 0 || ed2.m_NumPatches == 0 ) return 0;

    uint32 vec_nid1[4] = { p_solid1->T_VID(eid1,0), p_solid1->T_VID(eid1,1), p_solid1->T_VID(eid1,2), p_solid1->T_VID(eid1,3) };
    uint32 vec_nid2[4] = { p_solid2->T_VID(eid2,0), p_solid2->T_VID(eid2,1), p_solid2->T_VID(eid2,2), p_solid2->T_VID(eid2,3) };

    Vec3* vec_v1_2(0);
    Vec3* vec_v2_2(0);
    util::ScopedAllocator scoped_allocator( p_context->m_ScratchPad, "AddCFP_Intersection_DCR3_E_DCR3_E_BDT" );
#ifdef __ENABLE_DCR3_ELEMENT_CACHE
    if( !p_context->m_DCR2DCR_E2E_BDT_CacheTransformedV )
    {
        vec_v1_2 = scoped_allocator.NewArrayPOD<Vec3>(ed1.m_NumVertices);
        vec_v2_2 = scoped_allocator.NewArrayPOD<Vec3>(ed2.m_NumVertices);
    }
#else
    vec_v1_2 = scoped_allocator.NewArrayPOD<Vec3>(ed1.m_NumVertices);
    vec_v2_2 = scoped_allocator.NewArrayPOD<Vec3>(ed2.m_NumVertices);
#endif

    Mat4x4 Bs_invBm_1,Bs_invBm_2; // Only used if !Pretransform
    Mat3x3 F1,F2; //Ds * mal::Inverse(Dm) );
    if( p_context->m_DCR2DCR_E2E_BDT_Pretransform )
    {
        // Transform ALL vertices e1
        {
#ifdef __ENABLE_DCR3_ELEMENT_CACHE
            //\todo THIS IS UGLY, we do a const cast, ALSO, scoped
            //alloc is BROKEN here because scope is declared inside
            //Pretransform conditional, WORKS because we don't
            //overwrite, but scoped mem is accessed OUTSIDE the scope
            //it's reserved for
            if( p_context->m_DCR2DCR_E2E_BDT_CacheTransformedV )
            {
                const Mat4x4& Bs_invBm( element_cache1[eid1].Bs_invBm() );
                F1 = mal::GRange<1,1,3,3>(Bs_invBm);
                if( element_cache1[eid1].m_FirstTransformedVID == 0xFFFFFFFF )
                    GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed1.m_NumVertices );
                vec_v1_2 = const_cast<Vec3*>( element_cache1.GetTransformedV( eid1, tr1_2 ) );
            }
            else
            {
                const Mat4x4& Bs_invBm( element_cache1[eid1].Bs_invBm() );
                F1 = mal::GRange<1,1,3,3>(Bs_invBm);
                Mat3x3 tr1_2_F( tr1_2.Rot() * F1 );
                Vec3 tr1_2_T( tr1_2 * mal::GRange<1,3>( mal::GColumn<0>(Bs_invBm) ) );
                // Transform ed1 m_NumVertices vertices from m_First and store them in 0-based indices inside vec_v1_2
                for( unsigned int it_vie=0; it_vie<ed1.m_NumVertices; it_vie++ )
                    vec_v1_2[it_vie] = tr1_2_F * p_dcr1->m_vecV[ed1.m_FirstVID + it_vie] + tr1_2_T;
                element_cache1[eid1].m_Stats_NumTransform++;
                GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed1.m_NumVertices );
            }
#else
            Mat3x3 Ds( mal::GMat3x3_From_Columns(vec_sdof1[vec_nid1[1]]-vec_sdof1[vec_nid1[0]],
                                                 vec_sdof1[vec_nid1[2]]-vec_sdof1[vec_nid1[0]],
                                                 vec_sdof1[vec_nid1[3]]-vec_sdof1[vec_nid1[0]]) );
            const Vec3* default_sdof1( p_solid1->GetVecDefaultSDOF() );
            Mat3x3 Dm( mal::GMat3x3_From_Columns(default_sdof1[vec_nid1[1]]-default_sdof1[vec_nid1[0]],
                                                 default_sdof1[vec_nid1[2]]-default_sdof1[vec_nid1[0]],
                                                 default_sdof1[vec_nid1[3]]-default_sdof1[vec_nid1[0]]) ); //\todo Could be precomputed in DCR::ED
            F1 = Ds * mal::Inverse(Dm);
            Transform3 tr1_2_F1( tr1_2.m_Pos, tr1_2.m_Rot * F1 );
            Vec3 p0_2( tr1_2.m_Rot * vec_sdof1[vec_nid1[0]] );
            // Transform ed1 m_NumVertices vertices from m_First and store them in 0-based indices inside vec_v1_2
            for( unsigned int it_vie=0; it_vie<ed1.m_NumVertices; it_vie++ )
                vec_v1_2[it_vie] = tr1_2_F1 * (p_dcr1->m_vecV[ed1.m_FirstVID + it_vie] - default_sdof1[vec_nid1[0]]) + p0_2;
            GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed1.m_NumVertices );
#endif
        }
        // Transform ALL vertices e2
        {
#ifdef __ENABLE_DCR3_ELEMENT_CACHE
            //\todo THIS IS UGLY, we do a const cast, ALSO, scoped
            //alloc is BROKEN here because scope is declared inside
            //Pretransform conditional, WORKS because we don't
            //overwrite, but scoped mem is accessed OUTSIDE the scope
            //it's reserved for
            if( p_context->m_DCR2DCR_E2E_BDT_CacheTransformedV )
            {
                const Mat4x4& Bs_invBm( element_cache2[eid2].Bs_invBm() );
                F2 = mal::GRange<1,1,3,3>(Bs_invBm);
                if( element_cache2[eid2].m_FirstTransformedVID == 0xFFFFFFFF )
                    GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed2.m_NumVertices );
                vec_v2_2 = const_cast<Vec3*>( element_cache2.GetTransformedV( eid2, Transform3::Identity() ) );
            }
            else
            {
                const Mat4x4& Bs_invBm( element_cache2[eid2].Bs_invBm() );
                F2 = mal::GRange<1,1,3,3>(Bs_invBm);
                Vec3 T( mal::GRange<1,3>( mal::GColumn<0>(Bs_invBm) ) );
                // Transform ed1 m_NumVertices vertices from m_First and store them in 0-based indices inside vec_v1_2
                for( unsigned int it_vie=0; it_vie<ed2.m_NumVertices; it_vie++ )
                    vec_v2_2[it_vie] = F2 * p_dcr2->m_vecV[ed2.m_FirstVID + it_vie] + T;
                element_cache2[eid2].m_Stats_NumTransform++;
                GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed2.m_NumVertices );
            }
#else

            Mat3x3 B( mal::GMat3x3_From_Columns(vec_sdof2[vec_nid2[1]]-vec_sdof2[vec_nid2[0]],
                                                vec_sdof2[vec_nid2[2]]-vec_sdof2[vec_nid2[0]],
                                                vec_sdof2[vec_nid2[3]]-vec_sdof2[vec_nid2[0]]) );
            const Vec3* default_sdof2( p_solid2->GetVecDefaultSDOF() );
            Mat3x3 Dm( mal::GMat3x3_From_Columns(default_sdof2[vec_nid2[1]]-default_sdof2[vec_nid2[0]],
                                                 default_sdof2[vec_nid2[2]]-default_sdof2[vec_nid2[0]],
                                                 default_sdof2[vec_nid2[3]]-default_sdof2[vec_nid2[0]]) ); //\todo Could be precomputed in DCR::ED
            F2 = B * mal::Inverse(Dm);
            Vec3 p0_2( vec_sdof2[vec_nid2[0]] );
            for( unsigned int it_vie=0; it_vie<ed2.m_NumVertices; it_vie++ )
                vec_v2_2[it_vie] = F2 * (p_dcr2->m_vecV[ed2.m_FirstVID + it_vie] - default_sdof2[vec_nid2[0]]) + p0_2;
            GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed2.m_NumVertices );
#endif
        }
    }
    else
    {
#ifdef __ENABLE_DCR3_ELEMENT_CACHE
        Bs_invBm_1 = element_cache1[eid1].Bs_invBm();
        Bs_invBm_2 = element_cache2[eid2].Bs_invBm();
#else
        //\todo this is SLOW, try to reuse F1, F2 instead!
        const Vec3* default_sdof1( p_solid1->GetVecDefaultSDOF() );
        const Vec3* default_sdof2( p_solid2->GetVecDefaultSDOF() );
        Bs_invBm_1 = Mat4x4( 1, 1, 1, 1,
                             vec_sdof1[vec_nid1[0]].x(), vec_sdof1[vec_nid1[1]].x(), vec_sdof1[vec_nid1[2]].x(), vec_sdof1[vec_nid1[3]].x(),
                             vec_sdof1[vec_nid1[0]].y(), vec_sdof1[vec_nid1[1]].y(), vec_sdof1[vec_nid1[2]].y(), vec_sdof1[vec_nid1[3]].y(),
                             vec_sdof1[vec_nid1[0]].z(), vec_sdof1[vec_nid1[1]].z(), vec_sdof1[vec_nid1[2]].z(), vec_sdof1[vec_nid1[3]].z() )
                     * mal::Inverse( Mat4x4( 1, 1, 1, 1,
                                             default_sdof1[vec_nid1[0]].x(), default_sdof1[vec_nid1[1]].x(), default_sdof1[vec_nid1[2]].x(), default_sdof1[vec_nid1[3]].x(),
                                             default_sdof1[vec_nid1[0]].y(), default_sdof1[vec_nid1[1]].y(), default_sdof1[vec_nid1[2]].y(), default_sdof1[vec_nid1[3]].y(),
                                             default_sdof1[vec_nid1[0]].z(), default_sdof1[vec_nid1[1]].z(), default_sdof1[vec_nid1[2]].z(), default_sdof1[vec_nid1[3]].z() ) ); //\todo COULD BE PRECOMPUTED in DCR.E
        Bs_invBm_2 = Mat4x4( 1, 1, 1, 1,
                             vec_sdof2[vec_nid2[0]].x(), vec_sdof2[vec_nid2[1]].x(), vec_sdof2[vec_nid2[2]].x(), vec_sdof2[vec_nid2[3]].x(),
                             vec_sdof2[vec_nid2[0]].y(), vec_sdof2[vec_nid2[1]].y(), vec_sdof2[vec_nid2[2]].y(), vec_sdof2[vec_nid2[3]].y(),
                             vec_sdof2[vec_nid2[0]].z(), vec_sdof2[vec_nid2[1]].z(), vec_sdof2[vec_nid2[2]].z(), vec_sdof2[vec_nid2[3]].z() )
                     * mal::Inverse( Mat4x4( 1, 1, 1, 1,
                                             default_sdof2[vec_nid2[0]].x(), default_sdof2[vec_nid2[1]].x(), default_sdof2[vec_nid2[2]].x(), default_sdof2[vec_nid2[3]].x(),
                                             default_sdof2[vec_nid2[0]].y(), default_sdof2[vec_nid2[1]].y(), default_sdof2[vec_nid2[2]].y(), default_sdof2[vec_nid2[3]].y(),
                                             default_sdof2[vec_nid2[0]].z(), default_sdof2[vec_nid2[1]].z(), default_sdof2[vec_nid2[2]].z(), default_sdof2[vec_nid2[3]].z() ) ); //\todo COULD BE PRECOMPUTED in DCR.E
#endif
        F1 = mal::GRange<1,1,3,3>(Bs_invBm_1);
        F2 = mal::GRange<1,1,3,3>(Bs_invBm_2);
    }

    // Precomp nodes in 2-refsys for BV(BDOP)
    Vec3 vec_node_pos1_2[4] = { tr1_2*vec_sdof1[vec_nid1[0]], tr1_2*vec_sdof1[vec_nid1[1]], tr1_2*vec_sdof1[vec_nid1[2]], tr1_2*vec_sdof1[vec_nid1[3]] };
    Vec3 vec_node_pos2_2[4] = { vec_sdof2[vec_nid2[0]], vec_sdof2[vec_nid2[1]], vec_sdof2[vec_nid2[2]], vec_sdof2[vec_nid2[3]] };
    // Presort node indices along each BV axis (\note For a GKDOP<D,K> each array must allocate K/2 sub-arrays arrays of D node-in-element indices)
    int vec_sorted_nie1[12] = { 0,1,2,3, 0,1,2,3, 0,1,2,3 }; //x [0,4) y [4,8) z [8,12)
    int vec_sorted_nie2[12] = { 0,1,2,3, 0,1,2,3, 0,1,2,3 }; //x [0,4) y [4,8) z [8,12)
    bv::Compute_AABB3_From_BDOP4_PreSort( vec_node_pos1_2, vec_sorted_nie1 );
    bv::Compute_AABB3_From_BDOP4_PreSort( vec_node_pos2_2, vec_sorted_nie2 );

#ifdef __ENABLE_DCR3_BDT_EXTRA_AXIS_ELEMENT
    Vec3 axis_2(0);
    // Compute axis for the "flattest" DCR.E (the DCR.E with less DCR.P is expected to be the "flattest")
    if( ed1.m_NumPatches <= ed2.m_NumPatches )
    {
        for( unsigned int it_pie1=0; it_pie1<ed1.m_NumPatches; it_pie1++ )
            axis_2 += p_dcr1->m_vecP[ed1.m_FirstPID+it_pie1].m_VectorArea;
        axis_2 = tr1_2.m_Rot * (mal::Adjugate(F1) * axis_2); //\sum P1.N affine-transformed by F1^-T and rotated to solid2 ref sys
    }
    else
    {
        for( unsigned int it_pie2=0;it_pie2<ed2.m_NumPatches; it_pie2++ )
            axis_2 += p_dcr2->m_vecP[ed2.m_FirstPID+it_pie2].m_VectorArea;
        axis_2 = mal::Adjugate(F2) * axis_2; //\sum P2.N affine-transformed by F2^-T
    }
    /* Compute both DCR.E aggregate axis and choose "flattest" (for
       equal scalar areas, the flattest patch has the largest vector
       area)
       \note WORSE than ed1.m_NumPatches <= ed2.m_NumPatches criteria
    {
        Vec3 axis1_1(0);
        for( unsigned int it_pie1=0; it_pie1<ed1.m_NumPatches; it_pie1++ )
            axis1_1 += p_dcr1->m_vecP[ed1.m_FirstPID+it_pie1].m_VectorArea;
        axis1_1 = mal::Adjugate(F1) * axis1_1; //\sum P1.N affine-transformed by F1^-T and rotated to solid2 ref sys
        Vec3 axis2_2(0);
        for( unsigned int it_pie2=0;it_pie2<ed2.m_NumPatches; it_pie2++ )
            axis2_2 += p_dcr2->m_vecP[ed2.m_FirstPID+it_pie2].m_VectorArea;
        axis2_2 = mal::Adjugate(F2) * axis2_2; //\sum P2.N affine-transformed by F2^-T
        axis_2 = mal::NormSq(axis1_1) > mal::NormSq(axis2_2) ? tr1_2.m_Rot * axis1_1 : axis2_2;
    }
    */
    int vec_sorted_nie_axis1[4] = { 0,1,2,3 };
    int vec_sorted_nie_axis2[4] = { 0,1,2,3 };
    Real vec_node_proj_axis1[4] = { mal::Dot(vec_node_pos1_2[0],axis_2), mal::Dot(vec_node_pos1_2[1],axis_2), mal::Dot(vec_node_pos1_2[2],axis_2), mal::Dot(vec_node_pos1_2[3],axis_2) };
    Real vec_node_proj_axis2[4] = { mal::Dot(vec_node_pos2_2[0],axis_2), mal::Dot(vec_node_pos2_2[1],axis_2), mal::Dot(vec_node_pos2_2[2],axis_2), mal::Dot(vec_node_pos2_2[3],axis_2) };
    //We use Sort4_Bubbble, simpler code and slightly slower than "optimal" Sort4_HC
    for( int i=0; i<3; i++ )
        for( int j=0; j<3; j++ )
        {
            if( vec_node_proj_axis1[vec_sorted_nie_axis1[j]] > vec_node_proj_axis1[vec_sorted_nie_axis1[j+1]] )
                std::swap(vec_sorted_nie_axis1[j],vec_sorted_nie_axis1[j+1]);
            if( vec_node_proj_axis2[vec_sorted_nie_axis2[j]] > vec_node_proj_axis2[vec_sorted_nie_axis2[j+1]] )
                std::swap(vec_sorted_nie_axis2[j],vec_sorted_nie_axis2[j+1]);
        }
#endif

#ifdef __ENABLE_DCR3_BDT_EXTRA_AXIS_DOP3K14
    int vec_sorted_nie_extra1[4][4] = { { 0,1,2,3 }, { 0,1,2,3 }, { 0,1,2,3 }, { 0,1,2,3 } };
    int vec_sorted_nie_extra2[4][4] = { { 0,1,2,3 }, { 0,1,2,3 }, { 0,1,2,3 }, { 0,1,2,3 } };
    Real vec_node_proj_extra1[4][4] = { { bv::GKDOP_Dot<3,14,3>(vec_node_pos1_2[0]), bv::GKDOP_Dot<3,14,3>(vec_node_pos1_2[1]), bv::GKDOP_Dot<3,14,3>(vec_node_pos1_2[2]), bv::GKDOP_Dot<3,14,3>(vec_node_pos1_2[3]) },
                                        { bv::GKDOP_Dot<3,14,4>(vec_node_pos1_2[0]), bv::GKDOP_Dot<3,14,4>(vec_node_pos1_2[1]), bv::GKDOP_Dot<3,14,4>(vec_node_pos1_2[2]), bv::GKDOP_Dot<3,14,4>(vec_node_pos1_2[3]) },
                                        { bv::GKDOP_Dot<3,14,5>(vec_node_pos1_2[0]), bv::GKDOP_Dot<3,14,5>(vec_node_pos1_2[1]), bv::GKDOP_Dot<3,14,5>(vec_node_pos1_2[2]), bv::GKDOP_Dot<3,14,5>(vec_node_pos1_2[3]) },
                                        { bv::GKDOP_Dot<3,14,6>(vec_node_pos1_2[0]), bv::GKDOP_Dot<3,14,6>(vec_node_pos1_2[1]), bv::GKDOP_Dot<3,14,6>(vec_node_pos1_2[2]), bv::GKDOP_Dot<3,14,6>(vec_node_pos1_2[3]) } };
    Real vec_node_proj_extra2[4][4] = { { bv::GKDOP_Dot<3,14,3>(vec_node_pos2_2[0]), bv::GKDOP_Dot<3,14,3>(vec_node_pos2_2[1]), bv::GKDOP_Dot<3,14,3>(vec_node_pos2_2[2]), bv::GKDOP_Dot<3,14,3>(vec_node_pos2_2[3]) },
                                        { bv::GKDOP_Dot<3,14,4>(vec_node_pos2_2[0]), bv::GKDOP_Dot<3,14,4>(vec_node_pos2_2[1]), bv::GKDOP_Dot<3,14,4>(vec_node_pos2_2[2]), bv::GKDOP_Dot<3,14,4>(vec_node_pos2_2[3]) },
                                        { bv::GKDOP_Dot<3,14,5>(vec_node_pos2_2[0]), bv::GKDOP_Dot<3,14,5>(vec_node_pos2_2[1]), bv::GKDOP_Dot<3,14,5>(vec_node_pos2_2[2]), bv::GKDOP_Dot<3,14,5>(vec_node_pos2_2[3]) },
                                        { bv::GKDOP_Dot<3,14,6>(vec_node_pos2_2[0]), bv::GKDOP_Dot<3,14,6>(vec_node_pos2_2[1]), bv::GKDOP_Dot<3,14,6>(vec_node_pos2_2[2]), bv::GKDOP_Dot<3,14,6>(vec_node_pos2_2[3]) } };
    //We use Sort4_Bubbble, simpler code and slightly slower than "optimal" Sort4_HC
    for( int a=0; a<4; a++ ) //axis
        for( int i=0; i<3; i++ ) //bubble iter
            for( int j=0; j<3; j++ ) //bubble iter
            {
                if( vec_node_proj_extra1[a][vec_sorted_nie_extra1[a][j]] > vec_node_proj_extra1[a][vec_sorted_nie_extra1[a][j+1]] )
                    std::swap(vec_sorted_nie_extra1[a][j],vec_sorted_nie_extra1[a][j+1]);
                if( vec_node_proj_extra2[a][vec_sorted_nie_extra2[a][j]] > vec_node_proj_extra2[a][vec_sorted_nie_extra2[a][j+1]] )
                    std::swap(vec_sorted_nie_extra2[a][j],vec_sorted_nie_extra2[a][j+1]);
            }
#endif

    // For all patch pairs
    for( unsigned int it_pie1=0; it_pie1<ed1.m_NumPatches; it_pie1++ )
    {
        uint32 pid1( ed1.m_FirstPID+it_pie1 );
        const DCR_TetSolidShape3::Patch& pd1( p_dcr1->m_vecP[pid1] );

#ifdef __ENABLE_DCR3_BDT_EXTRA_AXIS_PATCH
        Vec3 axis_2( tr1_2.m_Rot * mal::Adjugate(F1) * pd1.m_VectorArea ); //P.N affine-transformed by Bs*Bm^-1 and rotated to solid2 ref sys
        int vec_sorted_nie_axis1[4] = { 0,1,2,3 };
        int vec_sorted_nie_axis2[4] = { 0,1,2,3 };
        Real vec_node_proj_axis1[4] = { mal::Dot(vec_node_pos1_2[0],axis_2), mal::Dot(vec_node_pos1_2[1],axis_2), mal::Dot(vec_node_pos1_2[2],axis_2), mal::Dot(vec_node_pos1_2[3],axis_2) };
        Real vec_node_proj_axis2[4] = { mal::Dot(vec_node_pos2_2[0],axis_2), mal::Dot(vec_node_pos2_2[1],axis_2), mal::Dot(vec_node_pos2_2[2],axis_2), mal::Dot(vec_node_pos2_2[3],axis_2) };
        //We use Sort4_Bubbble, simpler code and slightly slower than "optimal" Sort4_HC
        for( int i=0; i<3; i++ )
            for( int j=0; j<3; j++ )
            {
                if( vec_node_proj_axis1[vec_sorted_nie_axis1[j]] > vec_node_proj_axis1[vec_sorted_nie_axis1[j+1]] )
                    std::swap(vec_sorted_nie_axis1[j],vec_sorted_nie_axis1[j+1]);
                if( vec_node_proj_axis2[vec_sorted_nie_axis2[j]] > vec_node_proj_axis2[vec_sorted_nie_axis2[j+1]] )
                    std::swap(vec_sorted_nie_axis2[j],vec_sorted_nie_axis2[j+1]);
            }
#endif

        for( unsigned int it_pie2=0; it_pie2<ed2.m_NumPatches; it_pie2++ )
        {
            uint32 pid2( ed2.m_FirstPID+it_pie2 );
            const DCR_TetSolidShape3::Patch& pd2( p_dcr2->m_vecP[pid2] );
            // Init stack with pd1,pd2 root BDTN
            typedef std::pair< DCR_TetSolidShape3::Patch::bdt_node_tip_range,
                               DCR_TetSolidShape3::Patch::bdt_node_tip_range > stack_entry_type;
            std::vector< stack_entry_type > stackBDTN; //\todo USE scope-alloc stack with fixed maxsize    \todo Subst std::vector stack with std::array
            stackBDTN.push_back( stack_entry_type(DCR_TetSolidShape3::Patch::bdt_node_tip_range(0,pd1.m_NumTriangles),
                                                  DCR_TetSolidShape3::Patch::bdt_node_tip_range(0,pd2.m_NumTriangles) ) );
            while( !stackBDTN.empty() )
            {
                // Pop PS subarrays
                stack_entry_type se( stackBDTN.back() );
                stackBDTN.pop_back();
                const DCR_TetSolidShape3::Patch::bdt_node_tip_range node_tr1( se.first );
                const DCR_TetSolidShape3::Patch::bdt_node_tip_range node_tr2( se.second );
                const uint32 length_tr1( node_tr1.second - node_tr1.first );
                const uint32 length_tr2( node_tr2.second - node_tr2.first );
                if( length_tr1 * length_tr2 <= p_context->m_DCR2DCR_E2E_BDT_MaxLeafTests ) //use bruteforce below a given number of pairwise tests
                {
                    // Test all pairwise triangles in subtrees \todo Consider inverting loop if length_tr2 > length_tr1, OR CONSIDER caching transformed triangles to avoid performing length_tr1*length_tr2 => length_tr1 + length_tr2
                    if( !p_context->m_DCR2DCR_E2E_BDT_Pretransform )
                        GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, 3*(length_tr1 + length_tr1*length_tr2) );

                    for( unsigned int it_tin1=0; it_tin1<length_tr1; it_tin1++ )
                    {
                        uint32 tid1( pd1.m_FirstTID + node_tr1.first + it_tin1 );
                        const Vec3 a1_2( p_context->m_DCR2DCR_E2E_BDT_Pretransform
                                         ? vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(0) - ed1.m_FirstVID ]
                                         : tr1_2 * mal::GRange<1,3>(Bs_invBm_1 * mal::Concat(1, p_dcr1->m_vecV[ p_dcr1->m_vecT[tid1].GetVID(0) ])) );
                        const Vec3 b1_2( p_context->m_DCR2DCR_E2E_BDT_Pretransform
                                         ? vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(1) - ed1.m_FirstVID ]
                                         : tr1_2 * mal::GRange<1,3>(Bs_invBm_1 * mal::Concat(1, p_dcr1->m_vecV[ p_dcr1->m_vecT[tid1].GetVID(1) ])) );
                        const Vec3 c1_2( p_context->m_DCR2DCR_E2E_BDT_Pretransform
                                         ? vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(2) - ed1.m_FirstVID ]
                                         : tr1_2 * mal::GRange<1,3>(Bs_invBm_1 * mal::Concat(1, p_dcr1->m_vecV[ p_dcr1->m_vecT[tid1].GetVID(2) ])) );
                        for( unsigned int it_tin2=0; it_tin2<length_tr2; it_tin2++ )
                        {
                            uint32 tid2( pd2.m_FirstTID + node_tr2.first + it_tin2 );
                            const Vec3 a2_2( p_context->m_DCR2DCR_E2E_BDT_Pretransform
                                             ? vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(0) - ed2.m_FirstVID ]
                                             : mal::GRange<1,3>(Bs_invBm_2 * mal::Concat(1, p_dcr2->m_vecV[ p_dcr2->m_vecT[tid2].GetVID(0) ])) );
                            const Vec3 b2_2( p_context->m_DCR2DCR_E2E_BDT_Pretransform
                                             ? vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(1) - ed2.m_FirstVID ]
                                             : mal::GRange<1,3>(Bs_invBm_2 * mal::Concat(1, p_dcr2->m_vecV[ p_dcr2->m_vecT[tid2].GetVID(1) ])) );
                            const Vec3 c2_2( p_context->m_DCR2DCR_E2E_BDT_Pretransform
                                             ? vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(2) - ed2.m_FirstVID ]
                                             : mal::GRange<1,3>(Bs_invBm_2 * mal::Concat(1, p_dcr2->m_vecV[ p_dcr2->m_vecT[tid2].GetVID(2) ])) );
                            intersection_triangle3_triangle3_result_type ittr;
                            if( Intersection_Triangle3_Triangle3( a1_2, b1_2, c1_2, a2_2, b2_2, c2_2, ittr ) )
                            {
                                // save Crossing Feature Pair
                                ittr.m_Point0 = tr2_0 * ittr.m_Point0;
                                ittr.m_Point1 = tr2_0 * ittr.m_Point1;
                                vec_cfp.push_back( crossing_triangle_pair_type( ittr,
                                                                                eid1, pid1, tid1,
                                                                                tr2_0*a1_2, tr2_0*b1_2, tr2_0*c1_2,
                                                                                eid2, pid2, tid2,
                                                                                tr2_0*a2_2, tr2_0*b2_2, tr2_0*c2_2 ) );
                                num_crossings++;
                            }
                        }
                    }
                }
                else
                {
                    // Get node/triangle data
                    const DCR_TetSolidShape3::Triangle& bdtn1( p_dcr1->m_vecT[pd1.m_FirstTID+node_tr1.first] );
                    const DCR_TetSolidShape3::Triangle& bdtn2( p_dcr2->m_vecT[pd2.m_FirstTID+node_tr2.first] );
                    /* IMPORTANT!! WE MUST handle the length_tr1==1
                       and length_tr2==1 cases explicitly, because
                       nodes have a BDOP that bounds their whole
                       subtree, which can be MUCH LARGER than their
                       individual BDOP and will NOT cull overlaps at
                       all. This is a consequence of storing
                       nodes-as-triangles. Instead, we compute the
                       global BV of the deformed triangle.
                    */
                    bool bOverlap(false);
                    bool bDescend1(false);
                    if( p_context->m_DCR2DCR_E2E_BDT_Split3
                        && length_tr1 == 1 )
                    {
                        uint32 tid1( pd1.m_FirstTID + node_tr1.first );
                        Vec3 a1_2( p_context->m_DCR2DCR_E2E_BDT_Pretransform
                                   ? vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(0) - ed1.m_FirstVID ]
                                   : tr1_2 * mal::GRange<1,3>(Bs_invBm_1 * mal::Concat(1, p_dcr1->m_vecV[ p_dcr1->m_vecT[tid1].GetVID(0) ])) );
                        Vec3 b1_2( p_context->m_DCR2DCR_E2E_BDT_Pretransform
                                   ? vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(1) - ed1.m_FirstVID ]
                                   : tr1_2 * mal::GRange<1,3>(Bs_invBm_1 * mal::Concat(1, p_dcr1->m_vecV[ p_dcr1->m_vecT[tid1].GetVID(1) ])) );
                        Vec3 c1_2( p_context->m_DCR2DCR_E2E_BDT_Pretransform
                                   ? vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(2) - ed1.m_FirstVID ]
                                   : tr1_2 * mal::GRange<1,3>(Bs_invBm_1 * mal::Concat(1, p_dcr1->m_vecV[ p_dcr1->m_vecT[tid1].GetVID(2) ])) );
                        bv::AABB3 aabb1_2( bv::AABB3(a1_2).Merge(b1_2).Merge(c1_2) );
                        bv::BDOP4 bdop2_2( bdtn2.BDTN_BDOPq() );
                        bv::AABB3 aabb2_2( bv::Compute_AABB3_From_BDOP4_PreSorted( bdop2_2, vec_node_pos2_2, vec_sorted_nie2 ) );
                        bOverlap = bv::TestOverlap( aabb1_2, aabb2_2 );
                        bDescend1 = false;
                    }
                    else if( p_context->m_DCR2DCR_E2E_BDT_Split3
                             && length_tr2 == 1 )
                    {
                        uint32 tid2( pd2.m_FirstTID + node_tr2.first );
                        Vec3 a2_2( p_context->m_DCR2DCR_E2E_BDT_Pretransform
                                   ? vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(0) - ed2.m_FirstVID ]
                                   : mal::GRange<1,3>(Bs_invBm_2 * mal::Concat(1, p_dcr2->m_vecV[ p_dcr2->m_vecT[tid2].GetVID(0) ])) );
                        Vec3 b2_2( p_context->m_DCR2DCR_E2E_BDT_Pretransform
                                   ? vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(1) - ed2.m_FirstVID ]
                                   : mal::GRange<1,3>(Bs_invBm_2 * mal::Concat(1, p_dcr2->m_vecV[ p_dcr2->m_vecT[tid2].GetVID(1) ])) );
                        Vec3 c2_2( p_context->m_DCR2DCR_E2E_BDT_Pretransform
                                   ? vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(2) - ed2.m_FirstVID ]
                                   : mal::GRange<1,3>(Bs_invBm_2 * mal::Concat(1, p_dcr2->m_vecV[ p_dcr2->m_vecT[tid2].GetVID(2) ])) );
                        bv::AABB3 aabb2_2( bv::AABB3(a2_2).Merge(b2_2).Merge(c2_2) );
                        bv::BDOP4 bdop1_1( bdtn1.BDTN_BDOPq() );
                        bv::AABB3 aabb1_2( bv::Compute_AABB3_From_BDOP4_PreSorted( bdop1_1, vec_node_pos1_2, vec_sorted_nie1 ) );
                        bOverlap = bv::TestOverlap( aabb1_2, aabb2_2 );
                        bDescend1 = true;
                    }
                    else //both nodes are non-leaf
                    {
                        // Get cartesian BVs from BDOPs
#ifdef __ENABLE_TETSS_DCR_PATCH_BDT_UNPACKED //TEMP: NOT faster, disabled in cVersion=2 to reduce mem/cache usage
                        const bv::BDOP4 bdop1_1( p_context->m_DCR2DCR_E2E_BDT_Unpacked ? bdtn1.BDTN_BDOPq_UNPACKED() : bdtn1.BDTN_BDOPq() );
                        const bv::BDOP4 bdop2_2( p_context->m_DCR2DCR_E2E_BDT_Unpacked ? bdtn2.BDTN_BDOPq_UNPACKED() : bdtn2.BDTN_BDOPq() );
#else
                        const bv::BDOP4 bdop1_1( bdtn1.BDTN_BDOPq() );
                        const bv::BDOP4 bdop2_2( bdtn2.BDTN_BDOPq() );
#endif
                        Real vol1,vol2;
                        if (
#ifdef __ENABLE_DCR3_BDT_EXTRA_AXIS
                            // Early out, avoids 5-10% of BDOP4_BDOP4 full test and decreases BDOP-Tree false positives (??%)
                            bv::TestOverlap_BDOP4_BDOP4_Axis( bdop1_1, vec_node_proj_axis1, vec_sorted_nie_axis1,
                                                              bdop2_2, vec_node_proj_axis2, vec_sorted_nie_axis2 )
                            &&
#endif
#ifdef __ENABLE_DCR3_BDT_EXTRA_AXIS_DOP3K14
                            bv::TestOverlap_BDOP4_BDOP4_Axis( bdop1_1, vec_node_proj_extra1[0], vec_sorted_nie_extra1[0],
                                                              bdop2_2, vec_node_proj_extra2[0], vec_sorted_nie_extra2[0] )
                            && bv::TestOverlap_BDOP4_BDOP4_Axis( bdop1_1, vec_node_proj_extra1[1], vec_sorted_nie_extra1[1],
                                                                 bdop2_2, vec_node_proj_extra2[1], vec_sorted_nie_extra2[1] )
                            && bv::TestOverlap_BDOP4_BDOP4_Axis( bdop1_1, vec_node_proj_extra1[2], vec_sorted_nie_extra1[2],
                                                                 bdop2_2, vec_node_proj_extra2[2], vec_sorted_nie_extra2[2] )
                            && bv::TestOverlap_BDOP4_BDOP4_Axis( bdop1_1, vec_node_proj_extra1[3], vec_sorted_nie_extra1[3],
                                                                 bdop2_2, vec_node_proj_extra2[3], vec_sorted_nie_extra2[3] )
                            &&
#endif
                            bv::TestOverlap_BDOP4_BDOP4_ComputeVolume( bdop1_1, vec_node_pos1_2, vec_sorted_nie1,
                                                                       bdop2_2, vec_node_pos2_2, vec_sorted_nie2,
                                                                       vol1, vol2 ) )
                        {
                            bOverlap = true;
                            bDescend1 = ( length_tr2 == 1 //\todo If Split3, this is no longer possible
                                          || (!p_context->m_DCR2DCR_E2E_BDT_DescendLarger && length_tr1 > 1 && length_tr1 > length_tr2)
                                          || (p_context->m_DCR2DCR_E2E_BDT_DescendLarger && length_tr1 > 1 && vol1 >= vol2 ) );
                        }
                    }

                    // Test BDOP overlap and compute volume if true
                    if( bOverlap )
                    {
                        // Recurse \todo I'M SURE this logic can be
                        // simplified taking into account the previous
                        // length_tr1*length_tr2<maxleaftests branch
                        // (MAYBE req maxleaftests>0, which is
                        // REASONABLE)
                        if( bDescend1 )
                        {
                            // Recourse BDTN1.L/R + BDTN.S[0]
                            stackBDTN.push_back( stack_entry_type( DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr1.first, node_tr1.first+1 ),
                                                                   node_tr2 ) );
                            int remaining_tr1( length_tr1 - 1 );
                            if( remaining_tr1 > 1 )
                            {

                                stackBDTN.push_back( stack_entry_type( DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr1.first+1, node_tr1.first+1+remaining_tr1/2 ),
                                                                       node_tr2 ) );
                                stackBDTN.push_back( stack_entry_type( DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr1.first+1+remaining_tr1/2, node_tr1.second ),
                                                                       node_tr2 ) );
                            }
                            else if( remaining_tr1 == 1 )
                                stackBDTN.push_back( stack_entry_type( DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr1.first+1, node_tr1.first+2 ),
                                                                       node_tr2 ) );
                        }
                        else
                        {
                            // Recourse BDTN2.L/R + BDTN2.S[0]
                            stackBDTN.push_back( stack_entry_type( node_tr1,
                                                                   DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr2.first, node_tr2.first+1 ) ) );
                            int remaining_tr2( length_tr2 - 1 );
                            if( remaining_tr2 > 1 )
                            {
                                stackBDTN.push_back( stack_entry_type( node_tr1,
                                                                       DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr2.first+1, node_tr2.first+1+remaining_tr2/2 ) ) );
                                stackBDTN.push_back( stack_entry_type( node_tr1,
                                                                       DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr2.first+1+remaining_tr2/2, node_tr2.second ) ) );
                            }
                            else if( remaining_tr2 == 1 )
                                stackBDTN.push_back( stack_entry_type( node_tr1,
                                                                       DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr2.first+1, node_tr2.first+2 ) ) );
                        }
                    }
                }
            }
        }
    }
    return num_crossings;
}

#ifdef __ENABLE_SIMD_PRETRANSFORM
/* Intersect DCR3 Element vs DCR3 Element using b-DOP-Tree
   - SIMD optimizations
   \todo Subst std::vector stack with std::array
*/
unsigned int AddCFP_Intersection_DCR3_E_DCR3_E_BDT_SIMD( const TetSolidShape3* p_solid1, const Vec3* vec_sdof1, const DCR_TetSolidShape3* p_dcr1, uint32 eid1,
                                                         const TetSolidShape3* p_solid2, const Vec3* vec_sdof2, const DCR_TetSolidShape3* p_dcr2, uint32 eid2,
                                                         const Transform3& tr1_2, const Transform3& tr2_0,
                                                         const dcr3_element_cache_type& element_cache1, const dcr3_element_cache_type& element_cache2,
                                                         std::vector< crossing_triangle_pair_type >& vec_cfp,
                                                         const Context* p_context = g_pDefaultContext )
{
    uint32 num_crossings(0);
    const DCR_TetSolidShape3::Element& ed1( p_dcr1->m_vecE[eid1] );
    const DCR_TetSolidShape3::Element& ed2( p_dcr2->m_vecE[eid2] );

    /* This can happen if non-layer[0] elements are present in
     * the DCR and contain no geometry patches
     */
    if( ed1.m_NumPatches == 0 || ed2.m_NumPatches == 0 ) return 0;

    uint32 vec_nid1[4] = { p_solid1->T_VID(eid1,0), p_solid1->T_VID(eid1,1), p_solid1->T_VID(eid1,2), p_solid1->T_VID(eid1,3) };
    uint32 vec_nid2[4] = { p_solid2->T_VID(eid2,0), p_solid2->T_VID(eid2,1), p_solid2->T_VID(eid2,2), p_solid2->T_VID(eid2,3) };

    // Transform ALL vertices e1
    simd::V3f* vec_v1_2(0);
    {

        if( element_cache1[eid1].m_FirstTransformedVID == 0xFFFFFFFF )
            GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed1.m_NumVertices );
        vec_v1_2 = const_cast<simd::V3f*>( element_cache1.GetTransformedV_SIMD( eid1, tr1_2 ) );
    }
    // Transform ALL vertices e2
    simd::V3f* vec_v2_2(0);
    {
        if( element_cache2[eid2].m_FirstTransformedVID == 0xFFFFFFFF )
            GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed2.m_NumVertices );
        vec_v2_2 = const_cast<simd::V3f*>( element_cache2.GetTransformedV_SIMD( eid2, Transform3::Identity() ) );
    }

    const Mat3x3 F1( mal::GRange<1,1,3,3>( element_cache1[eid1].Bs_invBm() ) );
    const Mat3x3 F2( mal::GRange<1,1,3,3>( element_cache2[eid2].Bs_invBm() ) );

    // Precomp nodes in 2-refsys for BV(BDOP)
    Vec3 vec_node_pos1_2[4] = { tr1_2*vec_sdof1[vec_nid1[0]], tr1_2*vec_sdof1[vec_nid1[1]], tr1_2*vec_sdof1[vec_nid1[2]], tr1_2*vec_sdof1[vec_nid1[3]] };
    Vec3 vec_node_pos2_2[4] = { vec_sdof2[vec_nid2[0]], vec_sdof2[vec_nid2[1]], vec_sdof2[vec_nid2[2]], vec_sdof2[vec_nid2[3]] };
    // Presort node indices along each BV axis (\note For a GKDOP<D,K> each array must allocate K/2 sub-arrays arrays of D node-in-element indices)
    int vec_sorted_nie1[12] = { 0,1,2,3, 0,1,2,3, 0,1,2,3 }; //x [0,4) y [4,8) z [8,12)
    int vec_sorted_nie2[12] = { 0,1,2,3, 0,1,2,3, 0,1,2,3 }; //x [0,4) y [4,8) z [8,12)
    bv::Compute_AABB3_From_BDOP4_PreSort( vec_node_pos1_2, vec_sorted_nie1 );
    bv::Compute_AABB3_From_BDOP4_PreSort( vec_node_pos2_2, vec_sorted_nie2 );

#ifdef __ENABLE_DCR3_BDT_EXTRA_AXIS_ELEMENT
    Vec3 axis_2(0);
    // Compute axis for the "flattest" DCR.E (the DCR.E with less DCR.P is expected to be the "flattest")
    if( ed1.m_NumPatches <= ed2.m_NumPatches )
    {
        for( unsigned int it_pie1=0; it_pie1<ed1.m_NumPatches; it_pie1++ )
            axis_2 += p_dcr1->m_vecP[ed1.m_FirstPID+it_pie1].m_VectorArea;
        axis_2 = tr1_2.m_Rot * (mal::Adjugate(F1) * axis_2); //\sum P1.N affine-transformed by F1^-T and rotated to solid2 ref sys
    }
    else
    {
        for( unsigned int it_pie2=0;it_pie2<ed2.m_NumPatches; it_pie2++ )
            axis_2 += p_dcr2->m_vecP[ed2.m_FirstPID+it_pie2].m_VectorArea;
        axis_2 = mal::Adjugate(F2) * axis_2; //\sum P2.N affine-transformed by F2^-T
    }
    int vec_sorted_nie_axis1[4] = { 0,1,2,3 };
    int vec_sorted_nie_axis2[4] = { 0,1,2,3 };
    Real vec_node_proj_axis1[4] = { mal::Dot(vec_node_pos1_2[0],axis_2), mal::Dot(vec_node_pos1_2[1],axis_2), mal::Dot(vec_node_pos1_2[2],axis_2), mal::Dot(vec_node_pos1_2[3],axis_2) };
    Real vec_node_proj_axis2[4] = { mal::Dot(vec_node_pos2_2[0],axis_2), mal::Dot(vec_node_pos2_2[1],axis_2), mal::Dot(vec_node_pos2_2[2],axis_2), mal::Dot(vec_node_pos2_2[3],axis_2) };
    //We use Sort4_Bubbble, simpler code and slightly slower than "optimal" Sort4_HC
    for( int i=0; i<3; i++ )
        for( int j=0; j<3; j++ )
        {
            if( vec_node_proj_axis1[vec_sorted_nie_axis1[j]] > vec_node_proj_axis1[vec_sorted_nie_axis1[j+1]] )
                std::swap(vec_sorted_nie_axis1[j],vec_sorted_nie_axis1[j+1]);
            if( vec_node_proj_axis2[vec_sorted_nie_axis2[j]] > vec_node_proj_axis2[vec_sorted_nie_axis2[j+1]] )
                std::swap(vec_sorted_nie_axis2[j],vec_sorted_nie_axis2[j+1]);
        }
#endif

#ifdef __ENABLE_DCR3_BDT_EXTRA_AXIS_DOP3K14
    int vec_sorted_nie_extra1[4][4] = { { 0,1,2,3 }, { 0,1,2,3 }, { 0,1,2,3 }, { 0,1,2,3 } };
    int vec_sorted_nie_extra2[4][4] = { { 0,1,2,3 }, { 0,1,2,3 }, { 0,1,2,3 }, { 0,1,2,3 } };
    Real vec_node_proj_extra1[4][4] = { { bv::GKDOP_Dot<3,14,3>(vec_node_pos1_2[0]), bv::GKDOP_Dot<3,14,3>(vec_node_pos1_2[1]), bv::GKDOP_Dot<3,14,3>(vec_node_pos1_2[2]), bv::GKDOP_Dot<3,14,3>(vec_node_pos1_2[3]) },
                                        { bv::GKDOP_Dot<3,14,4>(vec_node_pos1_2[0]), bv::GKDOP_Dot<3,14,4>(vec_node_pos1_2[1]), bv::GKDOP_Dot<3,14,4>(vec_node_pos1_2[2]), bv::GKDOP_Dot<3,14,4>(vec_node_pos1_2[3]) },
                                        { bv::GKDOP_Dot<3,14,5>(vec_node_pos1_2[0]), bv::GKDOP_Dot<3,14,5>(vec_node_pos1_2[1]), bv::GKDOP_Dot<3,14,5>(vec_node_pos1_2[2]), bv::GKDOP_Dot<3,14,5>(vec_node_pos1_2[3]) },
                                        { bv::GKDOP_Dot<3,14,6>(vec_node_pos1_2[0]), bv::GKDOP_Dot<3,14,6>(vec_node_pos1_2[1]), bv::GKDOP_Dot<3,14,6>(vec_node_pos1_2[2]), bv::GKDOP_Dot<3,14,6>(vec_node_pos1_2[3]) } };
    Real vec_node_proj_extra2[4][4] = { { bv::GKDOP_Dot<3,14,3>(vec_node_pos2_2[0]), bv::GKDOP_Dot<3,14,3>(vec_node_pos2_2[1]), bv::GKDOP_Dot<3,14,3>(vec_node_pos2_2[2]), bv::GKDOP_Dot<3,14,3>(vec_node_pos2_2[3]) },
                                        { bv::GKDOP_Dot<3,14,4>(vec_node_pos2_2[0]), bv::GKDOP_Dot<3,14,4>(vec_node_pos2_2[1]), bv::GKDOP_Dot<3,14,4>(vec_node_pos2_2[2]), bv::GKDOP_Dot<3,14,4>(vec_node_pos2_2[3]) },
                                        { bv::GKDOP_Dot<3,14,5>(vec_node_pos2_2[0]), bv::GKDOP_Dot<3,14,5>(vec_node_pos2_2[1]), bv::GKDOP_Dot<3,14,5>(vec_node_pos2_2[2]), bv::GKDOP_Dot<3,14,5>(vec_node_pos2_2[3]) },
                                        { bv::GKDOP_Dot<3,14,6>(vec_node_pos2_2[0]), bv::GKDOP_Dot<3,14,6>(vec_node_pos2_2[1]), bv::GKDOP_Dot<3,14,6>(vec_node_pos2_2[2]), bv::GKDOP_Dot<3,14,6>(vec_node_pos2_2[3]) } };
    //We use Sort4_Bubbble, simpler code and slightly slower than "optimal" Sort4_HC
    for( int a=0; a<4; a++ ) //axis
        for( int i=0; i<3; i++ ) //bubble iter
            for( int j=0; j<3; j++ ) //bubble iter
            {
                if( vec_node_proj_extra1[a][vec_sorted_nie_extra1[a][j]] > vec_node_proj_extra1[a][vec_sorted_nie_extra1[a][j+1]] )
                    std::swap(vec_sorted_nie_extra1[a][j],vec_sorted_nie_extra1[a][j+1]);
                if( vec_node_proj_extra2[a][vec_sorted_nie_extra2[a][j]] > vec_node_proj_extra2[a][vec_sorted_nie_extra2[a][j+1]] )
                    std::swap(vec_sorted_nie_extra2[a][j],vec_sorted_nie_extra2[a][j+1]);
            }
#endif

    // For all patch pairs
    for( unsigned int it_pie1=0; it_pie1<ed1.m_NumPatches; it_pie1++ )
    {
        uint32 pid1( ed1.m_FirstPID+it_pie1 );
        const DCR_TetSolidShape3::Patch& pd1( p_dcr1->m_vecP[pid1] );

#ifdef __ENABLE_DCR3_BDT_EXTRA_AXIS_PATCH
        Vec3 axis_2( tr1_2.m_Rot * mal::Adjugate(F1) * pd1.m_VectorArea ); //P.N affine-transformed by Bs*Bm^-1 and rotated to solid2 ref sys
        int vec_sorted_nie_axis1[4] = { 0,1,2,3 };
        int vec_sorted_nie_axis2[4] = { 0,1,2,3 };
        Real vec_node_proj_axis1[4] = { mal::Dot(vec_node_pos1_2[0],axis_2), mal::Dot(vec_node_pos1_2[1],axis_2), mal::Dot(vec_node_pos1_2[2],axis_2), mal::Dot(vec_node_pos1_2[3],axis_2) };
        Real vec_node_proj_axis2[4] = { mal::Dot(vec_node_pos2_2[0],axis_2), mal::Dot(vec_node_pos2_2[1],axis_2), mal::Dot(vec_node_pos2_2[2],axis_2), mal::Dot(vec_node_pos2_2[3],axis_2) };
        //We use Sort4_Bubbble, simpler code and slightly slower than "optimal" Sort4_HC
        for( int i=0; i<3; i++ )
            for( int j=0; j<3; j++ )
            {
                if( vec_node_proj_axis1[vec_sorted_nie_axis1[j]] > vec_node_proj_axis1[vec_sorted_nie_axis1[j+1]] )
                    std::swap(vec_sorted_nie_axis1[j],vec_sorted_nie_axis1[j+1]);
                if( vec_node_proj_axis2[vec_sorted_nie_axis2[j]] > vec_node_proj_axis2[vec_sorted_nie_axis2[j+1]] )
                    std::swap(vec_sorted_nie_axis2[j],vec_sorted_nie_axis2[j+1]);
            }
#endif

        for( unsigned int it_pie2=0; it_pie2<ed2.m_NumPatches; it_pie2++ )
        {
            uint32 pid2( ed2.m_FirstPID+it_pie2 );
            const DCR_TetSolidShape3::Patch& pd2( p_dcr2->m_vecP[pid2] );
            // Init stack with pd1,pd2 root BDTN
            typedef std::pair< DCR_TetSolidShape3::Patch::bdt_node_tip_range,
                               DCR_TetSolidShape3::Patch::bdt_node_tip_range > stack_entry_type;
            std::vector< stack_entry_type > stackBDTN; //\todo USE scope-alloc stack with fixed maxsize    \todo Subst std::vector stack with std::array
            stackBDTN.push_back( stack_entry_type(DCR_TetSolidShape3::Patch::bdt_node_tip_range(0,pd1.m_NumTriangles),
                                                  DCR_TetSolidShape3::Patch::bdt_node_tip_range(0,pd2.m_NumTriangles) ) );
            while( !stackBDTN.empty() )
            {
                // Pop PS subarrays
                stack_entry_type se( stackBDTN.back() );
                stackBDTN.pop_back();
                const DCR_TetSolidShape3::Patch::bdt_node_tip_range node_tr1( se.first );
                const DCR_TetSolidShape3::Patch::bdt_node_tip_range node_tr2( se.second );
                const uint32 length_tr1( node_tr1.second - node_tr1.first );
                const uint32 length_tr2( node_tr2.second - node_tr2.first );
                if( length_tr1 * length_tr2 <= p_context->m_DCR2DCR_E2E_BDT_MaxLeafTests ) //use bruteforce below a given number of pairwise tests
                {
                    // Test all pairwise triangles in subtrees \todo Consider inverting loop if length_tr2 > length_tr1, OR CONSIDER caching transformed triangles to avoid performing length_tr1*length_tr2 => length_tr1 + length_tr2
                    if( !p_context->m_DCR2DCR_E2E_BDT_Pretransform )
                        GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, 3*(length_tr1 + length_tr1*length_tr2) );

                    for( unsigned int it_tin1=0; it_tin1<length_tr1; it_tin1++ )
                    {
                        uint32 tid1( pd1.m_FirstTID + node_tr1.first + it_tin1 );
                        const simd::V3f a1_2( vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(0) - ed1.m_FirstVID ] );
                        const simd::V3f b1_2( vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(1) - ed1.m_FirstVID ] );
                        const simd::V3f c1_2( vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(2) - ed1.m_FirstVID ] );
                        for( unsigned int it_tin2=0; it_tin2<length_tr2; it_tin2++ )
                        {
                            uint32 tid2( pd2.m_FirstTID + node_tr2.first + it_tin2 );
                            const simd::V3f a2_2( vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(0) - ed2.m_FirstVID ] );
                            const simd::V3f b2_2( vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(1) - ed2.m_FirstVID ] );
                            const simd::V3f c2_2( vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(2) - ed2.m_FirstVID ] );
                            intersection_triangle3_triangle3_result_type ittr;
                            if( Intersection_Triangle3_Triangle3_SIMD( a1_2, b1_2, c1_2, a2_2, b2_2, c2_2, ittr ) )
                            {
                                // save Crossing Feature Pair
                                ittr.m_Point0 = tr2_0 * ittr.m_Point0;
                                ittr.m_Point1 = tr2_0 * ittr.m_Point1;
                                vec_cfp.push_back( crossing_triangle_pair_type( ittr,
                                                                                eid1, pid1, tid1,
                                                                                tr2_0*a1_2.AsVec3f(), tr2_0*b1_2.AsVec3f(), tr2_0*c1_2.AsVec3f(),
                                                                                eid2, pid2, tid2,
                                                                                tr2_0*a2_2.AsVec3f(), tr2_0*b2_2.AsVec3f(), tr2_0*c2_2.AsVec3f() ) );
                                num_crossings++;
                            }
                        }
                    }
                }
                else
                {
                    // Get node/triangle data
                    const DCR_TetSolidShape3::Triangle& bdtn1( p_dcr1->m_vecT[pd1.m_FirstTID+node_tr1.first] );
                    const DCR_TetSolidShape3::Triangle& bdtn2( p_dcr2->m_vecT[pd2.m_FirstTID+node_tr2.first] );
                    /* IMPORTANT!! WE MUST handle the length_tr1==1
                       and length_tr2==1 cases explicitly, because
                       nodes have a BDOP that bounds their whole
                       subtree, which can be MUCH LARGER than their
                       individual BDOP and will NOT cull overlaps at
                       all. This is a consequence of storing
                       nodes-as-triangles. Instead, we compute the
                       global BV of the deformed triangle.
                    */
                    bool bOverlap(false);
                    bool bDescend1(false);
                    if( p_context->m_DCR2DCR_E2E_BDT_Split3
                        && length_tr1 == 1 )
                    {
                        uint32 tid1( pd1.m_FirstTID + node_tr1.first );
                        const Vec3 a1_2( vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(0) - ed1.m_FirstVID ].AsVec3f() );
                        const Vec3 b1_2( vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(1) - ed1.m_FirstVID ].AsVec3f() );
                        const Vec3 c1_2( vec_v1_2[ p_dcr1->m_vecT[tid1].GetVID(2) - ed1.m_FirstVID ].AsVec3f() );
                        bv::AABB3 aabb1_2( bv::AABB3(a1_2).Merge(b1_2).Merge(c1_2) );
                        bv::BDOP4 bdop2_2( bdtn2.BDTN_BDOPq() );
                        bv::AABB3 aabb2_2( bv::Compute_AABB3_From_BDOP4_PreSorted( bdop2_2, vec_node_pos2_2, vec_sorted_nie2 ) );
                        bOverlap = bv::TestOverlap( aabb1_2, aabb2_2 );
                        bDescend1 = false;
                    }
                    else if( p_context->m_DCR2DCR_E2E_BDT_Split3
                             && length_tr2 == 1 )
                    {
                        uint32 tid2( pd2.m_FirstTID + node_tr2.first );
                        const Vec3 a2_2( vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(0) - ed2.m_FirstVID ].AsVec3f() );
                        const Vec3 b2_2( vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(1) - ed2.m_FirstVID ].AsVec3f() );
                        const Vec3 c2_2( vec_v2_2[ p_dcr2->m_vecT[tid2].GetVID(2) - ed2.m_FirstVID ].AsVec3f() );
                        bv::AABB3 aabb2_2( bv::AABB3(a2_2).Merge(b2_2).Merge(c2_2) );
                        bv::BDOP4 bdop1_1( bdtn1.BDTN_BDOPq() );
                        bv::AABB3 aabb1_2( bv::Compute_AABB3_From_BDOP4_PreSorted( bdop1_1, vec_node_pos1_2, vec_sorted_nie1 ) );
                        bOverlap = bv::TestOverlap( aabb1_2, aabb2_2 );
                        bDescend1 = true;
                    }
                    else //both nodes are non-leaf
                    {
                        // Get cartesian BVs from BDOPs
#ifdef __ENABLE_TETSS_DCR_PATCH_BDT_UNPACKED //TEMP: NOT faster, disabled in cVersion=2 to reduce mem/cache usage
                        const bv::BDOP4 bdop1_1( p_context->m_DCR2DCR_E2E_BDT_Unpacked ? bdtn1.BDTN_BDOPq_UNPACKED() : bdtn1.BDTN_BDOPq() );
                        const bv::BDOP4 bdop2_2( p_context->m_DCR2DCR_E2E_BDT_Unpacked ? bdtn2.BDTN_BDOPq_UNPACKED() : bdtn2.BDTN_BDOPq() );
#else
                        const bv::BDOP4 bdop1_1( bdtn1.BDTN_BDOPq() );
                        const bv::BDOP4 bdop2_2( bdtn2.BDTN_BDOPq() );
#endif
                        Real vol1,vol2;
                        if (
#ifdef __ENABLE_DCR3_BDT_EXTRA_AXIS
                            // Early out, avoids 5-10% of BDOP4_BDOP4 full test and decreases BDOP-Tree false positives (??%)
                            bv::TestOverlap_BDOP4_BDOP4_Axis( bdop1_1, vec_node_proj_axis1, vec_sorted_nie_axis1,
                                                              bdop2_2, vec_node_proj_axis2, vec_sorted_nie_axis2 )
                            &&
#endif
#ifdef __ENABLE_DCR3_BDT_EXTRA_AXIS_DOP3K14
                            bv::TestOverlap_BDOP4_BDOP4_Axis( bdop1_1, vec_node_proj_extra1[0], vec_sorted_nie_extra1[0],
                                                              bdop2_2, vec_node_proj_extra2[0], vec_sorted_nie_extra2[0] )
                            && bv::TestOverlap_BDOP4_BDOP4_Axis( bdop1_1, vec_node_proj_extra1[1], vec_sorted_nie_extra1[1],
                                                                 bdop2_2, vec_node_proj_extra2[1], vec_sorted_nie_extra2[1] )
                            && bv::TestOverlap_BDOP4_BDOP4_Axis( bdop1_1, vec_node_proj_extra1[2], vec_sorted_nie_extra1[2],
                                                                 bdop2_2, vec_node_proj_extra2[2], vec_sorted_nie_extra2[2] )
                            && bv::TestOverlap_BDOP4_BDOP4_Axis( bdop1_1, vec_node_proj_extra1[3], vec_sorted_nie_extra1[3],
                                                                 bdop2_2, vec_node_proj_extra2[3], vec_sorted_nie_extra2[3] )
                            &&
#endif
                            bv::TestOverlap_BDOP4_BDOP4_ComputeVolume( bdop1_1, vec_node_pos1_2, vec_sorted_nie1,
                                                                       bdop2_2, vec_node_pos2_2, vec_sorted_nie2,
                                                                       vol1, vol2 ) )
                        {
                            bOverlap = true;
                            bDescend1 = ( length_tr2 == 1 //\todo If Split3, this is no longer possible
                                          || (!p_context->m_DCR2DCR_E2E_BDT_DescendLarger && length_tr1 > 1 && length_tr1 > length_tr2)
                                          || (p_context->m_DCR2DCR_E2E_BDT_DescendLarger && length_tr1 > 1 && vol1 >= vol2 ) );
                        }
                    }

                    // Test BDOP overlap and compute volume if true
                    if( bOverlap )
                    {
                        // Recurse \todo I'M SURE this logic can be
                        // simplified taking into account the previous
                        // length_tr1*length_tr2<maxleaftests branch
                        // (MAYBE req maxleaftests>0, which is
                        // REASONABLE)
                        if( bDescend1 )
                        {
                            // Recourse BDTN1.L/R + BDTN.S[0]
                            stackBDTN.push_back( stack_entry_type( DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr1.first, node_tr1.first+1 ),
                                                                   node_tr2 ) );
                            int remaining_tr1( length_tr1 - 1 );
                            if( remaining_tr1 > 1 )
                            {

                                stackBDTN.push_back( stack_entry_type( DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr1.first+1, node_tr1.first+1+remaining_tr1/2 ),
                                                                       node_tr2 ) );
                                stackBDTN.push_back( stack_entry_type( DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr1.first+1+remaining_tr1/2, node_tr1.second ),
                                                                       node_tr2 ) );
                            }
                            else if( remaining_tr1 == 1 )
                                stackBDTN.push_back( stack_entry_type( DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr1.first+1, node_tr1.first+2 ),
                                                                       node_tr2 ) );
                        }
                        else
                        {
                            // Recourse BDTN2.L/R + BDTN2.S[0]
                            stackBDTN.push_back( stack_entry_type( node_tr1,
                                                                   DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr2.first, node_tr2.first+1 ) ) );
                            int remaining_tr2( length_tr2 - 1 );
                            if( remaining_tr2 > 1 )
                            {
                                stackBDTN.push_back( stack_entry_type( node_tr1,
                                                                       DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr2.first+1, node_tr2.first+1+remaining_tr2/2 ) ) );
                                stackBDTN.push_back( stack_entry_type( node_tr1,
                                                                       DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr2.first+1+remaining_tr2/2, node_tr2.second ) ) );
                            }
                            else if( remaining_tr2 == 1 )
                                stackBDTN.push_back( stack_entry_type( node_tr1,
                                                                       DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr2.first+1, node_tr2.first+2 ) ) );
                        }
                    }
                }
            }
        }
    }
    return num_crossings;
}
#endif //__ENABLE_SIMD_PRETRANSFORM

unsigned int AddCFP_Intersection_DCR3_E_DCR3_E( const TetSolidShape3* p_solid1, const Vec3* vec_sdof1, const DCR_TetSolidShape3* p_dcr1, uint32 eid1,
                                                const TetSolidShape3* p_solid2, const Vec3* vec_sdof2, const DCR_TetSolidShape3* p_dcr2, uint32 eid2,
                                                const Transform3& tr1_2, const Transform3& tr2_0,
                                                const dcr3_element_cache_type& cache1, const dcr3_element_cache_type& cache2,
                                                std::vector< crossing_triangle_pair_type >& vec_cfp,
                                                const Context* p_context = g_pDefaultContext )
{
    switch( p_context->m_DCR2DCR_CFP_Method )
    {
    case Context::eDCR2DCR_CFP_BruteForce: return AddCFP_Intersection_DCR3_E_DCR3_E_BruteForce(p_solid1,vec_sdof1,p_dcr1,eid1,p_solid2,vec_sdof2,p_dcr2,eid2,tr1_2,tr2_0,cache1,cache2,vec_cfp,p_context); break;
#ifdef __ENABLE_SIMD_PRETRANSFORM
    case Context::eDCR2DCR_CFP_BDT: return AddCFP_Intersection_DCR3_E_DCR3_E_BDT_SIMD(p_solid1,vec_sdof1,p_dcr1,eid1,p_solid2,vec_sdof2,p_dcr2,eid2,tr1_2,tr2_0,cache1,cache2,vec_cfp,p_context); break;
#else
    case Context::eDCR2DCR_CFP_BDT: return AddCFP_Intersection_DCR3_E_DCR3_E_BDT(p_solid1,vec_sdof1,p_dcr1,eid1,p_solid2,vec_sdof2,p_dcr2,eid2,tr1_2,tr2_0,cache1,cache2,vec_cfp,p_context); break;
#endif
    default: return 0;
    }
}

/*IMPORTANT: See comments in Contact_X_Y_BVH TestOverlap_BSlabs in 2D

  \note This is used BEFORE the element_caches are available, consider
  declaring them before and using them save 2*Inverse(Mat3x3)
  computations \note cached invBs is in local coords, would still need
  to apply tr1/tr2
*/
bool TestOverlap_BSlabs( const TetSolidShape3* p_solid1, const DCR_TetSolidShape3* p_dcr1, const Transform3& tr1, const TetSolidShape3::sdof_type* vec_sdof1,
                         BVH_TetSolidShape3::entry_index_type eid1, const BVH_TetSolidShape3::bv_type& bv1,
                         const TetSolidShape3* p_solid2, const DCR_TetSolidShape3* p_dcr2, const Transform3& tr2, const TetSolidShape3::sdof_type* vec_sdof2,
                         BVH_TetSolidShape3::entry_index_type eid2, const BVH_TetSolidShape3::bv_type& bv2 )
{
    GEO_ASSERT( eid1 < p_dcr1->m_NumElements && eid2 < p_dcr2->m_NumElements );
    //------- Option (g)
    //--- Elements
    const DCR_TetSolidShape3::Element& ed1( p_dcr1->m_vecE[eid1] );
    const DCR_TetSolidShape3::Element& ed2( p_dcr2->m_vecE[eid2] );
    //--- global nodes
    uint32 vec_nid1[4] = { p_solid1->T_VID(eid1,0), p_solid1->T_VID(eid1,1), p_solid1->T_VID(eid1,2), p_solid1->T_VID(eid1,3) };
    uint32 vec_nid2[4] = { p_solid2->T_VID(eid2,0), p_solid2->T_VID(eid2,1), p_solid2->T_VID(eid2,2), p_solid2->T_VID(eid2,3) };
    Vec3 vec_node_pos1_0[4] = { tr1 * vec_sdof1[vec_nid1[0]], tr1 * vec_sdof1[vec_nid1[1]], tr1 * vec_sdof1[vec_nid1[2]], tr1 * vec_sdof1[vec_nid1[3]] };
    Vec3 vec_node_pos2_0[4] = { tr2 * vec_sdof2[vec_nid2[0]], tr2 * vec_sdof2[vec_nid2[1]], tr2 * vec_sdof2[vec_nid2[2]], tr2 * vec_sdof2[vec_nid2[3]] };
    //--- GDOP1 vs GSlab2_Vertices
    {
        // Barycentric transform Ds1 using node0 as the local origin
        Mat3x3 invDs1( mal::Inverse( mal::GMat3x3_From_Columns( vec_node_pos1_0[1]-vec_node_pos1_0[0],
                                                                vec_node_pos1_0[2]-vec_node_pos1_0[0],
                                                                vec_node_pos1_0[3]-vec_node_pos1_0[0] ) ) );
        // Transform eid2 nodes to Ds1 coords (u,v,t,s)
        Vec3 vec_node_vt[4] = { invDs1*(vec_node_pos2_0[0]-vec_node_pos1_0[0]),
                                invDs1*(vec_node_pos2_0[1]-vec_node_pos1_0[0]),
                                invDs1*(vec_node_pos2_0[2]-vec_node_pos1_0[0]),
                                invDs1*(vec_node_pos2_0[3]-vec_node_pos1_0[0]) };
        Vec4 vec_node_pos2_Ds1[4] = { Vec4( 1-vec_node_vt[0][0]-vec_node_vt[0][1]-vec_node_vt[0][2], vec_node_vt[0][0], vec_node_vt[0][1], vec_node_vt[0][2] ),
                                      Vec4( 1-vec_node_vt[1][0]-vec_node_vt[1][1]-vec_node_vt[1][2], vec_node_vt[1][0], vec_node_vt[1][1], vec_node_vt[1][2] ),
                                      Vec4( 1-vec_node_vt[2][0]-vec_node_vt[2][1]-vec_node_vt[2][2], vec_node_vt[2][0], vec_node_vt[2][1], vec_node_vt[2][2] ),
                                      Vec4( 1-vec_node_vt[3][0]-vec_node_vt[3][1]-vec_node_vt[3][2], vec_node_vt[3][0], vec_node_vt[3][1], vec_node_vt[3][2] ) };
        /* \todo HERE we could test BDOP[i] using the vec_node_pos2_B1
           as an early-out to avoid computing the bslab_pos2_B1, but
           seems overkill */
        // Compute BSlab2 vertices in B1 coords
        int bsi( ed2.m_BDOP_BestSlabIdx );
        Vec4 node_pos_i_times_min( ed2.m_BDOP[bsi].Min() * vec_node_pos2_Ds1[bsi] );
        Real one_minus_min( 1 - ed2.m_BDOP[bsi].Min() );
        Vec4 node_pos_i_times_max( ed2.m_BDOP[bsi].Max() * vec_node_pos2_Ds1[bsi] );
        Real one_minus_max( 1 - ed2.m_BDOP[bsi].Max() );
        Vec4 bslab_pos2_Ds1[6] = { node_pos_i_times_min + one_minus_min * vec_node_pos2_Ds1[(bsi+1)%4],
                                   node_pos_i_times_min + one_minus_min * vec_node_pos2_Ds1[(bsi+2)%4],
                                   node_pos_i_times_min + one_minus_min * vec_node_pos2_Ds1[(bsi+3)%4],
                                   node_pos_i_times_max + one_minus_max * vec_node_pos2_Ds1[(bsi+1)%4],
                                   node_pos_i_times_max + one_minus_max * vec_node_pos2_Ds1[(bsi+2)%4],
                                   node_pos_i_times_max + one_minus_max * vec_node_pos2_Ds1[(bsi+3)%4] };
        // Test BDOP[i], exit if no overlap
        for( int i=0; i<4; i++ )
        {
            Real interval2_Ds1_min( bslab_pos2_Ds1[0][i] );
            for( int j=1; j<6; j++ )
                if( interval2_Ds1_min > bslab_pos2_Ds1[j][i] )
                    interval2_Ds1_min = bslab_pos2_Ds1[j][i];
            if( interval2_Ds1_min > ed1.m_BDOP[i].Max() ) return false;
            Real interval2_Ds1_max( bslab_pos2_Ds1[0][i] );
            for( int j=1; j<6; j++ )
                if( interval2_Ds1_max < bslab_pos2_Ds1[j][i] )
                    interval2_Ds1_max = bslab_pos2_Ds1[j][i];
            if( interval2_Ds1_max < ed1.m_BDOP[i].Min() ) return false;
        }
    }
    //--- GDOP2 vs GSlab1_Vertices
    {
        // Barycentric transform Ds1 using node0 as the local origin
        Mat3x3 invDs2( mal::Inverse( mal::GMat3x3_From_Columns( vec_node_pos2_0[1]-vec_node_pos2_0[0],
                                                                vec_node_pos2_0[2]-vec_node_pos2_0[0],
                                                                vec_node_pos2_0[3]-vec_node_pos2_0[0] ) ) );
        // Transform eid2 nodes to Ds2 coords (u,v,t,s)
        Vec3 vec_node_vt[4] = { invDs2*(vec_node_pos1_0[0]-vec_node_pos2_0[0]),
                                invDs2*(vec_node_pos1_0[1]-vec_node_pos2_0[0]),
                                invDs2*(vec_node_pos1_0[2]-vec_node_pos2_0[0]),
                                invDs2*(vec_node_pos1_0[3]-vec_node_pos2_0[0]) };
        Vec4 vec_node_pos1_Ds2[4] = { Vec4( 1-vec_node_vt[0][0]-vec_node_vt[0][1]-vec_node_vt[0][2], vec_node_vt[0][0], vec_node_vt[0][1], vec_node_vt[0][2] ),
                                      Vec4( 1-vec_node_vt[1][0]-vec_node_vt[1][1]-vec_node_vt[1][2], vec_node_vt[1][0], vec_node_vt[1][1], vec_node_vt[1][2] ),
                                      Vec4( 1-vec_node_vt[2][0]-vec_node_vt[2][1]-vec_node_vt[2][2], vec_node_vt[2][0], vec_node_vt[2][1], vec_node_vt[2][2] ),
                                      Vec4( 1-vec_node_vt[3][0]-vec_node_vt[3][1]-vec_node_vt[3][2], vec_node_vt[3][0], vec_node_vt[3][1], vec_node_vt[3][2] ) };
        /* \todo HERE we could test BDOP[i] using the vec_node_pos2_B1
           as an early-out to avoid computing the bslab_pos2_B1, but
           seems overkill */
        // Compute BSlab2 vertices in B1 coords
        int bsi( ed1.m_BDOP_BestSlabIdx );
        Vec4 node_pos_i_times_min( ed1.m_BDOP[bsi].Min() * vec_node_pos1_Ds2[bsi] );
        Real one_minus_min( 1 - ed1.m_BDOP[bsi].Min() );
        Vec4 node_pos_i_times_max( ed1.m_BDOP[bsi].Max() * vec_node_pos1_Ds2[bsi] );
        Real one_minus_max( 1 - ed1.m_BDOP[bsi].Max() );
        Vec4 bslab_pos1_Ds2[6] = { node_pos_i_times_min + one_minus_min * vec_node_pos1_Ds2[(bsi+1)%4],
                                   node_pos_i_times_min + one_minus_min * vec_node_pos1_Ds2[(bsi+2)%4],
                                   node_pos_i_times_min + one_minus_min * vec_node_pos1_Ds2[(bsi+3)%4],
                                   node_pos_i_times_max + one_minus_max * vec_node_pos1_Ds2[(bsi+1)%4],
                                   node_pos_i_times_max + one_minus_max * vec_node_pos1_Ds2[(bsi+2)%4],
                                   node_pos_i_times_max + one_minus_max * vec_node_pos1_Ds2[(bsi+3)%4] };
        // Test BDOP[i], exit if no overlap
        for( int i=0; i<4; i++ )
        {
            Real interval1_Ds2_min( bslab_pos1_Ds2[0][i] );
            for( int j=1; j<6; j++ )
                if( interval1_Ds2_min > bslab_pos1_Ds2[j][i] )
                    interval1_Ds2_min = bslab_pos1_Ds2[j][i];
            if( interval1_Ds2_min > ed2.m_BDOP[i].Max() ) return false;
            Real interval1_Ds2_max( bslab_pos1_Ds2[0][i] );
            for( int j=1; j<6; j++ )
                if( interval1_Ds2_max < bslab_pos1_Ds2[j][i] )
                    interval1_Ds2_max = bslab_pos1_Ds2[j][i];
            if( interval1_Ds2_max < ed2.m_BDOP[i].Min() ) return false;
        }
    }
    return true;
}

/* DCR3 vs DCR3
   \todo Track BVH status (geomtimestamp) and Refit() only if out-of-date! OR assume up-to-date externally??
   \todo Get BVH geometry decoupled from static BVH topology
*/
bool Contact_DCRTS3_DCRTS3_BVH_TSS_Lazy_DCR( const TetSolidShape3* p_solid1, const Transform3& solid_tr1, const Real* vec_dof1,
                                             const TetSolidShape3* p_solid2, const Transform3& solid_tr2, const Real* vec_dof2,
                                             ContactData3& cd, ContactCache3* p_cc,
                                             const Context* p_context )
{
    bool bLog(p_context->m_DCR2DCR_Log_Enabled);
    bool bViz(p_context->m_DCR2DCR_Viz_Enabled);
    if( bLog ) { GEO_LOG("******************************** Contact_DCRTS3_DCRTS3_BVH_TSS_Lazy_DCR() ********************************"); }
    cd.Begin();

    //\todo THIS IS dangerous... Should have a generic SRV type with safe casts...
    const Vec3* vec_sdof1( reinterpret_cast<const Vec3*>(vec_dof1) );
    const Vec3* vec_sdof2( reinterpret_cast<const Vec3*>(vec_dof2) );

    const DCR_TetSolidShape3* pDCR1( p_solid1->GetDCR() );
    const DCR_TetSolidShape3* pDCR2( p_solid2->GetDCR() );
    if( 0 == pDCR1 || 0 == pDCR2 )
    {
        GEO_LOG_WARNING("Contact_DCRTS3_DCRTS3_BVH_MSS_Lazy_DCR() no DCR1 && DCR2, ignoring...");
        return false; //\todo Fallback to a different method?
    }

    //TEMP: By now, we'll just add the BVH if any is missing
    BVH_TetSolidShape3* pBVH1( p_solid1->GetBVH() );
    BVH_TetSolidShape3* pBVH2( p_solid2->GetBVH() );
    if( 0 == pBVH1 || 0 == pBVH2 )
    {
        if( 0 == pBVH1 )
        {
            pBVH1 = new BVH_TetSolidShape3;
            uint32 num_entries( pDCR1->m_NumElements );
            pBVH1->Rebuild_TopDown( num_entries,
                                    boost::bind<void>( &GEBV_TetSolidShape3_E<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
                                                       p_solid1, solid_tr1, vec_sdof1, //Transform3::Identity(), default_sdof1,
                                                       _1, _2 ) );
            GEO_LOG_WARNING("Contact_DCRTS3_DCRTS3_BVH_MSS_Lazy_DCR() no BVH1, creating on the fly, UGLY...");
            const_cast<TetSolidShape3*>(p_solid1)->SetBakedBVH_StrictlyNonshared_UglyHack( pBVH1 ); //ugly const cast...
        }
        if( 0 == pBVH2 )
        {
            pBVH2 = new BVH_TetSolidShape3;
            uint32 num_entries( pDCR2->m_NumElements );
            pBVH2->Rebuild_TopDown( num_entries,
                                    boost::bind<void>( &GEBV_TetSolidShape3_E<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
                                                       p_solid2, solid_tr2, vec_sdof2, //Transform3::Identity(), default_sdof2,
                                                       _1, _2 ) );
            GEO_LOG_WARNING("Contact_DCRTS3_DCRTS3_BVH_MSS_Lazy_DCR() no BVH2, creating on the fly, UGLY...");
            const_cast<TetSolidShape3*>(p_solid2)->SetBakedBVH_StrictlyNonshared_UglyHack( pBVH2 ); //ugly const cast...
        }
    }

    //--- Refit BVH \todo Do lazily
    /*\todo Instead of re-transforming nodes for each element they're
      involved in, consider pre-transforming them all in a single pass
      into a frame-allocated cache and access them by index in the
      GEBV functor by index.
      Ex:
        util::ScopedAllocator scoped_allocator( p_context->m_ScratchPad, "Name..." );
        Vec3* vec_node_pos1_0 = scoped_allocator.NewArrayPOD<Vec3>(pDCR1->m_NumVertices);
        mal::GTransform_Point_Array( vec_sdof1, tr1, vec_node_pos1_0 ); //\todo easily simd-able
      Then, invoke Refit() with pre-transformed nodes
        pBVH1->Refit( boost::bind<void>( GEBV_BV_E_transformed_nodes, *p_solid1, pDCR1, vec_node_pos1_0 ), _1, _2) );
    */
    GEO_NP_BEGIN_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_RefitBVH_ms );
    switch( p_context->m_BVH_Method )
    {
    case Context::eBVHM_BV_E:
        pBVH1->Refit( boost::bind<void>( &GEBV_TetSolidShape3_E<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
                                         p_solid1, solid_tr1, vec_sdof1, _1, _2) );
        pBVH2->Refit( boost::bind<void>( &GEBV_TetSolidShape3_E<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
                                         p_solid2, solid_tr2, vec_sdof2, _1, _2) );
        break;
        //\todo Other methods not yet available in 3D
    // case Context::eBVHM_BV_BSlab:
    //     pBVH1->Refit( boost::bind<void>( &GEBV_TetSolidShape3_BSlab<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
    //                                      p_solid1, pDCR1, solid_tr1, vec_sdof1, _1, _2) );
    //     pBVH2->Refit( boost::bind<void>( &GEBV_TetSolidShape3_BSlab<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
    //                                      p_solid2, pDCR2, solid_tr2, vec_sdof2, _1, _2) );
    //     break;
    case Context::eBVHM_BV_BDOP:
        pBVH1->Refit( boost::bind<void>( &GEBV_TetSolidShape3_BDOP<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
                                         p_solid1, pDCR1, solid_tr1, vec_sdof1, _1, _2) );
        pBVH2->Refit( boost::bind<void>( &GEBV_TetSolidShape3_BDOP<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
                                         p_solid2, pDCR2, solid_tr2, vec_sdof2, _1, _2) );
        break;
    case Context::eBVHM_BV_NoRefit: break;
    default: break;
    }
    GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_RefitBVH_ms );

    // Test overlap \todo BUILD INSIDE SCOPED ALLOCATOR!!
    std::vector< std::pair<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::entry_index_type> > vecOverlaps;

    GEO_NP_BEGIN_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_OverlapBVH_ms );
    bool bOverlap(false);
    switch( p_context->m_TestBVH_Method )
    {
    case Context::eTBVHM_BV_Vs_BV: bOverlap = pBVH1->Test( *pBVH2, vecOverlaps ); break;
    case Context::eTBVHM_BDOP_Vs_GSlabVertices:
        bOverlap = pBVH1->Test( *pBVH2,
                                [ p_solid1, p_solid2, pDCR1, pDCR2, solid_tr1, solid_tr2, vec_sdof1, vec_sdof2 ]
                                ( BVH_TetSolidShape3::entry_index_type eid1, BVH_TetSolidShape3::entry_index_type eid2,
                                  const BVH_TetSolidShape3::bv_type& bv1, const BVH_TetSolidShape3::bv_type& bv2 )
                                {
                                    GEO_ASSERT( eid1 < pDCR1->m_NumElements && eid2 < pDCR2->m_NumElements );
                                    return TestOverlap(bv1,bv2) //early-out, could be merged into TestBSlabs to reuse partial computations
                                    //&& TestTetrahedrons(eid1,eid2) //early-out, could be merged into TestBSlabs to reuse partial computations
                                    && TestOverlap_BSlabs( p_solid1, pDCR1, solid_tr1, vec_sdof1, eid1, bv1,
                                                           p_solid2, pDCR2, solid_tr2, vec_sdof2, eid2, bv2 );
                                },
                                vecOverlaps );
        break;
    default: break;
    }
    GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_OverlapBVH_ms );
    if( !bOverlap ) return false;

    dcr3_element_cache_type element_cache1( p_solid1, pDCR1, solid_tr1, vec_sdof1 );
    dcr3_element_cache_type element_cache2( p_solid2, pDCR2, solid_tr2, vec_sdof2 );

    Transform3 tr1_2( mal::Inverse(solid_tr2) * solid_tr1 );
    std::vector< crossing_triangle_pair_type > vecCTP;
    if( bLog ) { GEO_LOG("Computing CTP..."); }

    GEO_NP_BEGIN_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_FindITP_ms );
    for( auto overlap : vecOverlaps )
    {
        /*IMPORTANT!! If ANY object BVH includes internal elements
          NOT IN the DCR, this can happen, and we MUST AVOID
          accessing the inexistent DCR.E
          \todo CONSIDER detecting internal/internal overlaps and
          using/reporting them IF NO boundary overlap is detected.
        */
        if( overlap.first >= pDCR1->m_NumElements || overlap.second >= pDCR2->m_NumElements ) continue;

        AddCFP_Intersection_DCR3_E_DCR3_E( p_solid1, vec_sdof1, pDCR1, overlap.first,
                                           p_solid2, vec_sdof2, pDCR2, overlap.second,
                                           tr1_2, solid_tr2,
                                           element_cache1, element_cache2,
                                           vecCTP, p_context );
    }
    GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_FindITP_ms );
    GEO_NP_STAT_ACC( mig2015.m_Num_ITP, (uint32)vecCTP.size() );
    if( vecCTP.empty() ) return false;
    if( bLog ) { GEO_LOG("#CTP = %u",(uint32)vecCTP.size()); }
    if( bLog )
    {
        uint32 stats_total_tr_e_1(0), stats_acc_tr_e_1(0), stats_max_tr_e_1(0);
        uint32 stats_total_tr_e_2(0), stats_acc_tr_e_2(0), stats_max_tr_e_2(0);
        for( uint32 it_e=0; it_e<pDCR1->m_NumElements; it_e++ )
        {
            stats_acc_tr_e_1 += element_cache1[it_e].m_Stats_NumTransform;
            if( element_cache1[it_e].m_Stats_NumTransform > 0 )
                stats_total_tr_e_1++;
            if( stats_max_tr_e_1 < element_cache1[it_e].m_Stats_NumTransform )
                stats_max_tr_e_1 = element_cache1[it_e].m_Stats_NumTransform;
        }
        for( uint32 it_e=0; it_e<pDCR2->m_NumElements; it_e++ )
        {
            stats_acc_tr_e_2 += element_cache2[it_e].m_Stats_NumTransform;
            if( element_cache2[it_e].m_Stats_NumTransform > 0 )
                stats_total_tr_e_2++;
            if( stats_max_tr_e_2 < element_cache2[it_e].m_Stats_NumTransform )
                stats_max_tr_e_2 = element_cache2[it_e].m_Stats_NumTransform;
        }
        GEO_LOG("Transformed\n\tDCR1.E = total %u, acc %u, max %u, avg %f\n\tDCR2.E = total %u, acc %u, max %u, avg %f",
                stats_total_tr_e_1, stats_acc_tr_e_1, stats_max_tr_e_1, float(stats_acc_tr_e_1)/stats_total_tr_e_1,
                stats_total_tr_e_2, stats_acc_tr_e_2, stats_max_tr_e_2, float(stats_acc_tr_e_2)/stats_total_tr_e_2 );
    }

    //TEMP DRAW CTP as CP
    // for( auto cfp : vecCTP )
    // {
    //     cd.AddCP( cfp.m_ITTR.m_Point0, cfp.m_ITTR.m_Point1, mal::SafeNormalized(cfp.m_ITTR.m_Point1-cfp.m_ITTR.m_Point0), mal::Norm(cfp.m_ITTR.m_Point1-cfp.m_ITTR.m_Point0) );
    //     cd.AddPOF( PointOnFeature(), PointOnFeature() );
    // }

    /* Compute IntersectionCurves (IC)
       - FloodIB requires:
       - A seed CTP with unambiguous inwards/outwards flood direction
       - Any CTP in an IC is equally valid
       - A way to avoid "splilling" the flood outside the IB
       - Detect crossing of *any* CTP in the IC, regardless of IC ordering!
       - A way to avoid re-flooding wet areas
       - IC.CTP.IsOpen
       - All this needs to be done on both sides DCR1,DCR2
       => Thus, we do NOT NEED any special ordering for the CTP in an IC
       - ComputeIC only needs to gather all spatially pairwise-adjacent CTP into an unordered IC,
       - We'll use an MF-Set alg that merges CTP into IC as adjacencies are found

       => This will ACCEPT non-closed IC, which may fail to
       contain FloodIB... consider checking/enforcing closedness
       WITHOUT requiring sequantiallity (each CTP must have both
       ends "connected", but no strict order is required)
    */
    if( !p_context->m_DCR2DCR_IC_Enabled ) return false; //\note Enable/Disable to benchmark specific phases

    GEO_NP_BEGIN_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_MergeIC_ms );

    if( bLog ) { GEO_LOG("Creating IC..."); }
    std::vector< dcr3_intersection_curve_type > vecIC;
    // Compute AABB intersection between pDCR1 pDCR2
    bv::AABB3 aabb1( bv::GAABB_From( pBVH1->GetRoot().m_Geometry.m_BV ) );
    bv::AABB3 aabb2( bv::GAABB_From( pBVH2->GetRoot().m_Geometry.m_BV ) );
    bv::AABB3 aabb_intersection( bv::GAABB_Intersect(aabb1,aabb2) );
    /*
      GEO_LOG( "AABB1 = (%f,%f,%F) ; (%f,%f,%F)",
      aabb1.GetMin().x(), aabb1.GetMin().y(), aabb1.GetMin().z(),
      aabb1.GetMax().x(), aabb1.GetMax().y(), aabb1.GetMax().z() );
      GEO_LOG( "AABB2 = (%f,%f,%F) ; (%f,%f,%F)",
      aabb2.GetMin().x(), aabb2.GetMin().y(), aabb2.GetMin().z(),
      aabb2.GetMax().x(), aabb2.GetMax().y(), aabb2.GetMax().z() );
      GEO_LOG( "AABB vs AABB = (%f,%f,%F) ; (%f,%f,%F)",
      aabb_intersection.GetMin().x(), aabb_intersection.GetMin().y(), aabb_intersection.GetMin().z(),
      aabb_intersection.GetMax().x(), aabb_intersection.GetMax().y(), aabb_intersection.GetMax().z() );
    */
    // 1st SH pass: Pre-add both CTP endpoints
    const int cMaxDimension(32);
    uint32 dimensions[3] = {cMaxDimension,cMaxDimension,cMaxDimension};
    // Compute SH resolution for the largest extent
    Vec3 sizes( Real(2)*aabb_intersection.GetHalfSizes() );
    Real resolution = mal::Max(sizes) / cMaxDimension;
    for( int i=0; i<3; i++ )
        dimensions[i] = mal::Clamp<uint32>( mal::IntPart<uint32>( mal::Ceil( sizes[i] / resolution ) ), 1, cMaxDimension );

    //\todo THIS SHOULD BE scope-allocated with fixed resolution (cMaxDimension^3)
    GSimpleSpatialHash< uint32, crossing_triangle_pair_type::endpoint_id_type > ssh( 2*vecCTP.size(),
                                                                                     dimensions,
                                                                                     aabb_intersection.GetMin(),
                                                                                     aabb_intersection.GetMax() );
    // First pass: classify CTP midpoints
    ssh.BeginClassify();
    {
        for( unsigned int it_ctp=0; it_ctp<vecCTP.size(); it_ctp++ )
        {
            const crossing_triangle_pair_type& ctp( vecCTP[it_ctp] );
            ssh.Classify( crossing_triangle_pair_type::endpoint_id_type(it_ctp,0), ctp.m_ITTR.m_Point0 );
            ssh.Classify( crossing_triangle_pair_type::endpoint_id_type(it_ctp,1), ctp.m_ITTR.m_Point1 );
        }
    }
    ssh.EndClassify();

    // 2st SH pass: Query both CTP endpoints, merge with actual overlapping CTP
    // Reset MF-Set parents
    for( unsigned int it_ctp=0; it_ctp<vecCTP.size(); it_ctp++ ) vecCTP[it_ctp].m_ParentId = it_ctp;
    //\todo Prec 1e-6 WORKS BETTER but I think it's due to bAdjacent being TRUE when it shouldn't, which causes already connected matches due to shortcuts...
    const Real cEpsilonComputeIC(mal::Exp10(p_context->m_DCR2DCR_CFP_EpsilonIC_Log10)); //\todo Consider storing ITTR geometry in in solid2 ref sys instead of globally to improve prec far from the origin
    const Real cEpsilonComputeIC_Sq(mal::Sq(cEpsilonComputeIC));
    uint32 num_tests_sh(0);
    uint32 num_short_segments(0);
    uint32 num_incompatible_triangle_pair_matches(0);
    uint32 num_already_connected_matches(0);
    uint32 num_non_adjacent_matches(0);
    ssh.BeginTestAdd();
    {
        std::vector<crossing_triangle_pair_type::endpoint_id_type> vec_potential_overlaps;
        for( unsigned int it_ctp=0; it_ctp<vecCTP.size(); it_ctp++ )
        {
            crossing_triangle_pair_type& ctp( vecCTP[it_ctp] );
            for( int it_epictp=0; it_epictp<2; it_epictp++ )
            {
                Vec3 point_in_ctp_segment( it_epictp==0 ? ctp.m_ITTR.m_Point0 : ctp.m_ITTR.m_Point1 ); //\todo ctp.m_Segment[it_point_in_ctp_segment],
                crossing_triangle_pair_type::endpoint_id_type epid( it_ctp, it_epictp );
                if( ssh.TestAdd( epid, point_in_ctp_segment, cEpsilonComputeIC, vec_potential_overlaps ) )
                {
                    // Search for adjacent overlap with best matching (there should only be 1 adjacent overlap, but numeric fuckups happen (AlreadyConnectedMatches>0)
                    crossing_triangle_pair_type::endpoint_id_type best_oepid; //invalid
                    Real best_dist_sq(mal::Infinity<Real>());
                    for( auto oepid : vec_potential_overlaps )
                    {
                        unsigned int it_octp(oepid.CTPID());
                        crossing_triangle_pair_type& octp( vecCTP[it_octp] );
                        Vec3 point_in_octp_segment( oepid.EPICTP() == 0 ? octp.m_ITTR.m_Point0 : octp.m_ITTR.m_Point1 ); //\todo octp.m_Segment[it_point_in_octp_segment],
                        Real dist_sq( mal::NormSq(point_in_ctp_segment-point_in_octp_segment) );
                        if( dist_sq < cEpsilonComputeIC_Sq )
                        {
                            if( it_ctp != it_octp ) //TEMP DEBUG: epsilon-length ITTR segments yield this case, both VF/FV and EE
                            {
                                /* Check that matches are topologically adjacent through
                                   a shared edge and that octp is not already connected */
                                if( !CTP_AreAdjacent( ctp, epid.EPICTP(),
                                                      octp, oepid.EPICTP(),
                                                      pDCR1, pDCR2 ) )
                                    num_non_adjacent_matches++;
                                //TEMP: else if( octp.m_NumConnected[it_point_in_octp_segment] > 0 ) //This forces splitting self-crossing IC
                                else if( octp.m_NeighbourEPID[oepid.EPICTP()].IsValid() ) //This forces splitting self-crossing IC
                                    num_already_connected_matches++;
                                else if( dist_sq < best_dist_sq ) //includes first match
                                {
                                    best_oepid = oepid;
                                    best_dist_sq = dist_sq;
                                }
                            }
                            else
                                num_short_segments++;
                        }
                    }
                    // Merge with best_oepid if any
                    if( best_oepid.IsValid() )
                    {
                        unsigned int it_octp(best_oepid.CTPID());
                        crossing_triangle_pair_type& best_octp( vecCTP[it_octp] );
                        // path-compress IMPORTANT: Required because octp can already have a valid root, which MUST be merged into this ctp set
                        uint32 root_octpid( CTP_ComputeRootId(vecCTP,best_octp.m_ParentId) );
                        // merge previous root
                        vecCTP[root_octpid].m_ParentId = ctp.m_ParentId;
                        // count connections to detect open endpoints
                        ctp.m_NumConnected[epid.EPICTP()]++;
                        best_octp.m_NumConnected[best_oepid.EPICTP()]++;
                        // connect neighbours
                        ctp.m_NeighbourEPID[epid.EPICTP()] = best_oepid;
                        best_octp.m_NeighbourEPID[best_oepid.EPICTP()] = epid;
                    }
                    num_tests_sh += vec_potential_overlaps.size();
                }
            }
        }
    }
    ssh.EndTestAdd();
    // GEO_LOG( "Contact_DCRTS3_DCRTS3_BVH_MSS_Lazy_DCR() ComputeIC SH tests = %u << %u", num_tests_sh, uint32(vecCTP.size() * (vecCTP.size()-1)) );
    if( num_already_connected_matches > 0 )
    {
        GEO_LOG_WARNING( "Contact_DCRTS3_DCRTS3_BVH_MSS_Lazy_DCR() ComputeIC: #Short = %u (<%e), #ITP = %u ,#NAM = %u, #ACM = %u",
                         num_short_segments, cEpsilonComputeIC, num_incompatible_triangle_pair_matches, num_non_adjacent_matches, num_already_connected_matches );
    }

    // Create an IC for each unique CTP.m_ParentIC
    uint32 num_unconnected_endpoints(0);
    uint32 num_multiconnected_endpoints(0);
    for( unsigned int it_ctp=0; it_ctp<vecCTP.size(); it_ctp++ )
    {
        // Find root CTP
        uint32 root_ctpid( CTP_ComputeRootId(vecCTP,it_ctp) );
        crossing_triangle_pair_type& root_ctp( vecCTP[root_ctpid] );
        // Create root IC if none
        if( root_ctp.m_CurveId == 0xFFFFFFFF )
        {
            root_ctp.m_CurveId = vecIC.size();
            vecIC.push_back( dcr3_intersection_curve_type() );
            // GEO_LOG("Creating root %u curve %u",root_ctpid,root_ctp.m_CurveId);
        }
        // Add to root IC
        crossing_triangle_pair_type& ctp( vecCTP[it_ctp] );
        vecIC[ root_ctp.m_CurveId ].m_vecCTPID.push_back( it_ctp );
        ctp.m_CurveId = root_ctp.m_CurveId;
        // GEO_LOG("Adding ctp %u to root %u, curve %u",it_ctp,root_ctpid,root_ctp.m_CurveId);

        if( ctp.m_NumConnected[0] == 0 ) num_unconnected_endpoints++;
        else if( ctp.m_NumConnected[0] > 1 ) num_multiconnected_endpoints++;
        if( ctp.m_NumConnected[1] == 0 ) num_unconnected_endpoints++;
        else if( ctp.m_NumConnected[1] > 1 ) num_multiconnected_endpoints++;
    }
    if( bLog ) { GEO_LOG( "#IC = %u", (uint32)vecIC.size() ); }
    if( num_unconnected_endpoints > 0 || num_multiconnected_endpoints > 0 )
    {
        GEO_LOG_WARNING( "Contact_DCRTS3_DCRTS3_BVH_MSS_Lazy_DCR() ComputeIC #IC = %u, #Unconnected = %u #Multiconnected = %u",
                         (uint32)vecIC.size(), num_unconnected_endpoints, num_multiconnected_endpoints );
        //\todo If UNCONNECTED>0 but there's no other problems (such
        //as ACM>0), CONSIDER retrying merge using a larger epsilon,
        //as it's probably a precision problem (ex: 1 curve split into
        //2 open but consistent curves)
    }
    /* Enforce IC CCW orientation on side[0] (CW on side[1]), so that
       the clipped-triangles on side[0] have a CTP-contributed edge
       CCW as the unclipped T-edges
    */
    for( auto& ic : vecIC )
    {
        for( auto ctpid : ic.m_vecCTPID )
        {
            crossing_triangle_pair_type& ctp( vecCTP[ctpid] );
            Vec3 n0( mal::Cross( ctp.m_vecPos_0[0][1] - ctp.m_vecPos_0[0][0], ctp.m_vecPos_0[0][2] - ctp.m_vecPos_0[0][0] ) );
            Vec3 n1( mal::Cross( ctp.m_vecPos_0[1][1] - ctp.m_vecPos_0[1][0], ctp.m_vecPos_0[1][2] - ctp.m_vecPos_0[1][0] ) );
            if( mal::Dot(ctp.m_ITTR.m_Point1-ctp.m_ITTR.m_Point0, mal::Cross(n0,n1) ) < 0 )
            {
                // Swap 0/1 stuff
                std::swap(ctp.m_ITTR.m_Point0,ctp.m_ITTR.m_Point1); //TEMP: This BREAKS ITTR internal conventions...
                std::swap(ctp.m_NumConnected[0],ctp.m_NumConnected[1]);
                std::swap(ctp.m_NeighbourEPID[0],ctp.m_NeighbourEPID[1]);
                // Reconnect endpoints on neighbour CTP to the swapped endpoints
                for( int it_epictp=0; it_epictp<2; it_epictp++ )
                    if( ctp.m_NeighbourEPID[it_epictp].IsValid() )
                    {
                        crossing_triangle_pair_type& nctp( vecCTP[ctp.m_NeighbourEPID[it_epictp].CTPID()] );
                        for( int it_nepictp=0; it_nepictp<2; it_nepictp++ )
                            if( nctp.m_NeighbourEPID[it_nepictp].IsValid()
                                && nctp.m_NeighbourEPID[it_nepictp].CTPID() == ctpid )
                                nctp.m_NeighbourEPID[it_nepictp] = crossing_triangle_pair_type::endpoint_id_type(ctpid,it_epictp);
                    }
            }
        }
    }
    //\todo I'm almost SURE that the integrated normal can be computed
    //using exclusively the CCW edges (see DGP tutorial/external
    //calculus thing, they mention it)
    GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_MergeIC_ms );
    GEO_NP_STAT_ACC( mig2015.m_Num_IC, (uint32)vecIC.size() );

    // DEBUG IC
    if( bViz )
    {
        for( const auto& ic : vecIC )
        {
            cd.m_VD.m_vecDCR2DCR_vecIC.push_back( std::vector< std::pair<Vec3,Vec3> >() );
            /*
              for( auto ctpid : ic.m_vecCTPID )
              cd.m_VD.m_vecDCR2DCR_vecIC.back().push_back( std::make_pair(vecCTP[ctpid].m_ITTR.m_Point0,
              vecCTP[ctpid].m_ITTR.m_Point1) );
            */
            //topology iteration to test flipped CTP, seems OK
            uint32 ctpid = ic.m_vecCTPID[0];
            uint32 num_iter(0);
            bool bOpen(false);
            do
            {
                cd.m_VD.m_vecDCR2DCR_vecIC.back().push_back( std::make_pair(vecCTP[ctpid].m_ITTR.m_Point0,
                                                                            vecCTP[ctpid].m_ITTR.m_Point1) );
                bOpen = !vecCTP[ctpid].m_NeighbourEPID[1].IsValid();
                ctpid = vecCTP[ctpid].m_NeighbourEPID[1].CTPID();
                num_iter++;
            } while( !bOpen && ctpid != ic.m_vecCTPID[0]
                     && num_iter <= ic.m_vecCTPID.size() );
            if( num_iter > ic.m_vecCTPID.size() )
            {
                GEO_LOG_ERROR("INFINITE LOOP DETECTED and aborted... we observed alternating between CTPID 192 and 202 here, while starting at CTPID 0");
            }
        }
        for( const auto& ctp : vecCTP )
        {
            if( ctp.m_NumConnected[0] == 0 ) cd.m_VD.m_vecDCR2DCR_vecUnconnectedCTP.push_back( ctp.m_ITTR.m_Point0 );
            if( ctp.m_NumConnected[1] == 0 ) cd.m_VD.m_vecDCR2DCR_vecUnconnectedCTP.push_back( ctp.m_ITTR.m_Point1 );
        }
    }

    // if( p_context->m_DCR2DCR_IB_Method == Context::eDCR2DCR_IB_None ) return false;

    // 2) Flood IB \todo COULD BE DONE IN PARALLEL if vecCTP was duplicated and merged afterwards
    if( bLog ) { GEO_LOG("Creating IB..."); }
    std::vector<dcr3_intersection_boundary_type> vecIB[2];
    GEO_NP_BEGIN_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_FindCP_ms );
    if( !FloodIB( vecIB[0],
                  vecCTP, vecIC,
                  p_solid1, pDCR1, solid_tr1, vec_sdof1,
                  0, //side 0 => IB1
                  true, //inwards
                  element_cache1,
                  // debug
                  cd,
                  p_context ) ) return false;
    if( !FloodIB( vecIB[1],
                  vecCTP, vecIC,
                  p_solid2, pDCR2, solid_tr2, vec_sdof2,
                  1, //side 1 => IB2
                  true, //inwards
                  element_cache2,
                  // debug
                  cd,
                  p_context ) ) return false;
    GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_FindCP_ms );
    GEO_NP_STAT_ACC( mig2015.m_Num_IB, (uint32)(vecIB[0].size() + vecIB[1].size()) );
    if( bLog ) { GEO_LOG( "#IB1 = %u, #IB2 = %u", (uint32)vecIB[0].size(), (uint32)vecIB[1].size() ); }

    // 5) Generate CP from EC data
    if( p_context->m_DCR2DCR_RC_Method == Context::eDCR2DCR_RC_None ) return false;

    Real cDotThreshold( mal::Cos(mal::Deg2Rad(90.0f-p_context->m_DCR2DCR_RC_MinAngle_Deg)) ); //TEMP: Hack for MIG2015

    Real largest_extent_x_Side[2] = { bv::ComputeLargestExtent(pBVH1->GetRoot().m_Geometry.m_BV), bv::ComputeLargestExtent(pBVH2->GetRoot().m_Geometry.m_BV) };
    const TetSolidShape3* p_solid_x_Side[2] = { p_solid1, p_solid2 };
    const DCR_TetSolidShape3* pDCR_x_Side[2] = { pDCR1, pDCR2 };
    const BVH_TetSolidShape3* pBVH_x_Side[2] = { pBVH1, pBVH2 };
    Transform3 solid_tr_x_Side[2] = { solid_tr1, solid_tr2 };
    Transform3 inv_solid_tr_x_Side[2] = { mal::Inverse(solid_tr1), mal::Inverse(solid_tr2) };
    const dcr3_element_cache_type* p_element_cache_x_Side[2] = { &element_cache1, &element_cache2 };
    for( int side=0; side<2; side++ )
    {
        int oside(1-side);
        for( const auto& ib : vecIB[side] )
        {
            if( ib.m_NumTriangles > 0 )
            {
                Vec3 N_0( mal::SafeNormalized( ib.m_VectorArea_0 ) ); //\todo CAN BE 0, in this case we should probably ignore the IB or, better, subdivide it...
                for( const auto& patch : ib.m_vecP )
                {
                    if( patch.m_NumTriangles > 0 )
                    {
                        uint32 eid( pDCR_x_Side[side]->m_vecP[patch.m_PID].m_EID );
                        Vec3 p_0( patch.m_AvgPos_0 );
                        Vec3 n_0( mal::SafeNormalized( patch.m_VectorArea_0 ) ); //\todo CAN BE 0, in this case we should probably ignore the IB or, better, subdivide it...
                        Vec4 p_b( (*p_element_cache_x_Side[side])[eid].invBs() * mal::Concat(1, inv_solid_tr_x_Side[side] * p_0) ); //\todo THIS MAY BE WRONG...
                        GEO_LOG_ERROR_IF( p_b[0] < 0 || p_b[0] > 1
                                          || p_b[1] < 0 || p_b[1] > 1
                                          || p_b[2] < 0 || p_b[2] > 1
                                          || p_b[3] < 0 || p_b[3] > 1,
                                          "DCR[%d].E[%u] p_b = (%f,%f,%f,%f)",
                                          side, eid, p_b[0], p_b[1], p_b[2], p_b[3] );

                        switch( p_context->m_DCR2DCR_IM_Method )
                        {
                        case Context::eDCR2DCR_IM_Basic: break;
                        case Context::eDCR2DCR_IM_Projection: break;
                        case Context::eDCR2DCR_IM_Raycast: //c) Raycast: Raycast along EC normal, find other point, generate 2 CP with full info (barycentric coords wrt DCR.E, etc...)
                            {
                                // Select raycast direction
                                Vec3 ray_dir_0(0);
                                switch( p_context->m_DCR2DCR_RC_Method )
                                {
                                case Context::eDCR2DCR_RC_Global: //I
                                    /* \todo THIS WILL ONLY WORK for 1 IB1 <==> 1
                                     * IB2, as it assumes N1 = -N2... for complex
                                     * regions, a SINGLE IB1.N should be selected
                                     * and negated for its >1 correspondent IB2
                                     * (or the other way around!) */
                                    ray_dir_0 = -N_0;
                                    break;
                                case Context::eDCR2DCR_RC_Direct: //II
                                    ray_dir_0 = -n_0;
                                    break;
                                default: //\todo Inverse and Average not handled here
                                    return false;
                                }

                                // Perform raycast against IB/DCR, ignoring if non-inwards wrt normal
                                bool bHit( false );
                                RayHit3 rh;
                                if( p_context->m_DCR2DCR_RC_UseOnlyIB ) //1) IB Raycast
                                    return false;
                                else //2) DCR Raycast
                                {
                                    /* IMPORTANT: I've observed
                                     * raycast errors in
                                     * near-coplanarity cases, where
                                     * the ray cast from one side
                                     * fails to hit its almost
                                     * coincident counterpart triangle
                                     * on the oside and continues past
                                     * it, potentially hitting a
                                     * differen triangle much
                                     * farther. This induces a
                                     * spurious deep contact point
                                     * that causes a large contact
                                     * force/reaction discontinuity
                                     */
                                    bHit = mal::Dot(ray_dir_0,n_0) <= -cDotThreshold //\todo TESTING: Threshold allows ignoring patches near-coplanar with the global normal, even if backfacing
                                           && RaycastDCR( p_solid_x_Side[oside], pDCR_x_Side[oside], pBVH_x_Side[oside],
                                                          solid_tr_x_Side[oside], inv_solid_tr_x_Side[oside],
                                                          *p_element_cache_x_Side[oside],
                                                          //p_0 - p_context->m_DCR2DCR_RC_Thickness*ray_dir_0, ray_dir_0,
                                                          p_0, ray_dir_0,
                                                          largest_extent_x_Side[oside], rh, p_context );
                                }
                                //TEMP: Avoid, avoid other-hit farther than closest self-hit, for
                                if( bHit
                                    && rh.m_Interval.Min() > p_context->m_DCR2DCR_RC_Thickness )
                                {
                                    GEO_LOG( "DCR[%d].E[%u] thick-hit at %f > %f", side, eid, rh.m_Interval.Min(), p_context->m_DCR2DCR_RC_Thickness );
                                    RayHit3 self_rh;
                                    if( RaycastDCR( p_solid_x_Side[side], pDCR_x_Side[side], pBVH_x_Side[side],
                                                    solid_tr_x_Side[side], inv_solid_tr_x_Side[side],
                                                    *p_element_cache_x_Side[side],
                                                    p_0, ray_dir_0,
                                                    largest_extent_x_Side[side], self_rh, p_context )
                                        && self_rh.m_Interval.Min() < rh.m_Interval.Min() )
                                    {
                                        GEO_LOG( "DCR[%d].E[%u] found self-hit %f closer than thick-hit %f, ignoring CP ", side, eid, self_rh.m_Interval.Min(), rh.m_Interval.Min() );
                                        bHit = false;
                                    }
                                }
                                if( bHit )
                                {
                                    GEO_LOG_ERROR_IF( rh.m_Extra_BarycentricCoords[0] < 0 || rh.m_Extra_BarycentricCoords[0] > 1
                                                      || rh.m_Extra_BarycentricCoords[1] < 0 || rh.m_Extra_BarycentricCoords[1] > 1
                                                      || rh.m_Extra_BarycentricCoords[2] < 0 || rh.m_Extra_BarycentricCoords[2] > 1
                                                      || rh.m_Extra_BarycentricCoords[3] < 0 || rh.m_Extra_BarycentricCoords[3] > 1,
                                                      "DCR[%d].E[%u] RH_b = (%f,%f,%f,%f)",
                                                      side, eid, rh.m_Extra_BarycentricCoords[0], rh.m_Extra_BarycentricCoords[1], rh.m_Extra_BarycentricCoords[2], rh.m_Extra_BarycentricCoords[3] );

                                    //BUG: far hits are sometimes due
                                    //to coplanarity, both ray_pos
                                    //slightly behind triangle plane,
                                    //or ray_dir being perpendicular
                                    //to triangle normal, we HACK a
                                    //solution by throwing a backwards
                                    //ray and, if hit within small
                                    //threshold, accept it instead
                                    //IMPORTANT: Backwards ray
                                    //effectivelt inverts the normal
                                    //and thus can cause
                                    //sticking... we'll see
                                    if( rh.m_Interval.Min() > p_context->m_DCR2DCR_RC_FarHitThreshold )
                                    {
                                        GEO_LOG( "DCR[%d].E[%u] far-hit at %f > %f", side, eid, rh.m_Interval.Min(), p_context->m_DCR2DCR_RC_FarHitThreshold );
                                        RayHit3 backwards_rh;
                                        bool bHit2 = RaycastDCR( p_solid_x_Side[oside], pDCR_x_Side[oside], pBVH_x_Side[oside],
                                                                 solid_tr_x_Side[oside], inv_solid_tr_x_Side[oside],
                                                                 *p_element_cache_x_Side[oside],
                                                                 p_0, -ray_dir_0,
                                                                 largest_extent_x_Side[oside], backwards_rh, p_context );
                                        if( bHit2 && backwards_rh.m_Interval.Min() < p_context->m_DCR2DCR_RC_FarHitThreshold )
                                        {
                                            GEO_LOG( "DCR[%d].E[%u] found closer backwards hit at %f", side, eid, backwards_rh.m_Interval.Min() );
                                            rh = backwards_rh;
                                        }
                                    }

                                    if( side == 0 )
                                    {
                                        cd.AddCP( p_0,
                                                  rh.m_Point,
                                                  mal::SafeNormalized(rh.m_Point-p_0),
                                                  mal::Norm(rh.m_Point-p_0),
                                                  patch.m_Area );
                                        // mal::Sqrt(patch.m_AreaSq) );
                                        cd.AddPOF( PointOnFeature( feature_id(eFT_Tetrahedron,eid), p_b ),
                                                   PointOnFeature( feature_id(eFT_Tetrahedron,rh.m_FeatureId.AsTetrahedron()), rh.m_Extra_BarycentricCoords ) );
                                    }
                                    else
                                    {
                                        cd.AddCP( rh.m_Point,
                                                  p_0,
                                                  mal::SafeNormalized(p_0-rh.m_Point),
                                                  mal::Norm(p_0-rh.m_Point),
                                                  patch.m_Area );
                                        // mal::Sqrt(patch.m_AreaSq) );
                                        cd.AddPOF( PointOnFeature( feature_id(eFT_Tetrahedron,rh.m_FeatureId.AsTetrahedron()), rh.m_Extra_BarycentricCoords ),
                                                   PointOnFeature( feature_id(eFT_Tetrahedron,eid), p_b ) );
                                    }
                                    GEO_NP_STAT_INC( mig2015.m_Num_CP );
                                }
                            }
                            break;
                        default: break;
                        }
                    }
                    else
                    {
                        //Ignoring empty patch, should not happen!!
                    }
                }
            }
            else
            {
                //Ignoring empty IB, should NOT BE ROOT
            }
        }
    }
    if( bLog ) { GEO_LOG("#CP = %u",cd.Size() ); }

#ifdef __DISABLED_WHILE_DEVELOPING
    // 3) Build IB-graph TODO: 2D only by now

    // 4) Compute IR TODO: 2D first

    // 4) Global Untangling (GU) TODO: 2D first

    // 4) Correspondences TEMP: 2D only by now

    // 5) Generate CP from EC data
#endif //__DISABLED_WHILE_DEVELOPING

    cd.SetNumDisjointManifolds( 1 ); //\todo #IR??
    return cd.End();
}

}} //namespace geo::np
