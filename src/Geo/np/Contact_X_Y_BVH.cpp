#include "Contact_X_Y_BVH.h"
#include <Geo/shape/MeshSolidShape2.h>
#include "Intersection.h"
#include "RayCast.h"
#include "Overlap.h"
#include <Geo/bv/GBDOP.h>

#ifdef __GEO_ENABLE_NP_SCRATCHPAD
#  include <util/ScopedAllocator.h>
#endif

#include <boost/bind.hpp>

#include "Contact_X_Y_BruteForce.h" // Fallback methods if no BVH/DCR

//#define __USE_DCR_V_CP //\note This enables detailed per-DCR.V CP, as opposed to aggregated DCR.E CP, which is the preferred strategy

#define __ENABLE_FLOODIB_WHOLE_NONCROSSING_PATCHES //Optimization, should be enabled by default
#define __USE_RAYCAST_DCR_E //Optimization, should be enabled by default
#define __ENABLE_INTERSECTION_BDT_PRETRANSFORM //\todo SHOULD not be used, DCR.V should be transformed on the fly
#define __ENABLE_INTERSECTION_BDT_SPLIT3 //\todo THIS REQUIRES Context::m_DCR2DCR_E2E_BDT_MaxLeafTests > 0 or an infinite loop occurs
#define __ENABLE_INTERSECTION_BDT_BDOP3_PRESORTED //Optimization, pre-sorts DCR.E positions once before recursive intersection

namespace geo {
namespace np {

//! Transform a plane equation \todo Generalize to GTransformPlane<D> 2D/3D and move elsewhere
inline void TransformPlane( const Vec2& plane_n, Real plane_d,
                            const Transform2& tr,
                            Vec2& transformed_n, Real& transformed_d )
{
    Vec2 plane_p( -plane_d*plane_n );
    transformed_n = tr.m_Rot * plane_n;
    Vec2 transformed_p( tr * plane_p );
    transformed_d = -mal::Dot(transformed_n,transformed_p);
}

//--------------------------------------------------------------------------------------------------------------------------------
// MeshSolid2 Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

/*! Mesh Vs Plane */
bool Contact_DCRMS2_Plane2_BVH_MSS_Lazy_DCR( const MeshSolidShape2* p_mesh, const Transform2& mesh_tr, const Real* vec_dof,
                                             const Vec2& plane_normal, Real plane_coeff_d,
                                             ContactData2& cd, ContactCache2* p_cc,
                                             const Context* p_context )
{
    // GEO_LOG("Contact_DCRMS2_Plane2_BVH_MSS_Lazy_DCR()");
    cd.Begin();

    const DCR_MeshSolidShape2* pDCR( p_mesh->GetDCR() );
    if( 0 == pDCR )
    {
        GEO_LOG_WARNING("Contact_DCRMS2_Plane2_BVH_MSS_Lazy_DCR() with no DCR, falling back to Contact_MeshSolid2_Plane2_BruteForce");
        return Contact_MeshSolid2_Plane2_BruteForce( p_mesh, mesh_tr, vec_dof, plane_normal, plane_coeff_d, cd, p_cc, p_context );
    }

    //TEMP: By now, we'll just Refit the BVH unconditionally if present in the MSS, \todo See Contact_DCRMS2_DCRMS2_BVH_MSS_Lazy_DCR() for detailed discussion
    BVH_MeshSolidShape2* pBVH( p_mesh->GetBVH() );
    if( 0 == pBVH )
    {
        GEO_LOG_WARNING("Contact_DCRMS2_Plane2_BVH_MSS_Lazy_DCR() with no BVH, falling back to Contact_MeshSolid2_Plane2_BruteForce");
        return Contact_MeshSolid2_Plane2_BruteForce( p_mesh, mesh_tr, vec_dof, plane_normal, plane_coeff_d, cd, p_cc, p_context );
    }

    //\todo THIS IS dangerous... Should have a generic SRV type with safe casts...
    const Vec2* vec_sdof( reinterpret_cast<const Vec2* >(vec_dof) );
    const Vec2* default_sdof( p_mesh->GetVecDefaultSDOF() );

    /* Refit BVH
       TEMP: ASSUME UP TO DATE!
       \todo See Contact_DCRMS2_DCRMS2_BVH_MSS_Lazy_DCR() for detailed discussion
    switch( p_context->m_BVH_Method )
    {
    case Context::eBVHM_BV_E:
        pBVH->Refit( boost::bind<void>( &GEBV_MeshSolidShape2_E<BVH_MeshSolidShape2::entry_index_type,BVH_MeshSolidShape2::bv_type>,
                                        p_mesh, mesh_tr, vec_sdof, _1, _2) );
        break;
    case Context::eBVHM_BV_BSlab:
        pBVH->Refit( boost::bind<void>( &GEBV_MeshSolidShape2_BSlab<BVH_MeshSolidShape2::entry_index_type,BVH_MeshSolidShape2::bv_type>,
                                        p_mesh, pDCR, mesh_tr, vec_sdof, _1, _2) );
        break;
    case Context::eBVHM_BV_BDOP:
        pBVH->Refit( boost::bind<void>( &GEBV_MeshSolidShape2_BDOP<BVH_MeshSolidShape2::entry_index_type,BVH_MeshSolidShape2::bv_type>,
                                        p_mesh, pDCR, mesh_tr, vec_sdof, _1, _2) );
        break;
    default: break;
    }
    */

    /*\todo If the plane is a halfspace, we can avoid IM by testing
      all DCR-segments in overlapped elements against the halfspace and
      generating final CP directly according to segment-vertices plane
      equation sign (+,0,-), avoiding duplicated CP for vertices
      shared between adjacent segments.
      Elements classified Vs plane/halfspace
      - Outside => Discard
      - Inside => All DCR-features are in contact
      - Piercing => Must detail DCR-features

      Contact reduction:
      - The DCR is assumed >1 orders of magnitude larger than the embedding MSS, which will yield a LOT of CP.
      - CP reduction could be performed in a per-element basis.
      - Assuming a superposition-principle-like contact treatment, the DCR-CP could be averaged using their barycentric coords wrt DCR-E nodes, yielding a SINGLE CP per DCR-E (=> DCR-E-CP)
        - \todo Prove that superposition of barycentrically applied forces/impulses/displacements is equivalent to applying on the barycentric average...?

      - IMPORTANT: POF should be relative to MSS, NOT the DCR, for
        contact response to work directly on MSS.V, NOT on
        DCR.V... otherwise the dynamics code will need to reprocess
        all DCR-CP and compute their DCR-E-CP equivalents...

      - For completely interior DCR-E penetrating the plane/other object, this average could be precomputed and stored in the DCR-E description.
      - Other per-DCR-E aggregate magnitudes will be required, such as the contact length/area/volume represented by a CP (for pressure/traction computations)
      - Penetration depth will not be accurate if averaged independently of the separation normal. Consider instead averaging the separation displacement vectors, direction and magnitude at once.
      - COULD Intersection Mapping be performed on the REDUCED DCR-E-CP, instead of on the detailed DCR-CP?? This would be a MAJOR simplification...
        - If we implement IM by surface parametrization, surface topology may be required. In 3D the topology of DCR-E-CP will be the same as DCR-E topology, regardless of the detailed DCR-T topology...
      - Should intersection points be "always" CP? in 2D, there will be a few IP, but in 3D there may be a lot, consider reducing them too.

      Roadmap:
      1) DCR vs Plane
         1.1) DCR.T.CP => LOTS of CP
         1.2) DCR.P.CP => 1 CP per Patch
         1.3) Compare results visually and numerically
    */

    // Test overlap

    /* \todo We should compute the actual set of DCR.E that overlap
       the halfspace. It's NOT a good a idea to use the BVH here, as
       the planeBV will be often infinite (if unbounded and not
       axis-aligned). By now, we just consider that all DCR.E overlap
       the plane.
    */
    std::vector< BVH_MeshSolidShape2::entry_index_type > vecOverlaps;
    for( unsigned int it_e=0; it_e<pDCR->m_NumElements; it_e++ )
        vecOverlaps.push_back( it_e );

    // Transform plane to mesh-local coords
    Vec2 plane_normal_1;
    Real plane_coeff_d_1;
    TransformPlane( plane_normal, plane_coeff_d, mal::Inverse(mesh_tr), plane_normal_1, plane_coeff_d_1 );

    // Process overlaps
    for( auto eid : vecOverlaps )
    {
        const DCR_MeshSolidShape2::Element& ed( pDCR->m_vecE[eid] );
        // Get element nodes in mesh-local coords
        uint32 vec_nid[3];  //\todo NID could be stored in DCR::ED if not directly available in MSS
        p_mesh->P_VecVID(eid,vec_nid,3);
        // Test overlap DCR.E vs Plane in mesh-local coords
        // \todo THIS WOULD BE MUCH TIGHTER if we used BDOP/BSlab
        // vertices, and more DCR.E could be classified as interior
        // (much cheaper to handle than piercing DCR.E)
        Real vec_dist[3] = { mal::Dot( plane_normal_1, vec_sdof[vec_nid[0]] ) + plane_coeff_d_1,
                             mal::Dot( plane_normal_1, vec_sdof[vec_nid[1]] ) + plane_coeff_d_1,
                             mal::Dot( plane_normal_1, vec_sdof[vec_nid[2]] ) + plane_coeff_d_1 };
        if( vec_dist[0] <= 0 || vec_dist[1] <= 0 || vec_dist[2] <= 0 ) // interior | piercing
        {
            // Get barycentric transform
            Mat3x3 Bm( 1, 1, 1,
                       default_sdof[vec_nid[0]].x(), default_sdof[vec_nid[1]].x(), default_sdof[vec_nid[2]].x(),
                       default_sdof[vec_nid[0]].y(), default_sdof[vec_nid[1]].y(), default_sdof[vec_nid[2]].y() );
            Mat3x3 invBm( mal::Inverse( Bm ) );
            if( vec_dist[0] <= 0 && vec_dist[1] <= 0 && vec_dist[2] <= 0 ) //interior
            {
#ifdef __USE_DCR_V_CP
                // Add all DCR.E.V as CP, using their barycentric coords, as POF relative to the DCR.Element, NOT to the DCR.Segment
                for( unsigned int it_v=ed.m_FirstVID; it_v<ed.m_FirstVID+ed.m_NumVertices; it_v++ )
                {
                    // Compute DCR.V barycentric coords
                    Vec3 p_b( invBm * mal::Concat(1,pDCR->m_vecV[it_v]) ); //\todo could be precomputed!
                    // Compute DCR.V global coords from node positions
                    Vec2 p_1( p_b[0] * vec_sdof[vec_nid[0]]
                              + p_b[1] * vec_sdof[vec_nid[1]]
                              + p_b[2] * vec_sdof[vec_nid[2]] );
                    Vec2 p_0( mesh_tr * p_1 );
                    Real dist( mal::Dot(plane_normal,p_0) + plane_coeff_d );
                    GEO_ASSERT( dist <= 0 );
                    cd.AddCP( p_0, //on mesh
                              p_0 - dist*plane_normal, //on plane
                              plane_normal, -dist,
                              // DIST IS WRONG!! should be 1/2 DCR.Segment!
                              Real(0.5) * ( mal::Norm( p_1 - p_mesh->V_Pos( p_mesh->HE_OriginVID( p_mesh->HE_Prev(it_he) ), vec_sdof ) )
                                            + mal::Norm( p_1 - p_mesh->V_Pos( p_mesh->HE_FinalVID( it_he ), vec_sdof ) ) ) ); //\todo Radius is measured in DCR, not DCR.E, may be incorrect...
                    cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,eid), mal::Concat(p_b,0) ),
                               PointOnFeature() );
                }
#else //__USE_DCR_V_CP
                /* DCR-E-CP: Single CP with barycentrically
                   averaged pos/depth/normal from DCR.S and summed radius (could
                   be precomputed? NO because deformation changes it
                   nontrivially, but maybe precomputed pose patch
                   length would work fine... try it)
                */
                for( unsigned int it_patch=0; it_patch < ed.m_NumPatches; it_patch++ )
                {
                    Vec3 avg_p_b(0);
                    Vec2 avg_p_0(0);
                    Real total_length(0); //\todo consider using precomputed pose total_length
                    const geo::DCR_MeshSolidShape2::Patch& pd( pDCR->m_vecP[ed.m_FirstPID + it_patch] );
                    for( unsigned int it_sid=0; it_sid < pd.m_NumSegments; it_sid++ )
                    {
                        unsigned int sid( pd.m_FirstSID + it_sid );
                        uint32 vid0( pDCR->m_vecS[sid].GetVID(0) );
                        uint32 vid1( pDCR->m_vecS[sid].GetVID(1) );
                        // Compute barycentric coords for p,q
                        Vec3 p_b( invBm * mal::Concat(1,pDCR->m_vecV[vid0]) ); //\todo could be precomputed!
                        Vec3 q_b( invBm * mal::Concat(1,pDCR->m_vecV[vid1]) ); //\todo could be precomputed!
                        // Compute DCR.V global coords from node positions \todo IMPORTANT: These are NOT global, we have NOT YET applied mesh_tr
                        Vec2 p_0( p_b[0] * vec_sdof[vec_nid[0]] + p_b[1] * vec_sdof[vec_nid[1]] + p_b[2] * vec_sdof[vec_nid[2]] );
                        Vec2 q_0( q_b[0] * vec_sdof[vec_nid[0]] + q_b[1] * vec_sdof[vec_nid[1]] + q_b[2] * vec_sdof[vec_nid[2]] );
                        // Length-weighted average pos (\todo consider using precomputed pose lengths)
                        Real length( mal::Norm( p_0 - q_0 ) );
                        avg_p_b += length * Real(0.5) * (p_b + q_b);
                        avg_p_0 += length * Real(0.5) * (p_0 + q_0);
                        total_length += length;
                    }
                    if( total_length > 0 )
                    {
                        Real rcp_total_length( mal::Rcp(total_length) );
                        avg_p_b = avg_p_b * rcp_total_length;
                        avg_p_0 = mesh_tr * (avg_p_0 * rcp_total_length); //globalize using mesh_tr
                        Real avg_dist( mal::Dot(plane_normal,avg_p_0) + plane_coeff_d );
                        cd.AddCP( avg_p_0, //on mesh
                                  avg_p_0 - avg_dist*plane_normal, //on plane
                                  plane_normal, -avg_dist,
                                  total_length );
                        cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,eid), mal::Concat(avg_p_b,0) ),
                                   PointOnFeature() );
                    }
                }
#endif
            }
            else //piercing
            {
                //\todo CONSIDER transforming plane into element-local coords (barycentric transform and plane are affine, SHOULD WORK)
                // Transform all DCR.E.V
                // \todo THIS MAY BE a uselss optimization, as we'll need the barycentric coords anyway (for piercing segments, at least) and the global pos is easily computed from there...
                util::ScopedAllocator scoped_allocator( p_context->m_ScratchPad, "Contact_DCRMS2_DCRMS2_BruteForce_DCR" );
                Vec2* vec_v1_0 = scoped_allocator.NewArrayPOD<Vec2>(ed.m_NumVertices);
                {
                    Mat2x2 B( mal::GMat2x2_From_Columns(vec_sdof[vec_nid[1]]-vec_sdof[vec_nid[0]],
                                                        vec_sdof[vec_nid[2]]-vec_sdof[vec_nid[0]]) );
                    Mat2x2 Bm( mal::GMat2x2_From_Columns(default_sdof[vec_nid[1]]-default_sdof[vec_nid[0]],
                                                         default_sdof[vec_nid[2]]-default_sdof[vec_nid[0]]) ); //\todo Could be precomputed in DCR::ED
                    Mat2x2 B_invBm( B * mal::Inverse(Bm) );
                    Transform2 tr_B_invBm( mesh_tr.m_Pos, mesh_tr.m_Rot * B_invBm );
                    Vec2 p0_0( mesh_tr.m_Rot * vec_sdof[vec_nid[0]] );
                    // Transform ed1 m_NumVertices vertices from m_First and store them in 0-based indices inside vec_v1_0
                    for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
                        vec_v1_0[it_vie] = tr_B_invBm * (pDCR->m_vecV[ed.m_FirstVID + it_vie] - default_sdof[vec_nid[0]]) + p0_0;
                    GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed.m_NumVertices );
                }
#ifdef __USE_DCR_V_CP
                // Add all intersecting and interior DCR.E.V as CP
                for( unsigned int it_patch=0; it_patch < ed.m_NumPatches; it_patch++ )
                {
                    const geo::DCR_MeshSolidShape2::Patch& pd( pDCR->m_vecP[ed.m_FirstPID + it_patch] );
                    for( unsigned int it_sid=0; it_sid < pd.m_NumSegments; it_sid++ )
                    {
                        unsigned int sid( pd.m_FirstSID + it_sid );
                        uint32 vid0( pDCR->m_vecS[sid].GetVID(0) );
                        uint32 vid1( pDCR->m_vecS[sid].GetVID(1) );
                        Vec2 p_0( vec_v1_0[ vid0 - ed.m_FirstVID ] );
                        Vec2 q_0( vec_v1_0[ vid1 - ed.m_FirstVID ] );
                        Real lambda1, lambda2;
                        //determine inside/outside/piercing and add CP accordingly
                        RayHit2 rh;
                        if( GRayCast_HalfSpace<2>( p_0, q_0-p_0, Interval(0,1), plane_normal, plane_coeff_d, rh, 0 ) ) //\todo <2> required because GVec<int> but GRayHit<unsigned>... decide proper type for Dimension and BE CONSISTENT!
                        {
                            //IMPORTANT: FOR TOTALLY interior segments, we're adding endpoint TWICE!!! analyze min/max and avoid it! RADIUS will also depend on this case!
                            // Compute min/max contact points m,n \todo SHOULD NOT be coincident for a halfspace, but if ray is parallel to it may happen...
                            Vec2 m_0( p_0 + rh.m_Interval.Min() * (q_0-p_0) );
                            Vec2 n_0( p_0 + rh.m_Interval.Max() * (q_0-p_0) );
                            Real dist_m( mal::Dot(m_0,plane_normal)+plane_coeff_d ); //<0
                            Real dist_n( mal::Dot(n_0,plane_normal)+plane_coeff_d ); //<0
                            Real length( mal::Norm( m_0 - n_0 ) );
                            cd.AddCP( m_0,
                                      m_0 - dist_m*plane_normal,
                                      plane_normal, -dist_m,
                                      Real(0.5) * length ); //\todo COMPUTE CP RADIUS
                            cd.AddCP( n_0,
                                      n_0 - dist_n*plane_normal,
                                      plane_normal, -dist_n,
                                      Real(0.5) * length ); //\todo COMPUTE CP RADIUS
                            //IMPORTANT: POF are relative to the MSS, NOT the DCR
                            // Compute barycentric coords for p,q
                            Vec3 p_b( invBm * mal::Concat(1,pDCR->m_vecV[vid0]) ); //\todo could be precomputed!
                            Vec3 q_b( invBm * mal::Concat(1,pDCR->m_vecV[vid1]) ); //\todo could be precomputed!
                            // Add POF with linearly interpolated barycentric coords \todo CHECK IF CORRECT
                            Vec3 m_b( (1-rh.m_Interval.Min()) * p_b + rh.m_Interval.Min() * q_b );
                            Vec3 n_b( (1-rh.m_Interval.Max()) * p_b + rh.m_Interval.Max() * q_b );
                            //\todo m_b,n_b could also be computed using deformed barycentric transform invB as invB*m_0, invB*n_0
                            cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,eid), mal::Concat(m_b,0) ),
                                       PointOnFeature() );
                            cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,eid), mal::Concat(n_b,0) ),
                                       PointOnFeature() );
                        }
                    }
                }
#else
                // DCR-E-CP: Find intersecting DCR.P, compute weighted sums...
                for( unsigned int it_patch=0; it_patch < ed.m_NumPatches; it_patch++ )
                {
                    Vec3 avg_p_b(0);
                    Vec2 avg_p_0(0);
                    Real total_length(0);
                    const geo::DCR_MeshSolidShape2::Patch& pd( pDCR->m_vecP[ed.m_FirstPID + it_patch] );
                    for( unsigned int it_sid=0; it_sid < pd.m_NumSegments; it_sid++ )
                    {
                        unsigned int sid( pd.m_FirstSID + it_sid );
                        uint32 vid0( pDCR->m_vecS[sid].GetVID(0) );
                        uint32 vid1( pDCR->m_vecS[sid].GetVID(1) );
                        Vec2 p_0( vec_v1_0[ vid0 - ed.m_FirstVID ] );
                        Vec2 q_0( vec_v1_0[ vid1 - ed.m_FirstVID ] );
                        Real lambda1, lambda2;
                        //determine inside/outside/piercing and add CP accordingly
                        RayHit2 rh;
                        if( GRayCast_HalfSpace<2>( p_0, q_0-p_0, Interval(0,1), plane_normal, plane_coeff_d, rh, 0 ) ) //\todo <2> required because GVec<int> but GRayHit<unsigned>... decide proper type for Dimension and BE CONSISTENT!
                        {
                            // Compute min/max contact points m,n \todo SHOULD NOT be coincident for a halfspace, but if ray is parallel to it may happen...
                            Vec2 m_0( p_0 + rh.m_Interval.Min() * (q_0-p_0) );
                            Vec2 n_0( p_0 + rh.m_Interval.Max() * (q_0-p_0) );
                            // Compute barycentric coords for p,q
                            Vec3 p_b( invBm * mal::Concat(1,pDCR->m_vecV[vid0]) ); //\todo could be precomputed!
                            Vec3 q_b( invBm * mal::Concat(1,pDCR->m_vecV[vid1]) ); //\todo could be precomputed!
                            //\todo m_b,n_b could also be computed using deformed barycentric transform invB as invB*m_0, invB*n_0
                            Vec3 m_b( (1-rh.m_Interval.Min()) * p_b + rh.m_Interval.Min() * q_b );
                            Vec3 n_b( (1-rh.m_Interval.Max()) * p_b + rh.m_Interval.Max() * q_b );
                            // Length-weighted average pos
                            Real length( mal::Norm( m_0 - n_0 ) );
                            avg_p_b += length * Real(0.5) * (m_b + n_b);
                            avg_p_0 += length * Real(0.5) * (m_0 + n_0);
                            total_length += length;
                        }
                    }
                    if( total_length > 0 )
                    {
                        Real rcp_total_length( mal::Rcp(total_length) );
                        avg_p_b = avg_p_b * rcp_total_length;
                        avg_p_0 = avg_p_0 * rcp_total_length;
                        Real avg_dist( mal::Dot(plane_normal,avg_p_0) + plane_coeff_d );
                        cd.AddCP( avg_p_0, //on mesh
                                  avg_p_0 - avg_dist*plane_normal, //on plane
                                  plane_normal, -avg_dist,
                                  total_length );
                        cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,eid), mal::Concat(avg_p_b,0) ),
                                   PointOnFeature() );
                    }
                }
#endif
            }
        }
        // else exterior
    }
    if( cd.Size() > 0 ) cd.SetNumDisjointManifolds( 1 ); //\todo Count it properly, for a deformable object it may be > 1
    return cd.End();
}

bool TestOverlap_BSlabs( const MeshSolidShape2* p_mesh1, const DCR_MeshSolidShape2* p_dcr1, const Transform2& tr1, const MeshSolidShape2::sdof_type* vec_sdof1,
                         BVH_MeshSolidShape2::entry_index_type eid1, const BVH_MeshSolidShape2::bv_type& bv1,
                         const MeshSolidShape2* p_mesh2, const DCR_MeshSolidShape2* p_dcr2, const Transform2& tr2, const MeshSolidShape2::sdof_type* vec_sdof2,
                         BVH_MeshSolidShape2::entry_index_type eid2, const BVH_MeshSolidShape2::bv_type& bv2 )
{
    GEO_ASSERT( eid1 < p_dcr1->m_NumElements && eid2 < p_dcr2->m_NumElements );
    /* BSlab => GSlab
       - A BSlab defines a barycentric range BDOP[i] along a pose-local axis pose_axis[i] perpendicular to an element face in pose (pose_v0,pose_v1,pose_v2)
       - Pose-Slab (BSlab) => Global-Slab (GSlab)
         - The BSlab REMAINS a slab after barycentric transformation, because parallelism is preserved under affine transforms \todo PROVE/REFERENCE!
         - The global_axis[i] REMAINS perpendicular to the element face in global coords (global_v0,global_v1,global_v2) \todo PROVE/REFERENCE!
         - The global_range[i] changes... but how?!
           - It can always be computed from the 4/6 points that defime the global-slab as in GEBV
           - Due to parallelism, the range is the same for all near points and for all far points => only 2 global slab points need to be projected, 1 Min and 1 Max
           - MAYBE the range can be transformed "without projections" using some affine property that only depends on the global_axis[i] direction... and the pose-BDOP[i] range

       Potential tests:
       a) GSlab1 vs BV2
       b) GSlab1 vs Tri2/Tetra2
          - Requires projecting its 3/4 global vertices
       c) GSlab1 vs GSlab2_Vertices
          - Requires computing and projecting its 4/6 global vertices
       d) GSlab1 vs GBDOP2_Vertices
          - Requires computing and projecting ALL GBDOP global vertices
       e) Exact GSlab1 vs GSlab2:
          - Compute all global-slab vertices
          - Compute all SAT axis 3+3 in 2D, 4+4+6*6 in 3D
          - Test all SAT axis => project 3/4 vertices for each
          => TOO EXPENSIVE!!!!!
       f) Approx GSlab1 vs GSlab2:
          - Compute all global-slab vertices
          - Compute face-SAT axis 3+3 in 2D, 4+4 in 3D
          - Test face-SAT axis => project 3/4 vertices for each
       h) Exact GDOP1 vs GDOP2
          - Most exact test possible from BDOP
          - Requrires all BDOP edges (topology could be cached, but
            storage increases considerably)
          - Requires insane amount of SAT tests
          ==> TOO EXPENSIVE!!!
       g) BDOP1 vs GSlab2_Vertices:
          - Any pose-slab BDOP1[i] can ALSO be tested against
            GSlab2_Vertices[best_bslab_index] transformed into
            BaricentricCoords1... this avoids having to convert slabs to
            global coords!
          - Also, instead of transforming 2*(D-1)
            GSlab2_Vertices[best_bslab_index] we can just transform the
            D+1 nodes and compute the 2*(D-1) GSlab2_Vertices in
            BaricentricCoords1

       Analysis:
       - We assume that most elements in a DCR have Volume(slab) << Volume(element)
       - Discarded options
         - a) and b) are quite conservative and should be
           performed symmetrically (BV1|Tri1/Tetra1,GSlab2) to get significant
           volume culling.
         - d) will cull tightly, but computing GBDOP_Vertices is VERY expensive.
         - e) is VERY expensive, specially in 3D
         - h) is intractable...
       - Valid options:
         - c) may cull quite efficiently on its own, as it does consider GSlab2.
              - It may be repeated for EACH GSlab1 in BDOP1[0..D] reusing the GSlab2_Vertices[best_bslab_index]
              - This tests D+1 infinite slabs GSlab1[i] against a SINGLE finite volume GSlab2_Vertices[best_bslab_index]
         - f) is the global version of g) and exactly the same as doing c) symmetrically for all BDOP1 and BDOP2 axis
         - g) seems quite tight, as it uses exact BDOP1 in local
              coords and simplified GSlab2_Vertices. Avoids computing
              ANY GSlab axis, as opposed to a)..f) which requires
              cross-products in 3D, but requires inverse barycentric
              transforms...

         => Therefore, we choose either f) or g)
            f) works in Global coords and can leverage cached
               transformed nodes (may be available from Refit()!)
               - Requires all 3/4 GSlab axis, which are element face
                 normals (compute all edge vectors + 4 cross-products
                 in 3D)
               - Requires all 3/4 GSlab ranges, each from 2 projected
                 GSlab_Vertices that need to be computed from global nodes.
            g) Works with global node coords (may be available)
               - Requires inverse barycentric transforms (D vector sub + DxD matrix inversion)
               - GSlab_Vertices are directly computed in barycentric-local coords.

            => g) is more "elegant", it requires less low-level code,
               but f) allows incremental computation of the required
               magnitudes and can be aborted as soon as a low-level
               test fails. (eg: GSlab axis and ranges can be computed
               and tested one by one, as opposed to B^-1, which is
               computed atomically even if the first local BSlab range
               test fails... I'll implement g) but f) might be faster.
    */
    //------- Option (g)
    //--- Elements
    const DCR_MeshSolidShape2::Element& ed1( p_dcr1->m_vecE[eid1] );
    const DCR_MeshSolidShape2::Element& ed2( p_dcr2->m_vecE[eid2] );
    //--- global nodes
    uint32 vec_nid1[3];
    p_mesh1->P_VecVID( eid1, vec_nid1, 3 ); //\todo NID could be stored in DCR::ED if not directly available in MSS
    Vec2 vec_node_pos1_0[3] = { tr1 * vec_sdof1[vec_nid1[0]], tr1 * vec_sdof1[vec_nid1[1]], tr1 * vec_sdof1[vec_nid1[2]] };
    uint32 vec_nid2[3];
    p_mesh2->P_VecVID( eid2, vec_nid2, 3 ); //\todo NID could be stored in DCR::ED if not directly available in MSS
    Vec2 vec_node_pos2_0[3] = { tr2 * vec_sdof2[vec_nid2[0]], tr2 * vec_sdof2[vec_nid2[1]], tr2 * vec_sdof2[vec_nid2[2]] };
    //--- GDOP1 vs GSlab2_Vertices
    {
        // Barycentric transform B1 using node0 as the local origin
        Mat2x2 B1( mal::GMat2x2_From_Columns( vec_node_pos1_0[1]-vec_node_pos1_0[0], vec_node_pos1_0[2]-vec_node_pos1_0[0] ) );
        Mat2x2 invB1( mal::Inverse(B1) );
        // Transform eid2 nodes to B1 coords (u,v,t)
        Vec2 vec_node_vt[3] = { invB1*(vec_node_pos2_0[0]-vec_node_pos1_0[0]),
                                invB1*(vec_node_pos2_0[1]-vec_node_pos1_0[0]),
                                invB1*(vec_node_pos2_0[2]-vec_node_pos1_0[0]) };
        Vec3 vec_node_pos2_B1[3] = { Vec3( 1-vec_node_vt[0][0]-vec_node_vt[0][1], vec_node_vt[0][0], vec_node_vt[0][1] ),
                                     Vec3( 1-vec_node_vt[1][0]-vec_node_vt[1][1], vec_node_vt[1][0], vec_node_vt[1][1] ),
                                     Vec3( 1-vec_node_vt[2][0]-vec_node_vt[2][1], vec_node_vt[2][0], vec_node_vt[2][1] ) };
        /* \todo HERE we could test BDOP[i] using the vec_node_pos2_B1
           as an early-out to avoid computing the bslab_pos2_B1, but
           seems overkill */
        // Compute BSlab2 vertices in B1 coords
        int bsi( ed2.m_BDOP_BestSlabIdx );
        Vec3 node_pos_i_times_min( ed2.m_BDOP[bsi].Min() * vec_node_pos2_B1[bsi] );
        Real one_minus_min( 1 - ed2.m_BDOP[bsi].Min() );
        Vec3 node_pos_i_times_max( ed2.m_BDOP[bsi].Max() * vec_node_pos2_B1[bsi] );
        Real one_minus_max( 1 - ed2.m_BDOP[bsi].Max() );
        Vec3 bslab_pos2_B1[4] = { node_pos_i_times_min + one_minus_min * vec_node_pos2_B1[(bsi+1)%3],
                                  node_pos_i_times_min + one_minus_min * vec_node_pos2_B1[(bsi+2)%3],
                                  node_pos_i_times_max + one_minus_max * vec_node_pos2_B1[(bsi+1)%3],
                                  node_pos_i_times_max + one_minus_max * vec_node_pos2_B1[(bsi+2)%3] };
        // Test BDOP[i], exit if no overlap
        for( int i=0; i<3; i++ )
        {
            Interval interval2_B1( mal::Min( mal::Min( bslab_pos2_B1[0][i], bslab_pos2_B1[1][i] ),
                                             mal::Min( bslab_pos2_B1[2][i], bslab_pos2_B1[3][i] ) ),
                                   mal::Max( mal::Max( bslab_pos2_B1[0][i], bslab_pos2_B1[1][i] ),
                                             mal::Max( bslab_pos2_B1[2][i], bslab_pos2_B1[3][i] ) ) );
            if( !ed1.m_BDOP[i].TestOverlap( interval2_B1 ) )
                return false;
        }
    }
    //--- GDOP2 vs GSlab1_Vertices
    {
        // Barycentric transform B2 using node0 as the local origin
        Mat2x2 B2( mal::GMat2x2_From_Columns( vec_node_pos2_0[1]-vec_node_pos2_0[0], vec_node_pos2_0[2]-vec_node_pos2_0[0] ) );
        Mat2x2 invB2( mal::Inverse(B2) );
        // Transform eid1 nodes to B2 coords (u,v,t)
        Vec2 vec_node_vt[3] = { invB2*(vec_node_pos1_0[0]-vec_node_pos2_0[0]),
                                invB2*(vec_node_pos1_0[1]-vec_node_pos2_0[0]),
                                invB2*(vec_node_pos1_0[2]-vec_node_pos2_0[0]) };
        Vec3 vec_node_pos1_B2[3] = { Vec3( 1-vec_node_vt[0][0]-vec_node_vt[0][1], vec_node_vt[0][0], vec_node_vt[0][1] ),
                                     Vec3( 1-vec_node_vt[1][0]-vec_node_vt[1][1], vec_node_vt[1][0], vec_node_vt[1][1] ),
                                     Vec3( 1-vec_node_vt[2][0]-vec_node_vt[2][1], vec_node_vt[2][0], vec_node_vt[2][1] ) };
        /* \todo HERE we could test BDOP[i] using the vec_node_pos1_B2
           as an early-out to avoid computing the bslab_pos1_B2, but
           seems overkill */
        // Compute BSlab1 vertices in B2 coords
        int bsi( ed1.m_BDOP_BestSlabIdx );
        Vec3 node_pos_i_times_min( ed1.m_BDOP[bsi].Min() * vec_node_pos1_B2[bsi] );
        Real one_minus_min( 1 - ed1.m_BDOP[bsi].Min() );
        Vec3 node_pos_i_times_max( ed1.m_BDOP[bsi].Max() * vec_node_pos1_B2[bsi] );
        Real one_minus_max( 1 - ed1.m_BDOP[bsi].Max() );
        Vec3 bslab_pos1_B2[4] = { node_pos_i_times_min + one_minus_min * vec_node_pos1_B2[(bsi+1)%3],
                                  node_pos_i_times_min + one_minus_min * vec_node_pos1_B2[(bsi+2)%3],
                                  node_pos_i_times_max + one_minus_max * vec_node_pos1_B2[(bsi+1)%3],
                                  node_pos_i_times_max + one_minus_max * vec_node_pos1_B2[(bsi+2)%3] };
        // Test BDOP[i], exit if no overlap
        for( int i=0; i<3; i++ )
        {
            Interval interval1_B2( mal::Min( mal::Min( bslab_pos1_B2[0][i], bslab_pos1_B2[1][i] ),
                                             mal::Min( bslab_pos1_B2[2][i], bslab_pos1_B2[3][i] ) ),
                                   mal::Max( mal::Max( bslab_pos1_B2[0][i], bslab_pos1_B2[1][i] ),
                                             mal::Max( bslab_pos1_B2[2][i], bslab_pos1_B2[3][i] ) ) );
            if( !ed2.m_BDOP[i].TestOverlap( interval1_B2 ) )
                return false;
        }
    }
    return true;
}

struct crossing_segment_pair_type
{
    PointOnFeature m_POF[2];
    uint32 m_IBID[2]; //IB_1, IB_2
    uint32 m_EID[2]; //POF EID on each object
    uint32 m_PID[2]; //POF PID on each object
    Vec2 m_a_0[2]; //segment begin-points
    Vec2 m_b_0[2]; //segment end-points
    inline crossing_segment_pair_type( const PointOnFeature& pof1, uint32 eid1, uint32 pid1, const Vec2& a1_0, const Vec2& b1_0,
                                       const PointOnFeature& pof2, uint32 eid2, uint32 pid2, const Vec2& a2_0, const Vec2& b2_0 )
        { m_POF[0] = pof1; m_POF[1] = pof2;
          m_IBID[0] = -1; m_IBID[1] = -1; //unknown
          m_EID[0] = eid1; m_EID[1] = eid2;
          m_PID[0] = pid1; m_PID[1] = pid2;
          m_a_0[0] = a1_0; m_a_0[1] = a2_0;
          m_b_0[0] = b1_0; m_b_0[1] = b2_0; }
};

/* AddCFP_Intersection_DCR2_E_DCR2_E_BruteForce()

   Add crossing_segment_pair_type to given vec_cfp, found by
   bruteforce DCR.E vs DCR.E all Patch.Segment pairwise tests.

  \todo OPTIMIZATION: Consider caching LRU per-element
        transformed vertices to avoid retransforming its contents if
        it appears in several (e1,e2) pairs, however, caching is
        incompatible with ScopedAllocator...
        => Alternatively, we could sort pairs so that consecutive
           pairs may have the same e1 in <e1,e2> and keep e1 cached
           for each matching e2... this WOULD be compatible with the
           ScopedAllocator, but not with the bvh2bvh_overlap_iterator...
        => We could use the ScopedAllocator if the cache has a fixed
           size and we use it "linearly" through push/pop or flush it
           regularly instead of LRU lifetime that requires
           adding/removing any element in the cache at any time.
           => We can double-buffer the cache, so we fill a cache
              and, when it's full, any cache miss is counted and
              added to a secondary cache. Cache miss/hit are
              counted. When the secondary cache is becomes full, we
              flush (THE LEAST USED?) one and swap if necessary.
           - Given 2 scope-allocated memory blocs, each one can be
             managed as cache with *variable*-sized entries as
             follows:
             - Entries are stored sequentially from the beginning of
               the memory block.
             - Cached Entry descriptors
               (entry_id,entry_offset,entry_size) are stored in
               reverse order from the end of the memory block.
             - For each query, the cache needs to search cache
               descriptors and, if not found, add a new one and
               provide the required memory sub-block on cache-miss
               or failure if it does not fit.
           => Instead of filling the "primary" cache and moving to
              the secondary afterwards, it may be better to
              associate entries with caches according to their
              entry_id lower bit(s) and flush or realloc them
              directly on the first non-fitting miss... (consider
              partial-flush by popping only the required entries to
              fit the new missed entry.
              - A different full cache miss strategy would be to
                store the new entry in the LRU cached sub-block
                where it fits... with potential periodic compaction
                thanks to offsets (cannot save pointers, therefore)
                => this is becoming a bloody memory manager...
*/
unsigned int AddCFP_Intersection_DCR2_E_DCR2_E_BruteForce( const MeshSolidShape2* p_mesh1, const Vec2* vec_sdof1, const DCR_MeshSolidShape2* p_dcr1, uint32 eid1,
                                                           const MeshSolidShape2* p_mesh2, const Vec2* vec_sdof2, const DCR_MeshSolidShape2* p_dcr2, uint32 eid2,
                                                           const Transform2& tr1_2, const Transform2& tr2_0,
                                                           std::vector< crossing_segment_pair_type >& vec_cfp,
                                                           const Context* p_context = g_pDefaultContext )
{
    const DCR_MeshSolidShape2::Element& ed1( p_dcr1->m_vecE[eid1] );
    const DCR_MeshSolidShape2::Element& ed2( p_dcr2->m_vecE[eid2] );

    /* This can happen if non-layer[0] elements are present in
     * the DCR and contain no geometry patches
     */
    if( ed1.m_NumPatches == 0 || ed2.m_NumPatches == 0 ) return 0;

    const Vec2* default_sdof1( p_mesh1->GetVecDefaultSDOF() );
    const Vec2* default_sdof2( p_mesh2->GetVecDefaultSDOF() );

    // Transform ALL vertices e1 \todo and cache
    util::ScopedAllocator scoped_allocator_e1( p_context->m_ScratchPad, "Contact_DCRMS2_DCRMS2_BruteForce_MSS_Lazy_DCR_e1" );
    Vec2* vec_v1_2 = scoped_allocator_e1.NewArrayPOD<Vec2>(ed1.m_NumVertices);
    {
        uint32 vec_nid1[3];
        p_mesh1->P_VecVID( eid1, vec_nid1, 3 ); //\todo NID could be stored in DCR::ED if not directly available in MSS
        Mat2x2 B( mal::GMat2x2_From_Columns(vec_sdof1[vec_nid1[1]]-vec_sdof1[vec_nid1[0]],
                                            vec_sdof1[vec_nid1[2]]-vec_sdof1[vec_nid1[0]]) );
        Mat2x2 Bm( mal::GMat2x2_From_Columns(default_sdof1[vec_nid1[1]]-default_sdof1[vec_nid1[0]],
                                             default_sdof1[vec_nid1[2]]-default_sdof1[vec_nid1[0]]) ); //\todo Could be precomputed in DCR::ED
        Mat2x2 B_invBm( B * mal::Inverse(Bm) );
        Transform2 tr1_2_B_invBm( tr1_2.m_Pos, tr1_2.m_Rot * B_invBm );
        Vec2 p0_2( tr1_2.m_Rot * vec_sdof1[vec_nid1[0]] );
        // Transform ed1 m_NumVertices vertices from m_First and store them in 0-based indices inside vec_v1_2
        for( unsigned int it_vie=0; it_vie<ed1.m_NumVertices; it_vie++ )
            vec_v1_2[it_vie] = tr1_2_B_invBm * (p_dcr1->m_vecV[ed1.m_FirstVID + it_vie] - default_sdof1[vec_nid1[0]]) + p0_2;
        GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed1.m_NumVertices );
    }
    // Transform ALL vertices e2 \todo and cache
    util::ScopedAllocator scoped_allocator_e2( p_context->m_ScratchPad, "Contact_DCRMS2_DCRMS2_BruteForce_MSS_Lazy_DCR_e2" );
    Vec2* vec_v2_2 = scoped_allocator_e2.NewArrayPOD<Vec2>(ed2.m_NumVertices);
    {
        uint32 vec_nid2[3];
        p_mesh2->P_VecVID( eid2, vec_nid2, 3 ); //\todo NID could be stored in DCR::ED if not directly available in MSS
        Mat2x2 B( mal::GMat2x2_From_Columns(vec_sdof2[vec_nid2[1]]-vec_sdof2[vec_nid2[0]],
                                            vec_sdof2[vec_nid2[2]]-vec_sdof2[vec_nid2[0]]) );
        Mat2x2 Bm( mal::GMat2x2_From_Columns(default_sdof2[vec_nid2[1]]-default_sdof2[vec_nid2[0]],
                                             default_sdof2[vec_nid2[2]]-default_sdof2[vec_nid2[0]]) ); //\todo Could be precomputed in DCR::ED
        Mat2x2 B_invBm( B * mal::Inverse(Bm) );
        Vec2 p0_2( vec_sdof2[vec_nid2[0]] );
        for( unsigned int it_vie=0; it_vie<ed2.m_NumVertices; it_vie++ )
            vec_v2_2[it_vie] = B_invBm * (p_dcr2->m_vecV[ed2.m_FirstVID + it_vie] - default_sdof2[vec_nid2[0]]) + p0_2;
        GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed2.m_NumVertices );
    }
    //-- Test all patch pairs using cached transformed vertices
    // For each patch1 in e1
    uint32 num_crossings(0);
    for( unsigned int it_patch1=0; it_patch1 < ed1.m_NumPatches; it_patch1++ )
    {
        uint32 pid1( ed1.m_FirstPID + it_patch1 );
        const geo::DCR_MeshSolidShape2::Patch& pd1( p_dcr1->m_vecP[pid1] );
        // For each patch2 in e2
        for( unsigned int it_patch2=0; it_patch2 < ed2.m_NumPatches; it_patch2++ )
        {
            uint32 pid2( ed2.m_FirstPID + it_patch2 );
            const geo::DCR_MeshSolidShape2::Patch& pd2( p_dcr2->m_vecP[pid2] );
            // For each segment (a1,b1)
            for( unsigned int it_sid1=0; it_sid1 < pd1.m_NumSegments; it_sid1++ )
            {
                unsigned int sid1( pd1.m_FirstSID + it_sid1 );
                Vec2 a1_2( vec_v1_2[ p_dcr1->m_vecS[sid1].GetVID(0) - ed1.m_FirstVID ] );
                Vec2 b1_2( vec_v1_2[ p_dcr1->m_vecS[sid1].GetVID(1) - ed1.m_FirstVID ] );
                //\todo IMPORTANT! Early-out if (a1,b1) does NOT overlap bv2!! \todo use bslab2 here, it's FAST to test and provides the smallest volume among all BV types
                // for each segment (a2,b2)
                for( unsigned int it_sid2=0; it_sid2 < pd2.m_NumSegments; it_sid2++ )
                {
                    unsigned int sid2( pd2.m_FirstSID + it_sid2 );
                    Vec2 a2_2( vec_v2_2[ p_dcr2->m_vecS[sid2].GetVID(0) - ed2.m_FirstVID ] );
                    Vec2 b2_2( vec_v2_2[ p_dcr2->m_vecS[sid2].GetVID(1) - ed2.m_FirstVID ] );
                    Real lambda1, lambda2;
                    if( Intersection_Segment2_Segment2( a1_2, b1_2, a2_2, b2_2, lambda1, lambda2 ) )
                    {
                        // save Crossing Feature Pair
                        vec_cfp.push_back( crossing_segment_pair_type( PointOnFeature( feature_id( eFT_Segment, sid1 ), Vec4(Real(1)-lambda1,lambda1,0,0) ), //POF1
                                                                       eid1, pid1, //eid1,pid1
                                                                       tr2_0*a1_2, tr2_0*b1_2, //a1_0,b1_0
                                                                       PointOnFeature( feature_id( eFT_Segment, sid2 ), Vec4(Real(1)-lambda2,lambda2,0,0) ), //POD2
                                                                       eid2, pid2, //eid2,pid2
                                                                       tr2_0*a2_2, tr2_0*b2_2 ) ); //a2_0,b2_0
                        num_crossings++;
                    }
                }
            }
        }
    }
    return num_crossings;
}

/* Stack-based recursive BDT test
   - Computation in p_mesh2 ref sys, but results in global coords
   - BDT1 BDOPs are projected to BDT2 axis to test ranges
   - The overall dual-BDT recursion algorithm is similar to
     GBoundingVolumeHierarchy one, however there's two IMPORTANT
     differences:
     - The leaf case is generalized to #BDTN1.S * #BDTN2.S < max_leaf_tests
     - In a GBVH, recursion only needs to consume/split one of the
       nodes (eg: BDTN1) into L/R subtrees and recurse with the other
       one (BDTN2) unchanged. In a BDT, however, as nodes ARE
       segments, BDTN1 node consumption ACTUALLY consumes its
       representing segment BDTN1.S[0], which MUST be tested against
       the WHOLE BTDN2 subtree:
       a) An INEFFICIENT option would be to perform a BDTN1.S[0] *
          BDTN2.S brute-force test and effectively forget about
          BDTN1.S[0], however, this would NOT FILTER BDTN2 subtree at
          all despite having valid BDOP
       b) A BETTER option seems to split BDTN1 into 3 stacked queries,
          L/R as usual, and a BDTN1.S[0] x BDTN2 query (where the
          stacked BDT1 interval will contain strictly 1 node/segment
       => BOTH methods yield the same results, but (b) is much more
          elegant. Enabled with __ENABLE_INTERSECTION_BDT_SPLIT3

   - \todo Use uses DescendLarger(n1,n2) strategy in RTCD pg.257 (as in GBVH)
     - This requires comparing BDTNode sizes in GLOBAL coords, which
       may be expensive to compute. Instead, we could store BDTNode
       BDOP largest size (quantized) or axis multiply it by the DCR.E
       barycentric transform largest eigenvalue (or similar upper
       bound) to get a rough upper bound of BDTNode global size.

   - \todo OPTIMIZATIONS:
     - Consider pre-transform or not all DCR.E.V.
     - Consider passing DCR.E Bs,invBs,invBm as params instead of
       recomputing, and simplify DCR.E.V transform code accordingly
     - NO: We could pass the non-split BDTN cartesian BV through the
       stack to avoid recomputing it, at the cost of significantly
       increasing stack size, NOT WORTH IT
     - NO: Consider Barycentric test, instead of cartesian, though I
       don't think it'll be worth the effort/CPU
*/
unsigned int AddCFP_Intersection_DCR2_E_DCR2_E_BDT( const MeshSolidShape2* p_mesh1, const Vec2* vec_sdof1, const DCR_MeshSolidShape2* p_dcr1, uint32 eid1,
                                                    const MeshSolidShape2* p_mesh2, const Vec2* vec_sdof2, const DCR_MeshSolidShape2* p_dcr2, uint32 eid2,
                                                    const Transform2& tr1_2, const Transform2& tr2_0,
                                                    std::vector< crossing_segment_pair_type >& vec_cfp,
                                                    const Context* p_context = g_pDefaultContext )
{
    uint32 num_crossings(0);
    const DCR_MeshSolidShape2::Element& ed1( p_dcr1->m_vecE[eid1] );
    const DCR_MeshSolidShape2::Element& ed2( p_dcr2->m_vecE[eid2] );

    /* This can happen if non-layer[0] elements are present in
     * the DCR and contain no geometry patches
     */
    if( ed1.m_NumPatches == 0 || ed2.m_NumPatches == 0 ) return 0;

    const Vec2* default_sdof1( p_mesh1->GetVecDefaultSDOF() );
    const Vec2* default_sdof2( p_mesh2->GetVecDefaultSDOF() );

#ifdef __ENABLE_INTERSECTION_BDT_PRETRANSFORM
    // TEMPORAL: To test recursion we pre-transform ALL DCR.E.V, but
    // the final algorithm should transform them on the fly when
    // retrieving each DCR.S to test intersection!!

    // Transform ALL vertices e1 \todo and cache
    util::ScopedAllocator scoped_allocator_e1( p_context->m_ScratchPad, "Contact_DCRMS2_DCRMS2_BruteForce_MSS_Lazy_DCR_e1" );
    Vec2* vec_v1_2 = scoped_allocator_e1.NewArrayPOD<Vec2>(ed1.m_NumVertices);
    uint32 vec_nid1[3];
    p_mesh1->P_VecVID( eid1, vec_nid1, 3 ); //\todo NID could be stored in DCR::ED if not directly available in MSS
    {
        Mat2x2 B( mal::GMat2x2_From_Columns(vec_sdof1[vec_nid1[1]]-vec_sdof1[vec_nid1[0]],
                                            vec_sdof1[vec_nid1[2]]-vec_sdof1[vec_nid1[0]]) );
        Mat2x2 Bm( mal::GMat2x2_From_Columns(default_sdof1[vec_nid1[1]]-default_sdof1[vec_nid1[0]],
                                             default_sdof1[vec_nid1[2]]-default_sdof1[vec_nid1[0]]) ); //\todo Could be precomputed in DCR::ED
        Mat2x2 B_invBm( B * mal::Inverse(Bm) );
        Transform2 tr1_2_B_invBm( tr1_2.m_Pos, tr1_2.m_Rot * B_invBm );
        Vec2 p0_2( tr1_2.m_Rot * vec_sdof1[vec_nid1[0]] );
        // Transform ed1 m_NumVertices vertices from m_First and store them in 0-based indices inside vec_v1_2
        for( unsigned int it_vie=0; it_vie<ed1.m_NumVertices; it_vie++ )
            vec_v1_2[it_vie] = tr1_2_B_invBm * (p_dcr1->m_vecV[ed1.m_FirstVID + it_vie] - default_sdof1[vec_nid1[0]]) + p0_2;
        GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed1.m_NumVertices );
    }
    // Transform ALL vertices e2 \todo and cache
    util::ScopedAllocator scoped_allocator_e2( p_context->m_ScratchPad, "Contact_DCRMS2_DCRMS2_BruteForce_MSS_Lazy_DCR_e2" );
    Vec2* vec_v2_2 = scoped_allocator_e2.NewArrayPOD<Vec2>(ed2.m_NumVertices);
    uint32 vec_nid2[3];
    p_mesh2->P_VecVID( eid2, vec_nid2, 3 ); //\todo NID could be stored in DCR::ED if not directly available in MSS
    {
        Mat2x2 B( mal::GMat2x2_From_Columns(vec_sdof2[vec_nid2[1]]-vec_sdof2[vec_nid2[0]],
                                            vec_sdof2[vec_nid2[2]]-vec_sdof2[vec_nid2[0]]) );
        Mat2x2 Bm( mal::GMat2x2_From_Columns(default_sdof2[vec_nid2[1]]-default_sdof2[vec_nid2[0]],
                                             default_sdof2[vec_nid2[2]]-default_sdof2[vec_nid2[0]]) ); //\todo Could be precomputed in DCR::ED
        Mat2x2 B_invBm( B * mal::Inverse(Bm) );
        Vec2 p0_2( vec_sdof2[vec_nid2[0]] );
        for( unsigned int it_vie=0; it_vie<ed2.m_NumVertices; it_vie++ )
            vec_v2_2[it_vie] = B_invBm * (p_dcr2->m_vecV[ed2.m_FirstVID + it_vie] - default_sdof2[vec_nid2[0]]) + p0_2;
        GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed2.m_NumVertices );
    }
    // END TEMPORAL
#endif

    // Precomp nodes in 2-refsys for BV(BDOP)
    Vec2 vec_node_pos1_2[] = { tr1_2*vec_sdof1[vec_nid1[0]], tr1_2*vec_sdof1[vec_nid1[1]], tr1_2*vec_sdof1[vec_nid1[2]] };
    Vec2 vec_node_pos2_2[] = { vec_sdof2[vec_nid2[0]], vec_sdof2[vec_nid2[1]], vec_sdof2[vec_nid2[2]] };

#ifdef __ENABLE_INTERSECTION_BDT_BDOP3_PRESORTED
    // Presort node indices along each BV axis (\note For a GKDOP<D,K> each array must allocate K/2 sub-arrays arrays of D node-in-element indices)
    int vec_sorted_nie1[6] = { 0,1,2, 0,1,2 }; //x [0,3) y [3,6)
    int vec_sorted_nie2[6] = { 0,1,2, 0,1,2 }; //x [0,3) y [3,6)
    if( p_context->m_DCR2DCR_E2E_BDT_Presort )
    {
        bv::Compute_AABB2_From_BDOP3_PreSort( vec_node_pos1_2, vec_sorted_nie1 );
        bv::Compute_AABB2_From_BDOP3_PreSort( vec_node_pos2_2, vec_sorted_nie2 );
    }
#endif

    // For all patch pairs
    for( unsigned int it_pie1=0; it_pie1<ed1.m_NumPatches; it_pie1++ )
    {
        uint32 pid1( ed1.m_FirstPID+it_pie1 );
        const DCR_MeshSolidShape2::Patch& pd1( p_dcr1->m_vecP[pid1] );
        for( unsigned int it_pie2=0; it_pie2<ed2.m_NumPatches; it_pie2++ )
        {
            uint32 pid2( ed2.m_FirstPID+it_pie2 );
            const DCR_MeshSolidShape2::Patch& pd2( p_dcr2->m_vecP[pid2] );
            // Init stack with pd1,pd2 root BDTN
            typedef std::pair< DCR_MeshSolidShape2::Patch::bdt_node_sip_range,
                               DCR_MeshSolidShape2::Patch::bdt_node_sip_range > stack_entry_type;
            std::vector< stack_entry_type > stackBDTN;
            stackBDTN.push_back( stack_entry_type(DCR_MeshSolidShape2::Patch::bdt_node_sip_range(0,pd1.m_NumSegments),
                                                  DCR_MeshSolidShape2::Patch::bdt_node_sip_range(0,pd2.m_NumSegments) ) );
            while( !stackBDTN.empty() )
            {
                // Pop PS subarrays
                stack_entry_type se( stackBDTN.back() );
                stackBDTN.pop_back();
                const DCR_MeshSolidShape2::Patch::bdt_node_sip_range node_sr1( se.first );
                const DCR_MeshSolidShape2::Patch::bdt_node_sip_range node_sr2( se.second );
                const uint32 length_sr1( node_sr1.second - node_sr1.first );
                const uint32 length_sr2( node_sr2.second - node_sr2.first );
                if( length_sr1 * length_sr2 <= p_context->m_DCR2DCR_E2E_BDT_MaxLeafTests ) //use bruteforce below a given number of pairwise tests
                {
                    // Test all pairwise segments in subtrees \todo Consider inverting loop if length_sr2 > length_sr1
                    for( unsigned int it_sin1=0; it_sin1<length_sr1; it_sin1++ )
                    {
                        uint32 sid1( pd1.m_FirstSID + node_sr1.first + it_sin1 );
                        Vec2 a1_2( vec_v1_2[ p_dcr1->m_vecS[sid1].GetVID(0) - ed1.m_FirstVID ] );
                        Vec2 b1_2( vec_v1_2[ p_dcr1->m_vecS[sid1].GetVID(1) - ed1.m_FirstVID ] );
                        for( unsigned int it_sin2=0; it_sin2<length_sr2; it_sin2++ )
                        {
                            uint32 sid2( pd2.m_FirstSID + node_sr2.first + it_sin2 );
                            Vec2 a2_2( vec_v2_2[ p_dcr2->m_vecS[sid2].GetVID(0) - ed2.m_FirstVID ] );
                            Vec2 b2_2( vec_v2_2[ p_dcr2->m_vecS[sid2].GetVID(1) - ed2.m_FirstVID ] );
                            Real lambda1, lambda2;
                            if( Intersection_Segment2_Segment2( a1_2, b1_2, a2_2, b2_2, lambda1, lambda2 ) )
                            {
                                // save Crossing Feature Pair
                                vec_cfp.push_back( crossing_segment_pair_type( PointOnFeature( feature_id( eFT_Segment, sid1 ), Vec4(Real(1)-lambda1,lambda1,0,0) ), //POF1
                                                                               eid1, pid1, //eid1,pid1
                                                                               tr2_0*a1_2, tr2_0*b1_2, //a1_0,b1_0
                                                                               PointOnFeature( feature_id( eFT_Segment, sid2 ), Vec4(Real(1)-lambda2,lambda2,0,0) ), //POD2
                                                                               eid2, pid2, //eid2,pid2
                                                                               tr2_0*a2_2, tr2_0*b2_2 ) ); //a2_0,b2_0
                                num_crossings++;
                            }
                        }
                    }
                }
                else
                {
                    // Get node/segment data
                    const DCR_MeshSolidShape2::Segment& bdtn1( p_dcr1->m_vecS[pd1.m_FirstSID+node_sr1.first] );
                    const DCR_MeshSolidShape2::Segment& bdtn2( p_dcr2->m_vecS[pd2.m_FirstSID+node_sr2.first] );
                    const bv::BDOP3 bdop1_1( bdtn1.BDTN_BDOPq() );
                    const bv::BDOP3 bdop2_2( bdtn2.BDTN_BDOPq() );

                    // Get cartesian BVs from BDOPs
#ifdef __ENABLE_INTERSECTION_BDT_BDOP3_PRESORTED
                    bv::AABB2 bv1_2;
                    bv::AABB2 bv2_2;
                    if( p_context->m_DCR2DCR_E2E_BDT_Presort ) //\todo Use constructor instead of assignment when this param disappears
                    {
                        bv1_2 = bv::Compute_AABB2_From_BDOP3_PreSorted( bdop1_1, vec_node_pos1_2, vec_sorted_nie1 );
                        bv2_2 = bv::Compute_AABB2_From_BDOP3_PreSorted( bdop2_2, vec_node_pos2_2, vec_sorted_nie2 );
                    }
                    else
                    {
                        bv1_2 = bv::Compute_AABB2_From_BDOP3( bdop1_1, vec_node_pos1_2 );
                        bv2_2 = bv::Compute_AABB2_From_BDOP3( bdop2_2, vec_node_pos2_2 );
                    }
#else
                    bv::AABB2 bv1_2( bv::Compute_AABB2_From_BDOP3( bdop1_1, vec_node_pos1_2 ) );
                    bv::AABB2 bv2_2( bv::Compute_AABB2_From_BDOP3( bdop2_2, vec_node_pos2_2 ) );
#endif
                    // Test cartesian BV(BDOP)
                    if( bv::TestOverlap( bv1_2, bv2_2 ) )
                    {
                        // Recurse \todo I'M SURE this logic can be simplified taking into account the previous length_sr1*length_sr2<maxleaftests branch (MAYBE req maxleaftests>0, which is REASONABLE)
                        if( length_sr2 == 1 //\todo Compare BDOP sizes instead of length_sr to implement DescendLarger(n1,n2) strategy in RTCD pg.257 (as in GBVH)
                            || (length_sr1 > 1 && length_sr1 > length_sr2)  )
                        {
#ifdef __ENABLE_INTERSECTION_BDT_SPLIT3 //\todo I THINK the split logic can be simplified
                            // Recourse BDTN1.L/R + BDTN.S[0]
                            stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr1.first, node_sr1.first+1 ),
                                                                   node_sr2 ) );
                            int remaining_sr1( length_sr1 - 1 );
                            if( remaining_sr1 > 1 )
                            {

                                stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr1.first+1, node_sr1.first+1+remaining_sr1/2 ),
                                                                       node_sr2 ) );
                                stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr1.first+1+remaining_sr1/2, node_sr1.second ),
                                                                       node_sr2 ) );
                            }
                            else if( remaining_sr1 == 1 )
                                stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr1.first+1, node_sr1.first+2 ),
                                                                       node_sr2 ) );
#else
                            // Test segment1 vs all BDTN2.S
                            uint32 sid1( pd1.m_FirstSID + node_sr1.first );
                            Vec2 a1_2( vec_v1_2[ p_dcr1->m_vecS[sid1].GetVID(0) - ed1.m_FirstVID ] );
                            Vec2 b1_2( vec_v1_2[ p_dcr1->m_vecS[sid1].GetVID(1) - ed1.m_FirstVID ] );
                            for( unsigned int it_sin2=0; it_sin2<length_sr2; it_sin2++ )
                            {
                                uint32 sid2( pd2.m_FirstSID + node_sr2.first + it_sin2 );
                                Vec2 a2_2( vec_v2_2[ p_dcr2->m_vecS[sid2].GetVID(0) - ed2.m_FirstVID ] );
                                Vec2 b2_2( vec_v2_2[ p_dcr2->m_vecS[sid2].GetVID(1) - ed2.m_FirstVID ] );
                                Real lambda1, lambda2;
                                if( Intersection_Segment2_Segment2( a1_2, b1_2, a2_2, b2_2, lambda1, lambda2 ) )
                                {
                                    // save Crossing Feature Pair
                                    vec_cfp.push_back( crossing_segment_pair_type( PointOnFeature( feature_id( eFT_Segment, sid1 ), Vec4(Real(1)-lambda1,lambda1,0,0) ), //POF1
                                                                                   eid1, pid1, //eid1,pid1
                                                                                   tr2_0*a1_2, tr2_0*b1_2, //a1_0,b1_0
                                                                                   PointOnFeature( feature_id( eFT_Segment, sid2 ), Vec4(Real(1)-lambda2,lambda2,0,0) ), //POD2
                                                                                   eid2, pid2, //eid2,pid2
                                                                                   tr2_0*a2_2, tr2_0*b2_2 ) ); //a2_0,b2_0
                                    num_crossings++;
                                }
                            }
                            // Recourse BDTN1.L/R
                            int remaining_sr1( length_sr1 - 1 );
                            if( remaining_sr1 > 1 )
                            {
                                stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr1.first+1, node_sr1.first+1+remaining_sr1/2 ),
                                                                       node_sr2 ) );
                                stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr1.first+1+remaining_sr1/2, node_sr1.second ),
                                                                       node_sr2 ) );
                            }
                            else if( remaining_sr1 == 1 )
                                stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr1.first+1, node_sr1.first+2 ),
                                                                       node_sr2 ) );
#endif
                        }
                        else
                        {
#ifdef __ENABLE_INTERSECTION_BDT_SPLIT3
                            // Recourse BDTN2.L/R + BDTN2.S[0]
                            stackBDTN.push_back( stack_entry_type( node_sr1,
                                                                   DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr2.first, node_sr2.first+1 ) ) );
                            int remaining_sr2( length_sr2 - 1 );
                            if( remaining_sr2 > 1 )
                            {
                                stackBDTN.push_back( stack_entry_type( node_sr1,
                                                                       DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr2.first+1, node_sr2.first+1+remaining_sr2/2 ) ) );
                                stackBDTN.push_back( stack_entry_type( node_sr1,
                                                                       DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr2.first+1+remaining_sr2/2, node_sr2.second ) ) );
                            }
                            else if( remaining_sr2 == 1 )
                                stackBDTN.push_back( stack_entry_type( node_sr1,
                                                                       DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr2.first+1, node_sr2.first+2 ) ) );
#else
                            // Test segment1 vs all BDTN2.S
                            uint32 sid2( pd2.m_FirstSID + node_sr2.first );
                            Vec2 a2_2( vec_v2_2[ p_dcr2->m_vecS[sid2].GetVID(0) - ed2.m_FirstVID ] );
                            Vec2 b2_2( vec_v2_2[ p_dcr2->m_vecS[sid2].GetVID(1) - ed2.m_FirstVID ] );
                            for( unsigned int it_sin1=0; it_sin1<length_sr1; it_sin1++ )
                            {
                                uint32 sid1( pd1.m_FirstSID + node_sr1.first + it_sin1 );
                                Vec2 a1_2( vec_v1_2[ p_dcr1->m_vecS[sid1].GetVID(0) - ed1.m_FirstVID ] );
                                Vec2 b1_2( vec_v1_2[ p_dcr1->m_vecS[sid1].GetVID(1) - ed1.m_FirstVID ] );
                                Real lambda1, lambda2;
                                if( Intersection_Segment2_Segment2( a1_2, b1_2, a2_2, b2_2, lambda1, lambda2 ) )
                                {
                                    // save Crossing Feature Pair
                                    vec_cfp.push_back( crossing_segment_pair_type( PointOnFeature( feature_id( eFT_Segment, sid1 ), Vec4(Real(1)-lambda1,lambda1,0,0) ), //POF1
                                                                                   eid1, pid1, //eid1,pid1
                                                                                   tr2_0*a1_2, tr2_0*b1_2, //a1_0,b1_0
                                                                                   PointOnFeature( feature_id( eFT_Segment, sid2 ), Vec4(Real(1)-lambda2,lambda2,0,0) ), //POD2
                                                                                   eid2, pid2, //eid2,pid2
                                                                                   tr2_0*a2_2, tr2_0*b2_2 ) ); //a2_0,b2_0
                                    num_crossings++;
                                }
                            }
                            // Recourse BDTN2.L/R
                            int remaining_sr2( length_sr2 - 1 );
                            if( remaining_sr2 > 1 )
                            {
                                stackBDTN.push_back( stack_entry_type( node_sr1,
                                                                       DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr2.first+1, node_sr2.first+1+remaining_sr2/2 ) ) );
                                stackBDTN.push_back( stack_entry_type( node_sr1,
                                                                       DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr2.first+1+remaining_sr2/2, node_sr2.second ) ) );
                            }
                            else if( remaining_sr2 == 1 )
                                stackBDTN.push_back( stack_entry_type( node_sr1,
                                                                       DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr2.first+1, node_sr2.first+2 ) ) );
#endif
                        }
                    }
                }
            }
        }
    }
    return num_crossings;
}

unsigned int AddCFP_Intersection_DCR2_E_DCR2_E( const MeshSolidShape2* p_mesh1, const Vec2* vec_sdof1, const DCR_MeshSolidShape2* p_dcr1, uint32 eid1,
                                                const MeshSolidShape2* p_mesh2, const Vec2* vec_sdof2, const DCR_MeshSolidShape2* p_dcr2, uint32 eid2,
                                                const Transform2& tr1_2, const Transform2& tr2_0,
                                                std::vector< crossing_segment_pair_type >& vec_cfp,
                                                const Context* p_context = g_pDefaultContext )
{
    switch( p_context->m_DCR2DCR_CFP_Method )
    {
    case Context::eDCR2DCR_CFP_BruteForce: return AddCFP_Intersection_DCR2_E_DCR2_E_BruteForce(p_mesh1,vec_sdof1,p_dcr1,eid1,p_mesh2,vec_sdof2,p_dcr2,eid2,tr1_2,tr2_0,vec_cfp,p_context); break;
    case Context::eDCR2DCR_CFP_BDT: return AddCFP_Intersection_DCR2_E_DCR2_E_BDT(p_mesh1,vec_sdof1,p_dcr1,eid1,p_mesh2,vec_sdof2,p_dcr2,eid2,tr1_2,tr2_0,vec_cfp,p_context); break;
    default: return 0;
    }
}

#ifndef __ENABLE_MSS_DCR_PATCH_TOPOLOGY //\todo This methods MAY be useful anyway
std::pair<uint32,uint32> DCR_Find_EID_And_PID_From_SID( const DCR_MeshSolidShape2* p_dcr, uint32 sid )
{
    for( unsigned int it_e=0; it_e<p_dcr->m_NumElements; it_e++ )
    {
        const DCR_MeshSolidShape2::Element& ed( p_dcr->m_vecE[it_e] );
        for( unsigned int it_pie=0; it_pie<ed.m_NumPatches; it_pie++ )
        {
            const DCR_MeshSolidShape2::Patch& pd( p_dcr->m_vecP[ed.m_FirstPID+it_pie] );
            if( sid >= pd.m_FirstSID && sid < pd.m_FirstSID+pd.m_NumSegments )
                return std::make_pair(it_e,ed.m_FirstPID+it_pie);
        }
    }
    return std::make_pair(-1,-1);
}

bool DCR_Is_SID_Inside_EID( const DCR_MeshSolidShape2* p_dcr, uint32 sid, uint32 eid )
{
    const DCR_MeshSolidShape2::Element& ed( p_dcr->m_vecE[eid] );
    for( unsigned int it_pie=0; it_pie<ed.m_NumPatches; it_pie++ )
    {
        const DCR_MeshSolidShape2::Patch& pd( p_dcr->m_vecP[ed.m_FirstPID+it_pie] );
        if( sid >= pd.m_FirstSID && sid < pd.m_FirstSID+pd.m_NumSegments )
            return true;
    }
    return false;
}

bool DCR_Is_SID_Inside_PID( const DCR_MeshSolidShape2* p_dcr, uint32 sid, uint32 pid )
{
    const DCR_MeshSolidShape2::Patch& pd( p_dcr->m_vecP[pid] );
    if( sid >= pd.m_FirstSID && sid < pd.m_FirstSID+pd.m_NumSegments )
        return true;
    else
        return false;
}
#endif

struct intersection_boundary_type
{
    PointOnFeature m_BeginPOF;
    PointOnFeature m_EndPOF;

    //\note SID are NOT guaranteed to be globally consecutive, only per-DRC.E, we could store first,last per DCR.E, but we won't to keep it similar to 3D case where it's unfeasible
    //\todo Add per-DCR.S stuff if required
    std::vector<uint32> m_vecSID;

    // DCR.ElementPatch aggregate data
    struct EP
    {
        uint32 m_EID; uint32 m_PID;
        Vec2 m_AvgPos; Vec3 m_AvgBarycentricCoords; Vec2 m_AvgNormal; Real m_TotalLength; uint32 m_Valence;
        Mat3x3 m_invBs; Mat3x3 m_Bs_invBm;
        explicit inline EP( uint32 eid, uint32 pid )
        : m_EID(eid), m_PID(pid),
          m_AvgPos(0), m_AvgBarycentricCoords(0), m_AvgNormal(0), m_TotalLength(0), m_Valence(0)
            {}
    };
    std::vector<EP> m_vecEP;

    // IB aggregate data
    Vec2 m_AvgPos;
    Vec2 m_AvgNormal;
    Real m_TotalLength;
    explicit inline intersection_boundary_type( const PointOnFeature& begin_pof )
    : m_BeginPOF(begin_pof), m_AvgPos(0), m_AvgNormal(0), m_TotalLength(0) {}
};

/* Flood IB (on a given side object 0/1)
  - While unconsumed CFP
    - Find at an (entry) CFP_a where CCW order is INWARDS
    - Open a new IB_i
    - Add all DCR.F INWARDS until the next CFP_b is found //\todo If SID had global order, CFP_b could be found directly iterating over unconsumed CFP
      - Use CFP_b \lambda to sort CFP on the same feature CCW
      - Compute DCR.E( DCR.F ) and add if new
      - Can accumulate per DCR.E stuff
      - OPT: Can add whole DCR.E if they do not contain any CFP_b
    - Close IB
    - Consume both CFP_b
*/
bool FloodIB( std::vector<crossing_segment_pair_type>& vecSeeds,
              const MeshSolidShape2* p_mesh1, const DCR_MeshSolidShape2* p_dcr1, const Transform2& tr1, const Vec2* vec_sdof1,
              int side,
              const std::vector<bool>& vecIsCrossingP1,
              std::vector<intersection_boundary_type>& vecIB1,
              bool b_inwards, //otherwise,
              //TEMP, just for log/draw debug CP by now
              ContactData2& cd,
              bool b_log,
              bool b_draw )
{
    if( b_log ) GEO_LOG("FloodIB side=%d",side);
    unsigned int num_seeds( vecSeeds.size() );
    if( b_log ) GEO_LOG("#Seeds = %u", num_seeds);
    if( num_seeds % 2 != 0 )
    {
        // This fails due to geometric coincidences, silently return false
        GEO_LOG_ERROR( "Odd seed count, ignoring intersection, COINCIDENCES!?!?" );
        return false;
    }
    int side1( side );
    int side2( 1-side1 );
    const Vec2* default_sdof1( p_mesh1->GetVecDefaultSDOF() );
    for( unsigned int it_consumed=0; it_consumed<num_seeds; it_consumed+=2 )
    {
        // Find unconsumed CCW entry CFP \todo this REEVALUATES inwards/otwards seeds multiple times, it would be better to split into inwards/outwards lists in a single pre-pass!!
        unsigned int begin_seed_idx( num_seeds );
        for( unsigned int it_unconsumed=it_consumed; begin_seed_idx == num_seeds && it_unconsumed < num_seeds; it_unconsumed++ )
        {
            const crossing_segment_pair_type& cfp( vecSeeds[it_unconsumed] );
            // uint32 sid1( cfp.m_POF[side1].m_FeatureId.AsSegment() );
            // uint32 sid2( cfp.m_POF[side2].m_FeatureId.AsSegment() );
            Real d( mal::Dot( cfp.m_b_0[side1]-cfp.m_a_0[side2], mal::PerpendicularCW(cfp.m_b_0[side2]-cfp.m_a_0[side2]) ) );
            //\todo THIS TEST SEEMS REVERSED, inwards should be <0?, BUT WORKS... why??
            if( b_inwards && d > 0 )
                begin_seed_idx = it_unconsumed;
            if( !b_inwards && d < 0 )
                begin_seed_idx = it_unconsumed;
        }

        // GEO_ASSERT( begin_seed_idx < num_seeds );
        if( begin_seed_idx == num_seeds )
        {
            // This fails due to geometric coincidences, silently return false
            GEO_LOG_ERROR("INWARDS begin_seed_idx could not be found");
            return false;
        }

        // Consume BeginPOF
        std::swap( vecSeeds[it_consumed], vecSeeds[begin_seed_idx] );
        crossing_segment_pair_type& begin_cfp( vecSeeds[it_consumed] );
        // Open IB
        begin_cfp.m_IBID[side1] = vecIB1.size();
        vecIB1.push_back( intersection_boundary_type( begin_cfp.m_POF[side1] ) );
        intersection_boundary_type& ib( vecIB1.back() );

        // Add first segment/element
        unsigned int end_seed_idx( num_seeds );
        unsigned int sid( ib.m_BeginPOF.m_FeatureId.AsSegment() );
        unsigned int eid( begin_cfp.m_EID[side1] );
        unsigned int pid( begin_cfp.m_PID[side1] );
        ib.m_vecSID.push_back(sid);
        ib.m_vecEP.push_back( intersection_boundary_type::EP(eid,pid) );
        // Compute DCR.E barycentric deformation stuff
        intersection_boundary_type::EP& begin_ep( ib.m_vecEP.back() );
        {
            //IMPORTANT: Bs MAY be singular for degenerate
            //elements!!! In that case, the embedded segment MAY NOT BE
            //COLLAPSED, we could analyze the collapse configuration and
            //compute a valid "barycentric transform" that simply ignores
            //the collapsed vertices (eg: V-F collapsed tetrahedron would
            //ignore V barycentric weight and use only F vertices)
            uint32 vec_nid1[3];
            p_mesh1->P_VecVID( eid, vec_nid1, 3 ); //\todo NID could be stored in DCR::ED if not directly available in MSS
            Vec2 vec_pos_0[3] = { tr1*vec_sdof1[vec_nid1[0]], tr1*vec_sdof1[vec_nid1[1]], tr1*vec_sdof1[vec_nid1[2]] };
            Mat3x3 Bs( 1, 1, 1,
                       vec_pos_0[0].x(), vec_pos_0[1].x(), vec_pos_0[2].x(),
                       vec_pos_0[0].y(), vec_pos_0[1].y(), vec_pos_0[2].y() );
            Mat3x3 invBm( mal::Inverse( Mat3x3( 1, 1, 1,
                                                default_sdof1[vec_nid1[0]].x(), default_sdof1[vec_nid1[1]].x(), default_sdof1[vec_nid1[2]].x(),
                                                default_sdof1[vec_nid1[0]].y(), default_sdof1[vec_nid1[1]].y(), default_sdof1[vec_nid1[2]].y() ) ) ); //\todo invBm could be precomputed in DCR.ED!!
            begin_ep.m_invBs = mal::Inverse( Bs );
            begin_ep.m_Bs_invBm = Bs*invBm;
        }
        // Compute entry point from barycentrically deformed segment endpoints //\todo THIS is the CFP crossing point, could be saved into CFP on detection!
        Vec2 s0_0( mal::GRange<1,2>( begin_ep.m_Bs_invBm * mal::Concat(1,p_dcr1->m_vecV[p_dcr1->m_vecS[sid].GetVID(0)]) ) );
        Vec2 s1_0( mal::GRange<1,2>( begin_ep.m_Bs_invBm * mal::Concat(1,p_dcr1->m_vecV[p_dcr1->m_vecS[sid].GetVID(1)]) ) );
        Vec2 a1_0( ib.m_BeginPOF.m_BarycentricCoords[0] * s0_0 + ib.m_BeginPOF.m_BarycentricCoords[1] * s1_0 );

        /* Add CCW DCR.F to IB until first exit_cfp is found
           \note We'll track begin/end points for each added segment
           intervals. For completely enclosed segments, this
           will be the whole segment, but for entry/exit
           seed segments, we will only include the interior
           interval.
        */
        do
        {
            // Check if there's an unconsumed exit seed in the same SID
            if( true )// \todo IsCrossing(sid) OPTIMIZATION
            {
                /* Search ALL unconsumed looking for closest CCW potential end CFP.
                   \note The lambda parameter on a CCW segment
                   corresponds to m_BarycentricCoords[1], the
                   barycentric weight of the segment endpoint.
                */
                Real min_lambda( 1.01f );
                for( unsigned int it_unconsumed=it_consumed+1; it_unconsumed < num_seeds; it_unconsumed++ )
                {
                    const crossing_segment_pair_type& cfp( vecSeeds[it_unconsumed] );
                    if( cfp.m_POF[side1].m_FeatureId.AsSegment() == sid
                        && cfp.m_POF[side1].m_BarycentricCoords[1] < min_lambda )
                    {
                        min_lambda = cfp.m_POF[side1].m_BarycentricCoords[1];
                        if( min_lambda == ib.m_BeginPOF.m_BarycentricCoords[1] )
                            GEO_LOG_ERROR("Begin/End CFP are COINCIDENT!");
                        end_seed_idx = it_unconsumed;
                    }
                }
            }

            // If not found, advance CCW to neighbour SID
            if( end_seed_idx == num_seeds )
            {
                intersection_boundary_type::EP& ep( ib.m_vecEP.back() );
                // Compute overlap interval \note MAY have a1_0 = entry point (first point)
                Vec2 b1_0( mal::GRange<1,2>( ep.m_Bs_invBm * mal::Concat(1,p_dcr1->m_vecV[p_dcr1->m_vecS[sid].GetVID(1)]) ) );
                Real length( mal::Norm(a1_0-b1_0) );
                // Update IB aggregates
                ib.m_AvgPos += length*0.5*(a1_0+b1_0);
                ib.m_AvgNormal += mal::PerpendicularCW(a1_0-b1_0); //automatically length-weighted
                ib.m_TotalLength += length;
                // Update DCR.E aggregates
                ep.m_AvgPos += length*0.5*(a1_0+b1_0);
                ep.m_AvgNormal += mal::PerpendicularCW(a1_0-b1_0); //automatically length-weighted
                ep.m_TotalLength += length;
                // Advance SID
                sid = p_dcr1->m_vecS[sid].GetNSID(1);
                ib.m_vecSID.push_back(sid);
                //TEMP
                if( false )//b_draw )
                {
                    cd.AddCP( a1_0,b1_0, mal::SafeNormalized(b1_0-a1_0), mal::Norm(b1_0-a1_0) );
                }
                a1_0 = b1_0;
                // Advance DCR.EP if new
#ifdef __ENABLE_MSS_DCR_PATCH_TOPOLOGY
                if( sid == p_dcr1->m_vecP[p_dcr1->m_vecP[pid].m_vecNPID[1]].m_vecSID[0] ) // If sid == first_sid_in_next_pid, we're crossing a DCR.E/P border
#else
                if( !DCR_Is_SID_Inside_PID(p_dcr1,sid,pid) ) //\todo In 3D we won't have the cross-patch perimeter, thus we will detect DCR.S topology crossing to next DCR.P by inclusion (just a range-check, quite cheap)
#endif
                {
                    // Close previous DCR.E
                    if( ep.m_TotalLength > 0 )
                    {
                        ep.m_AvgPos = ep.m_AvgPos / ep.m_TotalLength;
                        ep.m_AvgNormal = mal::Normalized(ep.m_AvgNormal);
                        ep.m_AvgBarycentricCoords = ep.m_invBs * mal::Concat(1,ep.m_AvgPos);
                    }

#ifdef __ENABLE_MSS_DCR_PATCH_TOPOLOGY
                    // Advance pid/sid to the next CCW patch
                    pid = p_dcr1->m_vecP[pid].m_vecNPID[1];
                    GEO_ASSERT( sid == p_dcr1->m_vecP[pid].m_vecSID[0] ); //next SID from DCR.S topology must be the first SID from DCR.P topology
                    // sid = p_dcr1->m_vecP[pid].m_vecSID[0]; //BEGIN sid in next pid \note SID is already up-to-date
                    eid = p_dcr1->m_vecP[pid].m_EID; //save eid-from-pid for fast access
#else
                    // Open next DCR.E //\todo This SHOULD be accelerated by searching topologically from current pid instead of globally!!
                    std::pair<uint32,uint32> eid_pid( DCR_Find_EID_And_PID_From_SID(p_dcr1,sid) );
                    eid = eid_pid.first;
                    pid = eid_pid.second;
#endif

#ifdef __ENABLE_FLOODIB_WHOLE_NONCROSSING_PATCHES
                    // If DCR.P is not crossing, add it as a whole and advance to neighbour
                    while( !vecIsCrossingP1[pid] )
                    {
                        // GEO_LOG_WARNING( "Flooding whole patch (eid=%u,pid=%u,sid=%u)", eid, pid, sid );
                        ib.m_vecEP.push_back( intersection_boundary_type::EP(eid,pid) );
                        uint32 vec_nid1[3];
                        p_mesh1->P_VecVID( eid, vec_nid1, 3 ); //\todo NID could be stored in DCR::ED if not directly available in MSS
                        Vec2 vec_pos_0[3] = { tr1*vec_sdof1[vec_nid1[0]], tr1*vec_sdof1[vec_nid1[1]], tr1*vec_sdof1[vec_nid1[2]] };
                        Mat3x3 Bs( 1, 1, 1,
                                   vec_pos_0[0].x(), vec_pos_0[1].x(), vec_pos_0[2].x(),
                                   vec_pos_0[0].y(), vec_pos_0[1].y(), vec_pos_0[2].y() );
                        Mat3x3 invBm( mal::Inverse( Mat3x3( 1, 1, 1,
                                                            default_sdof1[vec_nid1[0]].x(), default_sdof1[vec_nid1[1]].x(), default_sdof1[vec_nid1[2]].x(),
                                                            default_sdof1[vec_nid1[0]].y(), default_sdof1[vec_nid1[1]].y(), default_sdof1[vec_nid1[2]].y() ) ) ); //\todo invBm could be precomputed in DCR.ED!!
                        intersection_boundary_type::EP& new_ep( ib.m_vecEP.back() );
                        new_ep.m_invBs = mal::Inverse( Bs );
                        new_ep.m_Bs_invBm = Bs*invBm;
                        // Add all aggregate WHOLE-patch stuff \todo TRY TO USE PRECOMPUTED undeformed data instead of recomputing it all...
                        const DCR_MeshSolidShape2::Patch& pd( p_dcr1->m_vecP[pid] );
#ifndef __ENABLE_MSS_DCR_PATCH_TOPOLOGY
                        GEO_ASSERT( pd.m_FirstSID == sid ); //Check that we're on the correct PID/SID
#endif
                        Vec2 new_ep_AvgNormal_Pose(0,0);
                        Vec2 new_ep_AvgPos_Pose(0,0);
                        Vec3 new_ep_AvgPos_Pose_b(0,0,0);
                        Real new_ep_TotalLength_Pose(0);
                        for( unsigned int it_sip=0; it_sip<pd.m_NumSegments; it_sip++ )
                        {
                            //\todo ib.m_vecSID.push_back(pd.m_FirstSID+it_sip); we should add this DCR.F, but we DO NOT USE THEM by now... we'll see
                            Vec2 s0_0( mal::GRange<1,2>( new_ep.m_Bs_invBm * mal::Concat(1,p_dcr1->m_vecV[p_dcr1->m_vecS[pd.m_FirstSID+it_sip].GetVID(0)]) ) );
                            Vec2 s1_0( mal::GRange<1,2>( new_ep.m_Bs_invBm * mal::Concat(1,p_dcr1->m_vecV[p_dcr1->m_vecS[pd.m_FirstSID+it_sip].GetVID(1)]) ) ); //\todo Swap s0/s1 to avoid recomputation!
                            Real length( mal::Norm(s0_0-s1_0) );
                            // Update IB aggregates
                            ib.m_AvgPos += length*0.5*(s0_0+s1_0); //\todo Consider different length-weighting to allow precomputation of whole-patch aggregates
                            ib.m_AvgNormal += mal::PerpendicularCW(s0_0-s1_0); //automatically length-weighted
                            ib.m_TotalLength += length;
                            // Update DCR.E aggregates
                            new_ep.m_AvgPos += length*0.5*(s0_0+s1_0);
                            new_ep.m_AvgNormal += mal::PerpendicularCW(s0_0-s1_0); //automatically length-weighted
                            new_ep.m_TotalLength += length;
#ifndef __ENABLE_MSS_DCR_PATCH_TOPOLOGY
                            a1_0 = s1_0; //IMPORTANT, leave current point at last used vertex!!
#endif
                            Vec2 s0_Pose( p_dcr1->m_vecV[p_dcr1->m_vecS[pd.m_FirstSID+it_sip].GetVID(0)] );
                            Vec2 s1_Pose( p_dcr1->m_vecV[p_dcr1->m_vecS[pd.m_FirstSID+it_sip].GetVID(1)] );
                            Vec2 diff_Pose( s0_Pose - s1_Pose );
                            Real length_Pose( mal::Norm(diff_Pose) );
                            // GEO_LOG( "d^0 error = %f", mal::Norm(diff_Pose-(s0_0-s1_0)) );
                            new_ep_AvgNormal_Pose += mal::PerpendicularCW( diff_Pose ); //==K*v
                            new_ep_AvgPos_Pose += length_Pose * 0.5 * (s0_Pose+s1_Pose);
                            Vec3 s0_Pose_b( invBm * mal::Concat(1,s0_Pose) );
                            Vec3 s1_Pose_b( invBm * mal::Concat(1,s1_Pose) );
                            new_ep_AvgPos_Pose_b += length_Pose * 0.5 * (s0_Pose_b+s1_Pose_b);
                            new_ep_TotalLength_Pose += length_Pose;
                        }
                        Vec2 new_ep_AvgPos_Global( mal::GRange<1,2>( new_ep.m_Bs_invBm * mal::Concat(1,new_ep_AvgPos_Pose) ) ); //\todo ALWAYS INCORRECT
                        Vec2 new_ep_AvgPos_Global_b( mal::GRange<1,2>( mal::Det(new_ep.m_Bs_invBm) * Bs * new_ep_AvgPos_Pose_b ) ); //\todo INCORRECT when actual deformation exists, OK for simple translations
                        /* todo If there is no way to transform
                           aggregated AvgPos, find Patch CoM instead
                           (using interior volume?) and transform it
                           barycentrically...

                           In any case, we are free to CHANGE the
                           definition of the "aggregated patch
                           position" so that it can be cheaply
                           transformed, as long as the modified
                           definition is consistent with "aggregation
                           of transformed geometry" for partially
                           overlapping patches.
                        */

                        //\todo See proof that N' = B^-T * N in http://stackoverflow.com/questions/13654401/what-is-the-logic-behind-transforming-normals-with-the-transpose-of-the-inverse
                        // Vec2 new_ep_AvgNormal_Global( mal::GRange<1,2>( mal::Inverse(mal::Transposed(new_ep.m_Bs_invBm)) * mal::Concat(1,new_ep_AvgNormal_Pose) ) ); //\todo INCORRECT if actual deformation exists, OK for translations
                        /* IMPORTANT: Derived on paper, THIS is the
                           correct way, AvgNormal must use 0 in the
                           extra coord as it comes from a point
                           substraction!

                           N^0 = K*v^0
                           N^1 = K*(M_0^1*v^0)
                               = L*(K*V^0)
                               = L*N^0
                           ==>
                           L = K*M*K^-1

                        Mat3x3 K( 1, 0, 0,
                                  0, 0, 1,
                                  0, -1, 0 ); //PerpendicularCW operator
                        Mat3x3 invK( mal::Inverse(K) );
                        Vec2 new_ep_AvgNormal_Global_OK( mal::GRange<1,2>( K * new_ep.m_Bs_invBm * invK * mal::Concat(0,new_ep_AvgNormal_Pose) ) );
                        */
                        // In 3D, for a general transform M, Mv1 x Mv2 = det M * M^-T * (v1 x v2), so that L = det M * M^-T. Applying the same L in 2D seems to work too, but I've got no detailed derivation...
                        // Vec2 new_ep_AvgNormal_Global_OK( mal::GRange<1,2>( mal::Det(new_ep.m_Bs_invBm) * mal::Inverse(mal::Transposed(new_ep.m_Bs_invBm)) * mal::Concat(0,new_ep_AvgNormal_Pose) ) ) ;
                        Vec2 new_ep_AvgNormal_Global_OK( mal::GRange<1,2>( mal::Transposed(mal::Adjugate(new_ep.m_Bs_invBm)) * mal::Concat(0,new_ep_AvgNormal_Pose) ) ); //OPT: det(M) M^-T = adj(M)^T (transposed adjugate)

                        /*\todo We NEED TO Transform TotalLength
                          too... let's hope it's possible... \see
                          http://en.wikipedia.org/wiki/Metric_tensor
                          for clues...

                          IMPORTANT: FUCKING SHIT... THIS
                          WORKS?!?!?!?!... MOSTLY... SAME rel error as
                          AvgPos_b, which I guess was to be
                          expected... Det(M) == Det(Range<1,1,2,2>(M))
                          because 1-coord in M is identity-transformed

                          - This is "approximately correct" for
                            deformations "aligned with the
                            patch"... as we're transforming LENGTH as
                            if it were AREA (det M * A is exact)... If
                            we find NO WAY to transform TotalLength in
                            O(1), we could consider using TotalArea
                            instead. Right now TotalLength is ONLY
                            strictly necessary to compute AvgPos, but
                            this will change if we use Avg of vertex
                            positions, instead of the patch CoM. We
                            were planning to use TotalLength for IB
                            comparison/weak order, BUT as long as
                            TotalArea is "somewhat" consistent with
                            that weak ordering, WE DO NOT CARE about
                            the actual TotalLength... In 3D, the same
                            will go for IB.TotalVolume instead of
                            TotalArea... Notice too that if we
                            actually aggregated area2d/volume3d of
                            interior elements, this approximation
                            would be EXACT (but we'd need to compute
                            interior SIM.E volumes that are not DCR.E,
                            which we'd like to avoid)

                          => IMPORTANT: derived on paper/mail with
                          Toni for D^2, he tested it and it
                          works. There seems to be NO WAY to transform
                          aggregate D, though, but D^2 is just fine.
                        */
                        Real new_ep_TotalLength_Global( mal::Det(mal::GRange<1,1,2,2>(new_ep.m_Bs_invBm)) * new_ep_TotalLength_Pose );
                        /*
                        GEO_LOG( "Error\n\tAbs L^0  %f  P^0 g%f b%f N^0 %f\n\tRel L^0  %f  P^0 g%f b%f N^0 %f",
                                 mal::Abs(new_ep_TotalLength_Global - new_ep.m_TotalLength),
                                 mal::Norm(new_ep_AvgPos_Global - new_ep.m_AvgPos),
                                 mal::Norm(new_ep_AvgPos_Global_b - new_ep.m_AvgPos),
                                 // mal::Norm(new_ep_AvgNormal_Global - new_ep.m_AvgNormal),
                                 mal::Norm(new_ep_AvgNormal_Global_OK - new_ep.m_AvgNormal),
                                 mal::Abs(new_ep_TotalLength_Global - new_ep.m_TotalLength) / new_ep.m_TotalLength,
                                 mal::Norm(new_ep_AvgPos_Global - new_ep.m_AvgPos) / mal::Norm(new_ep.m_AvgPos),
                                 mal::Norm(new_ep_AvgPos_Global_b - new_ep.m_AvgPos) / mal::Norm(new_ep.m_AvgPos),
                                 // mal::Norm(new_ep_AvgNormal_Global - new_ep.m_AvgNormal) / mal::Norm(new_ep.m_AvgNormal),
                                 mal::Norm(new_ep_AvgNormal_Global_OK - new_ep.m_AvgNormal) / mal::Norm(new_ep.m_AvgNormal) );
                        */
                        if( new_ep.m_TotalLength > 0 )
                        {
                            new_ep.m_AvgPos = new_ep.m_AvgPos / new_ep.m_TotalLength; //\todo THIS IS like a CoM, therefore, if m_AvgPos is precomputed, m_TotalLength is CONSTANT, as segment "mass" does NOT CHANGE when stretching, only DENSITY does, which we don't care about
                            new_ep.m_AvgNormal = mal::Normalized(new_ep.m_AvgNormal);
                            new_ep.m_AvgBarycentricCoords = new_ep.m_invBs * mal::Concat(1,new_ep.m_AvgPos);
                        }
#ifdef __ENABLE_MSS_DCR_PATCH_TOPOLOGY
                        // Advance sid to the next CCW patch
                        pid = pd.m_vecNPID[1];
                        sid = p_dcr1->m_vecP[pid].m_vecSID[0]; //BEGIN sid in next pid
                        eid = p_dcr1->m_vecP[pid].m_EID; //save eid-from-pid for fast access
                        //IMPORTANT: The last it_sip is NOT NECESSARILY the last CCW VID on this patch, thus we update it explicitly
                        a1_0 = mal::GRange<1,2>( new_ep.m_Bs_invBm * mal::Concat(1,p_dcr1->m_vecV[p_dcr1->m_vecS[sid].GetVID(0)]) );
#else
                        // Advance sid to the next CCW patch
                        sid = p_dcr1->m_vecS[ pd.m_FirstSID + pd.m_NumSegments - 1 ].GetNSID(1);
                        eid_pid = DCR_Find_EID_And_PID_From_SID(p_dcr1,sid);
                        eid = eid_pid.first;
                        pid = eid_pid.second;
#endif
                    }

                    // Open crossing patch (MAY or MAY NOT be the last, as crossing E may contain non-crossing P)
                    // Add DCR.E and compute its barycentric deformation stuff
                    {
                        ib.m_vecEP.push_back( intersection_boundary_type::EP(eid,pid) );
                        uint32 vec_nid1[3];
                        p_mesh1->P_VecVID( eid, vec_nid1, 3 ); //\todo NID could be stored in DCR::ED if not directly available in MSS
                        Vec2 vec_pos_0[3] = { tr1*vec_sdof1[vec_nid1[0]], tr1*vec_sdof1[vec_nid1[1]], tr1*vec_sdof1[vec_nid1[2]] };
                        Mat3x3 Bs( 1, 1, 1,
                                   vec_pos_0[0].x(), vec_pos_0[1].x(), vec_pos_0[2].x(),
                                   vec_pos_0[0].y(), vec_pos_0[1].y(), vec_pos_0[2].y() );
                        Mat3x3 invBm( mal::Inverse( Mat3x3( 1, 1, 1,
                                                            default_sdof1[vec_nid1[0]].x(), default_sdof1[vec_nid1[1]].x(), default_sdof1[vec_nid1[2]].x(),
                                                            default_sdof1[vec_nid1[0]].y(), default_sdof1[vec_nid1[1]].y(), default_sdof1[vec_nid1[2]].y() ) ) ); //\todo invBm could be precomputed in DCR.ED!!
                        ib.m_vecEP.back().m_invBs = mal::Inverse( Bs );
                        ib.m_vecEP.back().m_Bs_invBm = Bs*invBm;
                    }
#else
                    // Add DCR.E and compute its barycentric deformation stuff
                    {
                        ib.m_vecEP.push_back( intersection_boundary_type::EP(eid,-1) );
                        uint32 vec_nid1[3];
                        p_mesh1->P_VecVID( eid, vec_nid1, 3 ); //\todo NID could be stored in DCR::ED if not directly available in MSS
                        Vec2 vec_pos_0[3] = { tr1*vec_sdof1[vec_nid1[0]], tr1*vec_sdof1[vec_nid1[1]], tr1*vec_sdof1[vec_nid1[2]] };
                        Mat3x3 Bs( 1, 1, 1,
                                   vec_pos_0[0].x(), vec_pos_0[1].x(), vec_pos_0[2].x(),
                                   vec_pos_0[0].y(), vec_pos_0[1].y(), vec_pos_0[2].y() );
                        Mat3x3 invBm( mal::Inverse( Mat3x3( 1, 1, 1,
                                                            default_sdof1[vec_nid1[0]].x(), default_sdof1[vec_nid1[1]].x(), default_sdof1[vec_nid1[2]].x(),
                                                            default_sdof1[vec_nid1[0]].y(), default_sdof1[vec_nid1[1]].y(), default_sdof1[vec_nid1[2]].y() ) ) ); //\todo invBm could be precomputed in DCR.ED!!
                        ib.m_vecEP.back().m_invBs = mal::Inverse( Bs );
                        ib.m_vecEP.back().m_Bs_invBm = Bs*invBm;
                    }
#endif
                }
            }
            else // otherwise, add remaining segment interval to IB and close DCR.E at exit point
            {
                intersection_boundary_type::EP& ep( ib.m_vecEP.back() );
                // Compute overlap interval \note MAY have a1_0 = entry point (first point)
                const crossing_segment_pair_type& end_cfp( vecSeeds[end_seed_idx] );
                Vec2 s0_0( mal::GRange<1,2>( ep.m_Bs_invBm * mal::Concat(1,p_dcr1->m_vecV[p_dcr1->m_vecS[sid].GetVID(0)]) ) );
                Vec2 s1_0( mal::GRange<1,2>( ep.m_Bs_invBm * mal::Concat(1,p_dcr1->m_vecV[p_dcr1->m_vecS[sid].GetVID(1)]) ) );
                Vec2 b1_0( end_cfp.m_POF[side1].m_BarycentricCoords[0] * s0_0 + end_cfp.m_POF[side1].m_BarycentricCoords[1] * s1_0 );

                Real length( mal::Norm(a1_0-b1_0) );
                // Update IB aggregates
                ib.m_AvgPos += length*0.5*(a1_0+b1_0);
                ib.m_AvgNormal += mal::PerpendicularCW(a1_0-b1_0); //automatically length-weighted
                ib.m_TotalLength += length;
                // Update DCR.E aggregates
                ep.m_AvgPos += length*0.5*(a1_0+b1_0);
                ep.m_AvgNormal += mal::PerpendicularCW(a1_0-b1_0); //automatically length-weighted
                ep.m_TotalLength += length;
                // Close open DCR.E
                if( ep.m_TotalLength > 0 )
                {
                    ep.m_AvgPos = ep.m_AvgPos / ep.m_TotalLength;
                    ep.m_AvgNormal = mal::Normalized(ep.m_AvgNormal);
                    ep.m_AvgBarycentricCoords = ep.m_invBs * mal::Concat(1,ep.m_AvgPos);
                }
                //TEMP
                if( false )//b_draw )
                {
                    cd.AddCP( a1_0,b1_0, mal::SafeNormalized(b1_0-a1_0), mal::Norm(b1_0-a1_0) );
                }
            }
        } while( end_seed_idx == num_seeds ); //\todo Detect arriving back to entry_cfp, for safety?! \todo If IsCrossing(sid) were available, we would NOT NEED to check all CFP for each SID
        // Consume EndPOF
        std::swap( vecSeeds[it_consumed+1], vecSeeds[end_seed_idx] );
        crossing_segment_pair_type& end_cfp( vecSeeds[it_consumed+1] );
        // Close IB
        end_cfp.m_IBID[side1] = begin_cfp.m_IBID[side1];
        ib.m_EndPOF = end_cfp.m_POF[side1];
        if( ib.m_TotalLength > 0 )
        {
            ib.m_AvgPos = ib.m_AvgPos / ib.m_TotalLength;
            ib.m_AvgNormal = mal::Normalized(ib.m_AvgNormal); //\todo Comment-out to check divergence-theorem thing
        }
    }
    //TEMPORAL: debug
    for( unsigned int it_ib=0; it_ib<vecIB1.size(); it_ib++ )
    {
        if( b_log )
            GEO_LOG( "IB1[%u] with #SID = %d, #EP = %d",
                     it_ib, (int)vecIB1[it_ib].m_vecSID.size(), (int)vecIB1[it_ib].m_vecEP.size() );
        for( unsigned int it_ep=0; it_ep<vecIB1[it_ib].m_vecEP.size(); it_ep++ )
        {
            // GEO_LOG( "E[%u] = %u, P[%u], L=%f, P=(%f,%f), N=(%f,%f)",
            //          it_ep, vecIB1[it_ib].m_vecEP[it_ep].m_EID, vecIB1[it_ib].m_vecEP[it_ep].m_PID, vecIB1[it_ib].m_vecEP[it_ep].m_TotalLength,
            //          vecIB1[it_ib].m_vecEP[it_ep].m_AvgPos.x(), vecIB1[it_ib].m_vecEP[it_ep].m_AvgPos.y(),
            //          vecIB1[it_ib].m_vecEP[it_ep].m_AvgNormal.x(), vecIB1[it_ib].m_vecEP[it_ep].m_AvgNormal.y() );
            //Element aggregates
            if( false )//b_draw )
                cd.AddCP( vecIB1[it_ib].m_vecEP[it_ep].m_AvgPos,
                          vecIB1[it_ib].m_vecEP[it_ep].m_AvgPos + 0.1f*vecIB1[it_ib].m_vecEP[it_ep].m_AvgNormal,
                          vecIB1[it_ib].m_vecEP[it_ep].m_AvgNormal, 0.1f );
        }
        // IB pos/normal
        if( b_draw )
            cd.AddCP( vecIB1[it_ib].m_AvgPos, vecIB1[it_ib].m_AvgPos + 0.2f*vecIB1[it_ib].m_AvgNormal,
                      vecIB1[it_ib].m_AvgNormal, vecIB1[it_ib].m_TotalLength );//0.2f );
        //mal::Normalized(vecIB1[it_ib].m_AvgNormal), mal::Norm(vecIB1[it_ib].m_AvgNormal) ); //\todo Uncomment to check divergence-theorem thing
    }
    return true;
}

/* Find closest hit with any DCR.F in any pierced DCR.E present in the
   IB.

   \note We considered using up-to-date DCR.E BV from the BVH to cull
   raycasts, but we DO NOT KNOW which BVH.node contains a given DCR.E
   (eg: topdown build strategy reorders BVH entries), thus, we won't
   use BV or, if required, we'll recompute it from DCR.E geometry.

   + Consider using BSlab or BDOP instead of DCR.EP...
     - As we do Overlap_Segment2_Triangle2() in bcoords, the BSlab/BDOP
       can be used explicitly, there's no need to deform it!!
     - We can clip ray against all BDOP slabs and, if total interval is
       empty, there's no overlap
     => WORKS and culls >30% of RayCast_Segment2_SingleSided() calls

   \todo Consider using per-DCR.E static BVH (should be a specific,
   low-mem, BVH structure, ideally with implicit nodes, such as a ZBT
   (\see http://www.codercorner.com/ZeroByteBVH.pdf) or a simpler NMBVH
   (\see http://graphics.tu-bs.de/publications/Eisemann11NMH/)

   \todo For large IB the total number of IB.E and, therefore, of
   Overlap_Segment2_Triangle2() calls will be HUGE:
   - Consider some fast BVH/spatial hash to cull them
     globally...
     - If global parallel rays, a D-1 spatial hash in the directions
       perpendicular to the rays should be enough
   - Alternatively, consider using non-IB BVH but clipping the rays to
     the intersection of IB1,IB2 BV to reduce BVH access region.

   \note Ideally, RayCast_Segment2_XX should be single-sided but
   starting from a ray_pos INTERIOR to the object, thus, we'd reverse
   segment endpoints s0<->s1 and use RayCast_Segment2_SS. Using a
   DoubleSide query would return THE SAME point, but with reversed
   normal. However, ray_pos may NOT be interior for piercing IB.EP
   whose avgpos is OUTSIDE the other object => Thus, we'll call
   RayCast_DCR_E_DoubleSided_BruteForce(), which does not need to care
   about segment orientation (no need to reverse s0<->s1), but MAY
   REPORT OUTSIDE hits for initially exterior ray_pos, which MUST BE
   handled specifically (discarded?, we can check rh.m_Normal to
   detect it!)
*/
bool RaycastIB( const intersection_boundary_type& ib,
                const MeshSolidShape2* p_mesh, const DCR_MeshSolidShape2* p_dcr,
                const Vec2& p, const Vec2& n, Real max_length,
                RayHit2& rh )
{
    RayHit2 tmp_rh;
    rh.m_Interval.Set(max_length);
    //for each IB.E test ray and, if piercing, test all DCR.F in the DCR.E
    for( const auto& ib_ep : ib.m_vecEP )
    {
        // if( Overlap_Segment2_Triangle2( p, p+max_length*n, ib_ep.m_invBs ) ) //\todo Pass DCR.EP.BDOP and use it to cull ray
        if( Overlap_Segment2_Triangle2_BDOP( p, p+max_length*n, ib_ep.m_invBs, p_dcr->m_vecE[ib_ep.m_EID].m_BDOP ) )
        {
#ifdef __USE_RAYCAST_DCR_E
            //\todo WE NEED DCR.E vertices here, because the RayCast will be done in BCoords but must compute rh.m_Point and rh.m_Normal
            Mat3x3 Bs( mal::Inverse(ib_ep.m_invBs) ); //\todo SAVE THIS IN IB.EP
            //\todo SEGMENTS should be considered SINGLE-SIDED but REVERSED, as we're "inside" the object...
            if( RayCast_DCR_E_DoubleSided( p, n, Interval(0,rh.m_Interval.Min()), //clipped ray to optimize multiple rh candidate search
                                           p_dcr, ib_ep.m_EID,
                                           Bs, ib_ep.m_invBs, ib_ep.m_Bs_invBm,
                                           tmp_rh, 0 ) ) //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min()
                rh = tmp_rh; //\note This already contains rh data in global coords as well as feature_id=EID and barycentric coords wrt DCR.E
#else
            const DCR_MeshSolidShape2::Element& dcr_e( p_dcr->m_vecE[ib_ep.m_EID] );
            for( unsigned int it_pie=0; it_pie<dcr_e.m_NumPatches; it_pie++ )
            {
                const DCR_MeshSolidShape2::Patch& dcr_p( p_dcr->m_vecP[dcr_e.m_FirstPID+it_pie] );
                for( unsigned int it_sip=0; it_sip<dcr_p.m_NumSegments; it_sip++ )
                {
                    uint32 sid( dcr_p.m_FirstSID + it_sip );
                    Vec2 s0_0( mal::GRange<1,2>( ib_ep.m_Bs_invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecS[sid].GetVID(0) ]) ) );
                    Vec2 s1_0( mal::GRange<1,2>( ib_ep.m_Bs_invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecS[sid].GetVID(1) ]) ) );
                    //s0<->s1 swapped because p is assumed INSIDE p_dcr \todo This MAY NOT BE TRUE as p1 is an IB1.EP.m_AvgPos that MAY BE OUTSIDE DCR2
                    if( RayCast_Segment2_SingleSided( p, n, Interval(0,rh.m_Interval.Min()), s1_0, s0_0, tmp_rh, 0 ) ) //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min()
                    {
                        rh = tmp_rh;
                        rh.m_FeatureId = feature_id( eFT_Triangle, ib_ep.m_EID );
                        Vec3 b( ib_ep.m_invBs * mal::Concat(1,rh.m_Point) );
                        rh.m_Extra_BarycentricCoords = mal::Concat(b,0);
                    }
                }
            }
#endif
        }
    }
    return rh.m_Interval.Min() < max_length;
}

// TODO: Add RaycastDCR() that uses BVH.RayCast() internally, try to make it orthogonal to RC_Method
bool RaycastDCR( const intersection_boundary_type& ib,
                 const MeshSolidShape2* p_mesh, const Transform2& tr, const Vec2* vec_sdof,
                 const DCR_MeshSolidShape2* p_dcr, const BVH_MeshSolidShape2* p_bvh,
                 const Vec2& p, const Vec2& n, Real max_length,
                 RayHit2& rh )
{
    RayHit2 tmp_rh;
    rh.m_Interval.Set(max_length);
    const Vec2* default_sdof( p_mesh->GetVecDefaultSDOF() );
    std::vector< BVH_MeshSolidShape2::entry_index_type > vecOverlaps;
    if( p_bvh->RayCast( p, n, Interval(0,max_length), vecOverlaps ) )
    {
        for( auto eid : vecOverlaps )
        {
            if( eid < p_dcr->m_NumElements && p_dcr->m_vecE[eid].m_NumPatches > 0 ) //ignore non-DCR elements
            {
                //TEMP: invBs and Bs_invBm are EXPENSIVE to compute, but ARE NOT AVAILABLE for non-IB elements...
                uint32 vec_nid[3];
                p_mesh->P_VecVID( eid, vec_nid, 3 ); //\todo NID could be stored in DCR::ED if not directly available in MSS
                Vec2 vec_pos_0[3] = { tr*vec_sdof[vec_nid[0]], tr*vec_sdof[vec_nid[1]], tr*vec_sdof[vec_nid[2]] };
                Mat3x3 Bs( 1, 1, 1,
                           vec_pos_0[0].x(), vec_pos_0[1].x(), vec_pos_0[2].x(),
                           vec_pos_0[0].y(), vec_pos_0[1].y(), vec_pos_0[2].y() );
                Mat3x3 invBm( mal::Inverse( Mat3x3( 1, 1, 1,
                                                    default_sdof[vec_nid[0]].x(), default_sdof[vec_nid[1]].x(), default_sdof[vec_nid[2]].x(),
                                                    default_sdof[vec_nid[0]].y(), default_sdof[vec_nid[1]].y(), default_sdof[vec_nid[2]].y() ) ) ); //\todo invBm could be precomputed in DCR.ED!!
                Mat3x3 invBs( mal::Inverse( Bs ) );
                Mat3x3 Bs_invBm( Bs*invBm );
                if( Overlap_Segment2_Triangle2_BDOP( p, p+max_length*n, invBs, p_dcr->m_vecE[eid].m_BDOP ) )
                {
#ifdef __USE_RAYCAST_DCR_E
                    if( RayCast_DCR_E_DoubleSided( p, n, Interval(0,rh.m_Interval.Min()), //clipped ray to optimize multiple rh candidate search
                                                   p_dcr, eid,
                                                   Bs, invBs, Bs_invBm,
                                                   tmp_rh, 0 ) ) //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min()
                        rh = tmp_rh; //\note This already contains rh data in global coords as well as feature_id=EID and barycentric coords wrt DCR.E
#else
                    const DCR_MeshSolidShape2::Element& dcr_e( p_dcr->m_vecE[eid] );
                    for( unsigned int it_pie=0; it_pie<dcr_e.m_NumPatches; it_pie++ )
                    {
                        const DCR_MeshSolidShape2::Patch& dcr_p( p_dcr->m_vecP[dcr_e.m_FirstPID+it_pie] );
                        for( unsigned int it_sip=0; it_sip<dcr_p.m_NumSegments; it_sip++ )
                        {
                            uint32 sid( dcr_p.m_FirstSID + it_sip );
                            Vec2 s0_0( mal::GRange<1,2>( Bs_invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecS[sid].GetVID(0) ]) ) );
                            Vec2 s1_0( mal::GRange<1,2>( Bs_invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecS[sid].GetVID(1) ]) ) );
                            //s0<->s1 swapped because p is assumed INSIDE p_dcr \todo This MAY NOT BE TRUE as p1 is an IB1.EP.m_AvgPos that MAY BE OUTSIDE DCR2
                            if( RayCast_Segment2_SingleSided( p, n, Interval(0,rh.m_Interval.Min()), s1_0, s0_0, tmp_rh, 0 ) ) //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min()
                            {
                                rh = tmp_rh;
                                rh.m_FeatureId = feature_id( eFT_Triangle, eid );
                                Vec3 b( invBs * mal::Concat(1,rh.m_Point) );
                                rh.m_Extra_BarycentricCoords = mal::Concat(b,0);
                            }
                        }
                    }
#endif
                }
            }
        }
    }
    return rh.m_Interval.Min() < max_length;
}


//----- Heuristic "affinity" between boundariels (Larger ==> Better)
// typedef Real (*HeuristicFunc)( const Vec2 &p1, const Vec2 &n1, Real r1,
//                                const Vec2 &p2, const Vec2 &n2, Real r2 );

Real H_DistSq( const Vec2 &p1, const Vec2 &n1, Real r1,
               const Vec2 &p2, const Vec2 &n2, Real r2 )
{
    Vec2 v( p2 - p1 );
    Real norm_sq_v( mal::NormSq(v) );
    if( norm_sq_v < 1e-6 ) norm_sq_v = 1e-6;
    return mal::Rcp(norm_sq_v);
}

//\todo DOES NOT account for backfacing due to Sq(), useless...
inline Real H_DistSq_AngleSq( const Vec2& p1, const Vec2& n1, Real r1,
                              const Vec2& p2, const Vec2& n2, Real r2 )
{
    Vec2 v( p2 - p1 );
    Real larger_is_worse( mal::Sq( mal::Dot( v, mal::PerpendicularCW(n1) ) )
                          + mal::Sq( mal::Dot( v, mal::PerpendicularCW(n2) ) ) );
    if( larger_is_worse < 1e-6 ) larger_is_worse = 1e-6;
    return mal::Rcp(larger_is_worse);
}

/* Radiosity form factor (larger factor => larger energy transfer)
   \note Pi factor in divisor omitted
*/
inline Real H_RadiosityFormFactor( const Vec2& p1, const Vec2& n1, Real r1,
                                   const Vec2& p2, const Vec2& n2, Real r2 )
{
    Vec2 v( p2 - p1 );
    Real norm_sq_v( mal::NormSq(v) );
    if( norm_sq_v < 1e-6 ) return 0;
    Real form_factor( -mal::Dot( v, n1 ) * mal::Dot( v, n2 ) / mal::Sq(norm_sq_v) ); //\in [-1/norm_sq_r..1/norm_sq_r]
    return form_factor;
}

/* Affinity decreases:
   - Quadratically on distance |v|=|p2-p1|
   - On angle n1,n2
   - On angle n1,v
   - On angle n2,v
*/
inline Real H_Affinity( const Vec2& p1, const Vec2& n1, Real r1,
                        const Vec2& p2, const Vec2& n2, Real r2 )
{
    Vec2 v( p2 - p1 );
    Real norm_sq_v( mal::NormSq(v) ); //larger==>worse
    if( norm_sq_v < 1e-6 ) norm_sq_v = 1e-6;
    Real angle_n1n2( mal::Max<Real>(0,-mal::Dot(n1,n2)) ); //larger==>better, <=0 if perpendicular or backfacing
    Real angle_n1v_r( mal::Max<Real>(0,-mal::Dot(v,n1)) ); //larger==>better, <=0 if perpendicular or backfacing
    Real angle_n2v_r( mal::Max<Real>(0,mal::Dot(v,n2)) ); //larger==>better, <=0 if perpendicular or backfacing
    Real affinity( angle_n1n2*angle_n1v_r*angle_n2v_r / mal::Sq(norm_sq_v) ); //\in [0..1/norm_sq_v]
    return affinity;
}

inline Real Heuristic( const Vec2& p1, const Vec2& n1, Real r1,
                       const Vec2& p2, const Vec2& n2, Real r2,
                       const Context* p_context = g_pDefaultContext )
{
    switch( p_context->m_DCR2DCR_H_Method )
    {
    case Context::eDCR2DCR_H_DistSq: return H_DistSq(p1,n1,r1,p2,n2,r2); break;
    case Context::eDCR2DCR_H_DistSq_AngleSq: return H_DistSq_AngleSq(p1,n1,r1,p2,n2,r2); break;
    case Context::eDCR2DCR_H_FormFactor: return H_RadiosityFormFactor(p1,n1,r1,p2,n2,r2); break;
    case Context::eDCR2DCR_H_Affinity: return H_Affinity(p1,n1,r1,p2,n2,r2); break;
    default: return 0; break;
    }
}


/*\todo Assume that the BVH are up-to-date, get them through p_cc
  or add param p_bvh_geometry for each object and retrieve BVH
  "topology" from p_shape->GetBVH()
*/
bool Contact_DCRMS2_DCRMS2_BVH_MSS_Lazy_DCR( const MeshSolidShape2* p_mesh1, const Transform2& mesh_tr1, const Real* vec_dof1,
                                             const MeshSolidShape2* p_mesh2, const Transform2& mesh_tr2, const Real* vec_dof2,
                                             ContactData2& cd, ContactCache2* p_cc,
                                             const Context* p_context )
{
    cd.Begin();

    //\todo THIS IS dangerous... Should have a generic SRV type with safe casts...
    const Vec2* default_sdof1( p_mesh1->GetVecDefaultSDOF() );
    const Vec2* default_sdof2( p_mesh2->GetVecDefaultSDOF() );
    const Vec2* vec_sdof1( reinterpret_cast<const Vec2*>(vec_dof1) );
    const Vec2* vec_sdof2( reinterpret_cast<const Vec2*>(vec_dof2) );

    const DCR_MeshSolidShape2* pDCR1( p_mesh1->GetDCR() );
    const DCR_MeshSolidShape2* pDCR2( p_mesh2->GetDCR() );
    if( 0 == pDCR1 || 0 == pDCR2 )
    {
        GEO_LOG_WARNING("Contact_DCRMS2_DCRMS2_BVH_MSS_Lazy_DCR() no DCR1 && DCR2, ignoring...");
        return false; //\todo Fallback to a different method?
    }

    /*\todo We need to receive the up-to-date BVH SOMEHOW! this is an option...
    class MeshSolidS2MeshSolid_SPC2: public ISpecificPairwiseCache2
    {
        const BVH_MeshSolidShape2* m_pBVH1;
        const BVH_MeshSolidShape2* m_pBVH2;
    };
    const BVH_MeshSolidShape2* pBVH1( p_cc->m_pBVH1 );
    const BVH_MeshSolidShape2* pBVH2( p_cc->m_pBVH2 );
    if( 0 == pBVH1 || 0 == pBVH2 ) return false; //\todo Fallback to a different method?
    */

    /*\todo CONSIDER this... p_mesh are CONST and their BVH SHOULD NOT
      change here, as it should only contain the TOPOLOGY. Each object
      should store the BVH GEOMETRY and use it to Refit (either
      internally or on collision detection test if not up to date)
    if( pBVH1->GetGeometryTimestamp()m_GeometryTimeStamp != object1->m_TimeStamp )
    {
      pBVH1->Refit( mesh_tr1, vec_dof1 );
      pBVH1->SetGeometryTimestamp( object1->m_TimeStamp );
    }
    */

    //TEMP: By now, we'll just Refit the BVH unconditionally if present in the MSS
    BVH_MeshSolidShape2* pBVH1( p_mesh1->GetBVH() );
    BVH_MeshSolidShape2* pBVH2( p_mesh2->GetBVH() );
    if( 0 == pBVH1 || 0 == pBVH2 )
    {
        if( 0 == pBVH1 )
        {
            pBVH1 = new BVH_MeshSolidShape2;
            uint32 num_entries( pDCR1->m_NumElements );
            pBVH1->Rebuild_TopDown( num_entries,
                                    boost::bind<void>( &GEBV_MeshSolidShape2_E<BVH_MeshSolidShape2::entry_index_type,BVH_MeshSolidShape2::bv_type>,
                                                       p_mesh1, mesh_tr1, vec_sdof1, //Transform2::Identity(), default_sdof1,
                                                       _1, _2 ) );
            GEO_LOG_WARNING("Contact_DCRMS2_DCRMS2_BVH_MSS_Lazy_DCR() no BVH1, creating on the fly, UGLY...");
            const_cast<MeshSolidShape2*>(p_mesh1)->SetBakedBVH_StrictlyNonshared_UglyHack( pBVH1 ); //ugly const cast...
        }
        if( 0 == pBVH2 )
        {
            pBVH2 = new BVH_MeshSolidShape2;
            uint32 num_entries( pDCR2->m_NumElements );
            pBVH2->Rebuild_TopDown( num_entries,
                                    boost::bind<void>( &GEBV_MeshSolidShape2_E<BVH_MeshSolidShape2::entry_index_type,BVH_MeshSolidShape2::bv_type>,
                                                       p_mesh2, mesh_tr2, vec_sdof2, //Transform2::Identity(), default_sdof2,
                                                       _1, _2 ) );
            GEO_LOG_WARNING("Contact_DCRMS2_DCRMS2_BVH_MSS_Lazy_DCR() no BVH2, creating on the fly, UGLY...");
            const_cast<MeshSolidShape2*>(p_mesh2)->SetBakedBVH_StrictlyNonshared_UglyHack( pBVH2 ); //ugly const cast...
        }
        // GEO_LOG_WARNING("Contact_DCRMS2_DCRMS2_BVH_MSS_Lazy_DCR() no BVH1 && BVH2, ignoring...");
        // return false; //\todo Fallback to a different method?
    }

    //--- Refit BVH
    /*\todo Instead of re-transforming nodes for each element
      they're involved in, pre-transform them all in a single pass
      into a frame-allocated cache and access them by index in the
      GEBV functor by index.

      Ex:
      util::ScopedAllocator scoped_allocator( p_context->m_ScratchPad, "Name..." );
      Vec2* vec_node_pos1_0 = scoped_allocator.NewArrayPOD<Vec2>(pDCR1->m_NumVertices);
      mal::GTransform_Point_Array( vec_sdof1, tr1, vec_node_pos1_0 ); //\todo easily simd-able
    */

    // Refit both BVH
    /*\todo Refit using DCR BDOP/BSlab... there should be a predefined
      method in MeshSolidS2MeshSolid or DCR_MeshSolidShape2 or a free
      function to compute the GEBV in this cases... with a
      separate version that accepts cached transformed nodes that
      should be bound to the GEBV functor passed to Refit()
      pBVH1->Refit( boost::bind<void>( GEBV_BV_E_transformed_nodes, *p_mesh1, pDCR1, vec_node_pos1_0 ), _1, _2) );
      pBVH2->Refit( boost::bind<void>( GEBV_BV_E_transformed_nodes, *p_mesh2, pDCR2, vec_node_pos2_0 ), _1, _2) );

      \todo Why stop here, if KDOP is used, we could transform AND
      project all nodes ONCE and store their K/2 axis projections in
      an array passed to the BSlab/BDOP GEBV methods to avoid
      re-projecting nodes per-element therein...
    */
    switch( p_context->m_BVH_Method )
    {
    case Context::eBVHM_BV_E:
        pBVH1->Refit( boost::bind<void>( &GEBV_MeshSolidShape2_E<BVH_MeshSolidShape2::entry_index_type,BVH_MeshSolidShape2::bv_type>,
                                         p_mesh1, mesh_tr1, vec_sdof1, _1, _2) );
        pBVH2->Refit( boost::bind<void>( &GEBV_MeshSolidShape2_E<BVH_MeshSolidShape2::entry_index_type,BVH_MeshSolidShape2::bv_type>,
                                         p_mesh2, mesh_tr2, vec_sdof2, _1, _2) );
        break;
    case Context::eBVHM_BV_BSlab:
        pBVH1->Refit( boost::bind<void>( &GEBV_MeshSolidShape2_BSlab<BVH_MeshSolidShape2::entry_index_type,BVH_MeshSolidShape2::bv_type>,
                                         p_mesh1, pDCR1, mesh_tr1, vec_sdof1, _1, _2) );
        pBVH2->Refit( boost::bind<void>( &GEBV_MeshSolidShape2_BSlab<BVH_MeshSolidShape2::entry_index_type,BVH_MeshSolidShape2::bv_type>,
                                         p_mesh2, pDCR2, mesh_tr2, vec_sdof2, _1, _2) );
        break;
    case Context::eBVHM_BV_BDOP:
        pBVH1->Refit( boost::bind<void>( &GEBV_MeshSolidShape2_BDOP<BVH_MeshSolidShape2::entry_index_type,BVH_MeshSolidShape2::bv_type>,
                                         p_mesh1, pDCR1, mesh_tr1, vec_sdof1, _1, _2) );
        pBVH2->Refit( boost::bind<void>( &GEBV_MeshSolidShape2_BDOP<BVH_MeshSolidShape2::entry_index_type,BVH_MeshSolidShape2::bv_type>,
                                         p_mesh2, pDCR2, mesh_tr2, vec_sdof2, _1, _2) );
        break;
    case Context::eBVHM_BV_NoRefit: break;
    default: break;
    }

    // Test overlap \todo BUILD INSIDE SCOPED ALLOCATOR!!
    std::vector< std::pair<BVH_MeshSolidShape2::entry_index_type,BVH_MeshSolidShape2::entry_index_type> > vecOverlaps;

    bool bOverlap(false);
    switch( p_context->m_TestBVH_Method )
    {
    case Context::eTBVHM_BV_Vs_BV: bOverlap = pBVH1->Test( *pBVH2, vecOverlaps ); break;
    case Context::eTBVHM_BDOP_Vs_GSlabVertices:
        bOverlap = pBVH1->Test( *pBVH2,
                                [ p_mesh1, p_mesh2, pDCR1, pDCR2, mesh_tr1, mesh_tr2, vec_sdof1, vec_sdof2 ]
                                ( BVH_MeshSolidShape2::entry_index_type eid1, BVH_MeshSolidShape2::entry_index_type eid2,
                                  const BVH_MeshSolidShape2::bv_type& bv1, const BVH_MeshSolidShape2::bv_type& bv2 )
                                {
                                    GEO_ASSERT( eid1 < pDCR1->m_NumElements && eid2 < pDCR2->m_NumElements );
                                    return TestOverlap(bv1,bv2) //early-out, could be merged into TestBSlabs to reuse partial computations
                                    //&& TestTriangles(eid1,eid2) //early-out, could be merged into TestBSlabs to reuse partial computations
                                    && TestOverlap_BSlabs( p_mesh1, pDCR1, mesh_tr1, vec_sdof1, eid1, bv1,
                                                           p_mesh2, pDCR2, mesh_tr2, vec_sdof2, eid2, bv2 );
                                },
                                vecOverlaps );
        break;
    default: break;
    }
    //\todo This if should test pBVH1->Test() directly if there were no different m_TestBVH_Method
    if( bOverlap )
    {
        //GEO_LOG_WARNING("BVH2BVH %d overlaps", (int)vecOverlaps.size() );
        Transform2 tr1_2( mal::Inverse(mesh_tr2) * mesh_tr1 );
        std::vector< crossing_segment_pair_type > vecCFP;
        unsigned int num_crossings(0);
        for( auto overlap : vecOverlaps )
        {
            /*IMPORTANT!! If ANY object BVH includes internal elements
              NOT IN the DCR, this can happen, and we MUST AVOID
              accessing the inexistent DCR.E
              \todo CONSIDER detecting internal/internal overlaps and
              using/reporting them IF NO boundary overlap is detected.
            */
            if( overlap.first >= pDCR1->m_NumElements || overlap.second >= pDCR2->m_NumElements ) continue;

            /* \todo Consider filling an overlapped-DCR.E cache
               with Bs, Bs_invBm, etc... to avoid recomputation during
               Flood() and Raycast()
            */
            num_crossings += AddCFP_Intersection_DCR2_E_DCR2_E( p_mesh1, vec_sdof1, pDCR1, overlap.first,
                                                                p_mesh2, vec_sdof2, pDCR2, overlap.second,
                                                                tr1_2, mesh_tr2,
                                                                vecCFP, p_context );
        }

        /* Contact Determination
           - 2D:
             - Alg:
               1) Compute all CFP = (F_1, F_2)
               2) For each object O_1, O_2
                  - While unconsumed CFP
                   - Find at an (entry) CFP_a where CCW order is INWARDS
                   - Open a new IB_i
                   - Add all DCR.F INWARDS until the next CFP_b is found //\todo If SID had global order, CFP_b could be found directly iterating over unconsumed CFP
                     - Use CFP_b \lambda to sort CFP on the same feature CCW
                     - Compute DCR.E( DCR.F ) and add if new
                     - Can accumulate per DCR.E stuff
                     - OPT: Can add whole DCR.E if they do not contain any CFP_b
                   - Close IB
                   - Consume both CFP_b
               3) Build IB-graph
                 - Add one node per IB_1 and IB_2 (bipartite)
                 - For each CFP
                   - Add edge between (IB_1,IB2) connected
               4) Correspondences
                 - Analyze IB-graph
                 - ...
        */

        // Flag all crossing DCR.P
        std::vector<bool> vecIsCrossingP1(pDCR1->m_NumPatches,false);
        std::vector<bool> vecIsCrossingP2(pDCR2->m_NumPatches,false);
        for( const auto& cfp : vecCFP )
        {
            vecIsCrossingP1[ cfp.m_PID[0] ] = true;
            vecIsCrossingP2[ cfp.m_PID[1] ] = true;
        }

        // 2) Flood IB \todo COULD BE DONE IN PARALLEL if vecCFP was duplicated and merged afterwards
        std::vector<intersection_boundary_type> vecIB1;
        std::vector<intersection_boundary_type> vecIB2;
        if( !FloodIB( vecCFP,
                      p_mesh1, pDCR1, mesh_tr1, vec_sdof1,
                      0, //side 0 => IB1
                      vecIsCrossingP1,
                      vecIB1,
                      true, //inwards
                      cd, false, false ) ) return false;
        if( !FloodIB( vecCFP,
                      p_mesh2, pDCR2, mesh_tr2, vec_sdof2,
                      1, //side 1 => IB2
                      vecIsCrossingP2,
                      vecIB2,
                      true, //inwards
                      cd, false, false ) ) return false;

        // 3) Build IB-graph
        unsigned int num_ib1( vecIB1.size() );
        unsigned int num_ib2( vecIB2.size() );
        std::vector<bool> adj_mat_IB1xIB2( num_ib1*num_ib2, false ); //dense IB1xIB2 adjacency
        /* \todo MAY be required if EXTERNAL topology is considered
        unsigned int num_ib( num_ib1 + num_ib2 );
        std::vector<bool> symm_adj_mat_IB1_IB1( num_ib1 * (num_ib-1) / 2, false ); //LowerDiagonal IB1xIB1 adjacency
        */
        for( const auto& cfp : vecCFP )
        {
            uint32 ibid1( cfp.m_IBID[0] );
            uint32 ibid2( cfp.m_IBID[1] );
            adj_mat_IB1xIB2[ ibid1*num_ib2 + ibid2 ] = true;
            // Draw Viz correspondences
            // {
            //     cd.AddCP( vecIB1[ibid1].m_AvgPos, vecIB2[ibid2].m_AvgPos,
            //               mal::SafeNormalized(vecIB2[ibid2].m_AvgPos - vecIB1[ibid1].m_AvgPos),
            //               mal::Norm(vecIB2[ibid2].m_AvgPos - vecIB1[ibid1].m_AvgPos) );
            //     //TEMP: MUST ADD 1CP <=> 1POF...
            //     cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,0), Vec4::Zero() ),
            //                PointOnFeature( feature_id(eFT_Triangle,0), Vec4::Zero() ) );
            // }
        }

        // Compute IR
        struct intersection_region_type
        {
            std::vector<uint32> m_vecIBID1;
            std::vector<uint32> m_vecIBID2;
            // optional/incremental External Boundaries, parallel to respective vecIBID
            std::vector<intersection_region_type> m_vecEB1;
            std::vector<intersection_region_type> m_vecEB2;
        };
        //\todo Find vecIR, each one with its bounding IB, optional EB, aggregate data, etc...
        std::vector<uint32> vecUnassignedIBID1;
        for( unsigned int it_ib1=0; it_ib1<vecIB1.size(); it_ib1++ )
            vecUnassignedIBID1.push_back(it_ib1);
        while( !vecUnassignedIBID1.empty() )
        {
            uint32 ibid1( vecUnassignedIBID1.back() );
            vecUnassignedIBID1.pop_back();
            // Try to assign IBID1 to existing IR using IB1xIB2 adj, assign all xIB2 if match
            // Create new IR and assign IB and its xIB2 otherwise
        }

        /* 4) Global Untangling (GU)
        Compute IB-graph Connected Components
        -
        For each IB-graph CC
          Find global "untangling" and compute correspondences according to it
            Use IB aggregate data and IB-graph topology to decide untangling strategy
            - a) CUT some IB-CC adjacencies
            - b) COMPUTE unique IB-CC untengling direction (wrt Obj2)
                 - ComputeExtSizeIB()
                   - several approximations/implementations possible
                   - ComputeExtSizeIB_Length() Similar to FloodIB outwards but only compute DCR.S lengths
                   - ComputeExtSizeIB_Area() Similar to FloodIB outwards but compute DCR.S ext-area contributions using tri-surface from common point.
                     - MUST also consider opposed IB on the other object to close ext-region

                 - IMPORTANT: THIS will require backface-only
                   RayCasts, as we'll match exterior IB on one object
                   to interior IB on the other and we do NOT want to
                   detect "entrance" hits
        */
        /*\todo Call UntangleIR_XXX_YYY( ) for each IR
        for( const auto& ir : vecIR )
        {
            switch( p_context->m_DCR2DCR_Untangle_Method )
            {
            case eDCR2DCR_Untangle_Correspondences: UntangleIR_Correspondences( ir, cd ); break;
            case eDCR2DCR_Untangle_Min_ExtSize: UntangleIR_Min_ExtSize( ir, cd ); break;
            case eDCR2DCR_Untangle_Avg_ExtSize: UntangleIR_Avg_ExtSize( ir, cd ); break;
                //\todo OTHER strategies possible
            default: break;
            }
        }
        */

        // todo The FOLLOWING untangling code should moved to UntangleIR_Correspondences()

        // 4) Correspondences
        struct element_correspondence
        {
            uint32 m_IBID1;
            uint32 m_IBID2;
            uint32 m_EIB1; //EIB1,EIB2 relative to IB1,IB2
            uint32 m_EIB2;
            Real m_Heuristic;

            inline element_correspondence() {}
            inline element_correspondence( uint32 ibid1, uint32 ibid2, uint32 eib1, uint32 eib2, Real heuristic )
            : m_IBID1(ibid1), m_IBID2(ibid2), m_EIB1(eib1), m_EIB2(eib2), m_Heuristic(heuristic) {}
        };

        /* 4.a) Local \todo after GU
           For each IB1
             ECP: Find best ECP2ECP correspondences between (IB1,IB2) \note (E1,E2) may repeat for DIFFERENT, but will contain different Feature sets (independent IB patches)
               FCP: Find best FCP2FCP correspondences \note (F1,F2) CANNOT repeat for repeated (E1,E2) because they must refer disjoint feature patchese
        */
        std::vector< element_correspondence > vecEC;
        for( unsigned int it_ib1=0; it_ib1<num_ib1; it_ib1++ )
        {
            intersection_boundary_type& ib1( vecIB1[it_ib1] );
            // Gather IB1-->IB2 adjacency list
            std::vector<uint32> vec_adj_ibid2;
            for( unsigned int it_ib2=0; it_ib2<num_ib2; it_ib2++ )
                if( adj_mat_IB1xIB2[ it_ib1*num_ib2 + it_ib2 ] )
                    vec_adj_ibid2.push_back(it_ib2);
            // For each E1, find best E2
            for( unsigned int eib1=0; eib1<ib1.m_vecEP.size(); eib1++ )
            {
                const intersection_boundary_type::EP& ep1( ib1.m_vecEP[eib1] );
                uint32 best_ibid2(0);
                uint32 best_eib2(0);
                Real best_h(mal::MinusInfinity<Real>());
                for( auto ibid2 : vec_adj_ibid2 )
                {
                    const intersection_boundary_type& ib2( vecIB2[ibid2] );
                    for( unsigned int eib2=0; eib2<ib2.m_vecEP.size(); eib2++ )
                    {
                        const intersection_boundary_type::EP& ep2( ib2.m_vecEP[eib2] );
                        //\todo align/depth heuristics MAY NEED TO be prioritized or weighted... MUST BE symmetric in E1,E2
                        Real h = Heuristic( ep1.m_AvgPos, ep1.m_AvgNormal, ep1.m_TotalLength,
                                            ep2.m_AvgPos, ep2.m_AvgNormal, ep2.m_TotalLength,
                                            p_context );
                        if( h > best_h )
                        {
                            best_h = h;
                            best_ibid2 = ibid2;
                            best_eib2 = eib2;
                        }
                        //\todo Find E1.F1,E2.F2 sub-matchings here, if required
                    }
                }
                // Add correspondence if significant
                if( best_h > 0 )
                {
                    //\todo Create E1,E2 pair, with all required data ("contact confidence", etc...)
                    vecEC.push_back( element_correspondence( it_ib1, best_ibid2, eib1, best_eib2, best_h ) );
                    // Increase global DCR.E valences \todo CONSIDER IB-LOCAL DCR.E valences instead...
                    ib1.m_vecEP[eib1].m_Valence++;
                    vecIB2[best_ibid2].m_vecEP[best_eib2].m_Valence++;
                }
            }
        }
        unsigned int num_primary_correspondences( vecEC.size() );
        // Process unmatched E2 (valence==0), which were NOT the best-match for any E1 in an adjacent IB pair
        for( unsigned int it_ib2=0; it_ib2<num_ib2; it_ib2++ )
        {
            intersection_boundary_type& ib2( vecIB2[it_ib2] );
            // Gather IB2-->IB1 adjacency list
            std::vector<uint32> vec_adj_ibid1;
            for( unsigned int it_ib1=0; it_ib1<num_ib1; it_ib1++ )
                if( adj_mat_IB1xIB2[ it_ib1*num_ib2 + it_ib2 ] )
                    vec_adj_ibid1.push_back(it_ib1);
            // For each E1, find best E2
            for( unsigned int eib2=0; eib2<ib2.m_vecEP.size(); eib2++ )
            {
                const intersection_boundary_type::EP& ep2( ib2.m_vecEP[eib2] );
                if( ep2.m_Valence == 0 )
                {
                    uint32 best_ibid1(0);
                    uint32 best_eib1(0);
                    Real best_h(mal::MinusInfinity<Real>());
                    for( auto ibid1 : vec_adj_ibid1 )
                    {
                        const intersection_boundary_type& ib1( vecIB1[ibid1] );
                        for( unsigned int eib1=0; eib1<ib1.m_vecEP.size(); eib1++ )
                        {
                            const intersection_boundary_type::EP& ep1( ib1.m_vecEP[eib1] );
                            //\todo align/depth heuristics MAY NEED TO be prioritized or weighted... MUST BE symmetric in E1,E2
                            Real h = Heuristic( ep1.m_AvgPos, ep1.m_AvgNormal, ep1.m_TotalLength,
                                                ep2.m_AvgPos, ep2.m_AvgNormal, ep2.m_TotalLength,
                                                p_context );
                            if( h > best_h )
                            {
                                best_h = h;
                                best_ibid1 = ibid1;
                                best_eib1 = eib1;
                            }
                            //\todo Find E1.F1,E2.F2 sub-matchings here, if required
                        }
                    }
                    // Add correspondence if significant
                    if( best_h > 0 )
                    {
                        //\todo Create E1,E2 pair, with all required data ("contact confidence", etc...)
                        vecEC.push_back( element_correspondence( best_ibid1, it_ib2, best_eib1, eib2, best_h ) );
                        // Increase global DCR.E valences \todo CONSIDER IB-LOCAL DCR.E valences instead...
                        ib2.m_vecEP[eib2].m_Valence++;
                        vecIB1[best_ibid1].m_vecEP[best_eib1].m_Valence++;
                    }
                }
            }
        }

        // if( vecEC.size() > num_primary_correspondences )
        // GEO_LOG( "#EC = %d (%d primary + %d secondary)",
        //          (int)vecEC.size(), (int)num_primary_correspondences, (int)(vecEC.size() - num_primary_correspondences) );

        /*\todo At this point, correspondences with VERY SMALL
          heuristics may exist, and denote bad correspondences without
          a better alternative. Consider normalizing the heuristic
          range to [0..1] for each disjoint contact manifold, and
          using the heuristic as a "confidence" parameter
        */

        // 5) Generate CP from EC data
        Real largest_extent1( bv::ComputeLargestExtent(pBVH1->GetRoot().m_Geometry.m_BV) );
        Real largest_extent2( bv::ComputeLargestExtent(pBVH2->GetRoot().m_Geometry.m_BV) ); //\todo Actually, we could use largest extent of the intersection of BV1,BV2
        for( const auto& ec : vecEC )
        {
            uint32 eid1( vecIB1[ec.m_IBID1].m_vecEP[ec.m_EIB1].m_EID );
            uint32 eid2( vecIB2[ec.m_IBID2].m_vecEP[ec.m_EIB2].m_EID );
            Vec2 p1_0( vecIB1[ec.m_IBID1].m_vecEP[ec.m_EIB1].m_AvgPos );
            Vec2 p2_0( vecIB2[ec.m_IBID2].m_vecEP[ec.m_EIB2].m_AvgPos );
            Vec2 n1_0( vecIB1[ec.m_IBID1].m_vecEP[ec.m_EIB1].m_AvgNormal );
            Vec2 n2_0( vecIB2[ec.m_IBID2].m_vecEP[ec.m_EIB2].m_AvgNormal );
            Vec3 p1_b( vecIB1[ec.m_IBID1].m_vecEP[ec.m_EIB1].m_AvgBarycentricCoords );
            Vec3 p2_b( vecIB2[ec.m_IBID2].m_vecEP[ec.m_EIB2].m_AvgBarycentricCoords );
            switch( p_context->m_DCR2DCR_IM_Method )
            {
            case Context::eDCR2DCR_IM_Basic: //a) Basic correspondences
                cd.AddCP( p1_0, p2_0,
                          mal::SafeNormalized(p2_0-p1_0),
                          mal::Norm(p2_0-p1_0) );
                // \note If added, POF and CP must be in parallel arrays!
                cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,eid1), mal::Concat(p1_b,0) ),
                           PointOnFeature( feature_id(eFT_Triangle,eid2), mal::Concat(p2_b,0) ) );
                break;
            case Context::eDCR2DCR_IM_Projection: //b) Projection: Project along EC normal, find other point, generate 2 CP with full info (barycentric coords wrt DCR.E, etc...)
                //TEMP => Projection complicates computing barycentric coords (must compute enclosing DCR.E and its bcoords of an arbitrary point in space)
                break;
            case Context::eDCR2DCR_IM_Raycast: //c) Raycast: Raycast along EC normal, find other point, generate 2 CP with full info (barycentric coords wrt DCR.E, etc...)
                {
                    //IMPORTANT: E1,E2 are GUARANTEED to be, at least,
                    //partially INSIDE DCR2,DCR1 respectively, HOWEVER,
                    //E1.m_AvgPos MAY NOT be inside DCR2 and, therefore,
                    //any raycast from it may FAIL TO HIT DCR2/IB2. We DO
                    //NOT WANT to handle non-interior E1.m_AvgPos
                    //specifically, so we'll use a normal alignment
                    //criterion to avoid reversed CP normals in this case.

                    //I) Global RC: ECP1,IB2.Raycast(ECP1,IB2.N) => This should work well for "mostly convex" intersection regions //\todo CAN BE OPTIMIZED due to parallel raycasts! //\todo DOES NOT REQUIRE CORRESPONDENCES
                    //II) Direct RC: ECP1,IB2.Raycast(ECP1,-n1) => This should work well for similar n1,n2 \todo DOES NOT REQUIRE CORRESPONDENCES, should be similar to "raycasted collision detection" paper, ugly
                    //III) Inverse RC: ECP1,IB2.Raycast(ECP1,n2) => This should also work well for very different n1,n2 (eg: orthogonal)
                    //IV) Average RC: avg(n1,n2) => This tries to combine the strengths of II) and III), but the results are not very good...

                    //1) IB.Raycast => SOME points will not have RC correspondence... use basic correspondence | closest point on IB this cases? RC ideally cheaper, only IB.E need to be tested
                    // => IB raycast using all involved DCR.E in the IB, regardless of their piercing or completely internal status (=>This allows to iterate over IB.E instead of IB.F, and use DCR.E BV to speedup tests
                    //2) DCR.Raycast => All points will have a RC correspondence, but RC more expensive (can use BVH, though)

                    // Select raycast directions I,II,III,IV
                    Vec2 ray_dir1(0,0);
                    Vec2 ray_dir2(0,0);
                    switch( p_context->m_DCR2DCR_RC_Method )
                    {
                    case Context::eDCR2DCR_RC_Global: //I
                        /* \todo THIS WILL ONLY WORK for 1 IB1 <==> 1
                         * IB2, as it assumes N1 = -N2... for complex
                         * regions, a SINGLE IB1.N should be selected
                         * and negated for its >1 correspondent IB2
                         * (or the other way around!) */
                        ray_dir1 = vecIB2[ec.m_IBID2].m_AvgNormal;
                        ray_dir2 = vecIB1[ec.m_IBID1].m_AvgNormal;
                        break;
                    case Context::eDCR2DCR_RC_Direct: //II
                        ray_dir1 = -n1_0;
                        ray_dir2 = -n2_0;
                        break;
                    case Context::eDCR2DCR_RC_Inverse: //III
                        ray_dir1 = n2_0;
                        ray_dir2 = n1_0;
                        break;
                    case Context::eDCR2DCR_RC_Average: //IV
                        ray_dir1 = mal::SafeNormalized(n2_0-n1_0);
                        ray_dir2 = mal::SafeNormalized(n1_0-n2_0);
                        break;
                    default: break;
                    }

                    // Perform raycast against IB/DCR, ignoring if non-inwards wrt normal
                    bool bHit1( false );
                    bool bHit2( false );
                    RayHit2 rh1,rh2;
                    if( p_context->m_DCR2DCR_RC_UseOnlyIB ) //1) IB Raycast
                    {
                        bHit1 = mal::Dot(ray_dir1,n1_0) < 0 && RaycastIB( vecIB2[ec.m_IBID2], p_mesh2, pDCR2, p1_0, ray_dir1, largest_extent2, rh1 );
                        bHit2 = mal::Dot(ray_dir2,n2_0) < 0 && RaycastIB( vecIB1[ec.m_IBID1], p_mesh1, pDCR1, p2_0, ray_dir2, largest_extent1, rh2 );
                    }
                    else //2) DCR Raycast
                    {
                        bHit1 = mal::Dot(ray_dir1,n1_0) < 0 && RaycastDCR( vecIB2[ec.m_IBID2], p_mesh2, mesh_tr2, vec_sdof2, pDCR2, pBVH2, p1_0, ray_dir1, largest_extent2, rh1 );
                        bHit2 = mal::Dot(ray_dir2,n2_0) < 0 && RaycastDCR( vecIB1[ec.m_IBID1], p_mesh1, mesh_tr1, vec_sdof1, pDCR1, pBVH1, p2_0, ray_dir2, largest_extent1, rh2 );
                    }

                    // Perform global raycast if required to choose closest
                    if( p_context->m_DCR2DCR_RC_UseClosestOrGlobal ) //TEMP: This is an experiment to solve the Concave Vs Convex problem (single-region box in a hole) and does NOT seem to be any better than Global normals
                    {
                        RayHit2 rh1_global, rh2_global;
                        Vec2 ray_dir1_global( vecIB2[ec.m_IBID2].m_AvgNormal );
                        Vec2 ray_dir2_global( vecIB1[ec.m_IBID1].m_AvgNormal );
                        bool bHit1_global = mal::Dot(ray_dir1_global,n1_0) < 0 && RaycastDCR( vecIB2[ec.m_IBID2], p_mesh2, mesh_tr2, vec_sdof2, pDCR2, pBVH2, p1_0, ray_dir1_global, largest_extent2, rh1_global );
                        bool bHit2_global = mal::Dot(ray_dir2_global,n2_0) < 0 && RaycastDCR( vecIB1[ec.m_IBID1], p_mesh1, mesh_tr1, vec_sdof1, pDCR1, pBVH1, p2_0, ray_dir2_global, largest_extent1, rh2_global );

                        // Choose closest hit1
                        if( bHit1 && bHit1_global )
                        {
                            Real dist( mal::Norm(rh1.m_Point-p1_0) );
                            Real dist_global( mal::Norm(rh1_global.m_Point-p1_0) );
                            if( dist < dist_global )
                            {
                                cd.AddCP( p1_0, rh1.m_Point,
                                          mal::SafeNormalized(rh1.m_Point-p1_0),
                                          dist );
                                // rh contains feature_id=EID with valid extra barycentric coords wrt the DCR.E... we do NOT CARE about the actual DCR.F hit
                                cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,eid1), mal::Concat(p1_b,0) ),
                                           PointOnFeature( feature_id(eFT_Triangle,rh1.m_FeatureId.AsTriangle()), rh1.m_Extra_BarycentricCoords ) );
                            }
                            else
                            {
                                cd.AddCP( p1_0, rh1_global.m_Point,
                                          mal::SafeNormalized(rh1_global.m_Point-p1_0),
                                          dist_global );
                                // rh contains feature_id=EID with valid extra barycentric coords wrt the DCR.E... we do NOT CARE about the actual DCR.F hit
                                cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,eid1), mal::Concat(p1_b,0) ),
                                           PointOnFeature( feature_id(eFT_Triangle,rh1_global.m_FeatureId.AsTriangle()), rh1_global.m_Extra_BarycentricCoords ) );
                            }
                        }
                        else if( bHit1 )
                        {
                            cd.AddCP( p1_0, rh1.m_Point,
                                      mal::SafeNormalized(rh1.m_Point-p1_0),
                                      mal::Norm(rh1.m_Point-p1_0) );
                            // rh contains feature_id=EID with valid extra barycentric coords wrt the DCR.E... we do NOT CARE about the actual DCR.F hit
                            cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,eid1), mal::Concat(p1_b,0) ),
                                       PointOnFeature( feature_id(eFT_Triangle,rh1.m_FeatureId.AsTriangle()), rh1.m_Extra_BarycentricCoords ) );
                        }
                        else if( bHit1_global )
                        {
                            cd.AddCP( p1_0, rh1_global.m_Point,
                                      mal::SafeNormalized(rh1_global.m_Point-p1_0),
                                      mal::Norm(rh1_global.m_Point-p1_0) );
                            // rh contains feature_id=EID with valid extra barycentric coords wrt the DCR.E... we do NOT CARE about the actual DCR.F hit
                            cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,eid1), mal::Concat(p1_b,0) ),
                                       PointOnFeature( feature_id(eFT_Triangle,rh1_global.m_FeatureId.AsTriangle()), rh1_global.m_Extra_BarycentricCoords ) );
                        }

                        // Choose closest hit2
                        if( bHit2 && bHit2_global )
                        {
                            Real dist( mal::Norm(p2_0-rh2.m_Point) );
                            Real dist_global( mal::Norm(p2_0-rh2_global.m_Point) );
                            if( dist < dist_global )
                            {
                                cd.AddCP( rh2.m_Point, p2_0,
                                          mal::SafeNormalized(p2_0-rh2.m_Point),
                                          dist );
                                // rh contains feature_id=EID with valid extra barycentric coords wrt the DCR.E... we do NOT CARE about the actual DCR.F hit
                                cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,rh2.m_FeatureId.AsTriangle()), rh2.m_Extra_BarycentricCoords ),
                                           PointOnFeature( feature_id(eFT_Triangle,eid2), mal::Concat(p2_b,0) ) );
                            }
                            else
                            {
                                cd.AddCP( rh2_global.m_Point, p2_0,
                                          mal::SafeNormalized(p2_0-rh2_global.m_Point),
                                          dist_global );
                                // rh contains feature_id=EID with valid extra barycentric coords wrt the DCR.E... we do NOT CARE about the actual DCR.F hit
                                cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,rh2_global.m_FeatureId.AsTriangle()), rh2_global.m_Extra_BarycentricCoords ),
                                           PointOnFeature( feature_id(eFT_Triangle,eid2), mal::Concat(p2_b,0) ) );
                            }
                        }
                        else if( bHit2 )
                        {
                            cd.AddCP( rh2.m_Point, p2_0,
                                      mal::SafeNormalized(p2_0-rh2.m_Point),
                                      mal::Norm(p2_0-rh2.m_Point) );
                            // rh contains feature_id=EID with valid extra barycentric coords wrt the DCR.E... we do NOT CARE about the actual DCR.F hit
                            cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,rh2.m_FeatureId.AsTriangle()), rh2.m_Extra_BarycentricCoords ),
                                       PointOnFeature( feature_id(eFT_Triangle,eid2), mal::Concat(p2_b,0) ) );
                        }
                        else if( bHit2_global )
                        {
                            cd.AddCP( rh2_global.m_Point, p2_0,
                                      mal::SafeNormalized(p2_0-rh2_global.m_Point),
                                      mal::Norm(p2_0-rh2_global.m_Point) );
                            // rh contains feature_id=EID with valid extra barycentric coords wrt the DCR.E... we do NOT CARE about the actual DCR.F hit
                            cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,rh2_global.m_FeatureId.AsTriangle()), rh2_global.m_Extra_BarycentricCoords ),
                                       PointOnFeature( feature_id(eFT_Triangle,eid2), mal::Concat(p2_b,0) ) );
                        }
                    }
                    else
                    {
                        // Add CP symmetrized correspondences \todo CHECK IF VERY CLOSE AND ADD ONLY ONCE!!
                        if( bHit1 )
                        {
                            cd.AddCP( p1_0, rh1.m_Point,
                                      mal::SafeNormalized(rh1.m_Point-p1_0),
                                      mal::Norm(rh1.m_Point-p1_0) );
                            // rh contains feature_id=EID with valid extra barycentric coords wrt the DCR.E... we do NOT CARE about the actual DCR.F hit
                            cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,eid1), mal::Concat(p1_b,0) ),
                                       PointOnFeature( feature_id(eFT_Triangle,rh1.m_FeatureId.AsTriangle()), rh1.m_Extra_BarycentricCoords ) );
                        }
                        if( bHit2 )
                        {
                            cd.AddCP( rh2.m_Point, p2_0,
                                      mal::SafeNormalized(p2_0-rh2.m_Point),
                                      mal::Norm(p2_0-rh2.m_Point) );
                            // rh contains feature_id=EID with valid extra barycentric coords wrt the DCR.E... we do NOT CARE about the actual DCR.F hit
                            cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,rh2.m_FeatureId.AsTriangle()), rh2.m_Extra_BarycentricCoords ),
                                       PointOnFeature( feature_id(eFT_Triangle,eid2), mal::Concat(p2_b,0) ) );
                        }
                    }
                    // If ANY did not hit, add basic correspondence (\todo CONSIDER also adding closest point on IB, instead)
                    // \note THIS MAY happen if E1.m_AvgPos is OUTSIDE DCR2, even if E1 is partially inside it
                    if( p_context->m_DCR2DCR_RC_UseBasicIfNoHit && (!bHit1 || !bHit2) )
                    {
                        cd.AddCP( p1_0, p2_0,
                                  mal::SafeNormalized(p2_0-p1_0),
                                  mal::Norm(p2_0-p1_0) );
                        cd.AddPOF( PointOnFeature( feature_id(eFT_Triangle,eid1), mal::Concat(p1_b,0) ),
                                   PointOnFeature( feature_id(eFT_Triangle,eid2), mal::Concat(p2_b,0) ) );
                    }
                }
                break;
            default: break;
            }
        }
        cd.SetNumDisjointManifolds( num_crossings / 2 );
        return cd.End();
    }
    else
        return false;
}

/* TEMP: 2D and 3D floods are quite different, this are
   OLD RAMBLINGS, newest PhD.mm may invalidate any of this.

  - 3D:
    - Each Intersection Volume is bound by an Intersection
      Boundary on each object that intersect at one or
      several closed Intersection Curves formed by CFP.
      - Usually 1 closed IC, but in general an external IC
        can contain any number of internal IC that define
        "holes" (eg: a donut intersecting a plane will yield
        2 concentric IC)
   - Alg:
     - Find all IC by flooding of CFP
     - From each seed IC (may be an external IC or a hole IC)
     - Flood both obj INWARDS until Crossing features are
       found, marking them as flooded.
     - Remove ANY other seed IC that MAY have been flooded
       (ALL of their constituent CFP) from the explored
       seed
       => For simple external IC, no other IC will be removed,
       however, for complex IV bound by several IC, the flooding may
       start from ANY of them but will reach all other external/holes
       during flood, and will need to remmove them from open seeds.

  IMPORTANT: When redacting, emphasize that several opt make the
  algorithm practical in CPU & MEM usage, and that the naive "flood
  everything" approach would be very inefficient (global SID,
  etc...). Try to justify output-sensitivity.
*/


}} //namespace geo::np
