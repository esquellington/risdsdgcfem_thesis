#include "Contact_X_Y_BruteForce.h"
#include <Geo/shape/MeshSolidShape2.h>
#include <Geo/shape/TetSolidShape3.h>
#include "Intersection.h"
#include "RayCast.h"
#include "Overlap.h"

#ifdef __GEO_ENABLE_NP_SCRATCHPAD
#  include <util/ScopedAllocator.h>
#endif

#include <boost/bind.hpp>

namespace geo {
namespace np {

//--------------------------------------------------------------------------------------------------------------------------------
// MeshSolid2 Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

/*! Mesh Vs Plane */
bool Contact_MeshSolid2_Plane2_BruteForce( const MeshSolidShape2* p_mesh, const Transform2& mesh_tr, const Real* vec_dof,
                                           const Vec2& plane_normal, Real plane_coeff_d,
                                           ContactData2& cd, ContactCache2* p_cc,
                                           const np::Context* p_context )
{
    cd.Begin();
    //GEO_LOG_WARNING("Contact_MeshSolid2_Plane2_BruteForce %llx", (machine_uint_type)p_mesh_dof);

    //\todo THIS IS dangerous... Should have a generic SRV type with safe casts...
    const Vec2* vec_sdof( reinterpret_cast<const Vec2* >(vec_dof) );

    // For each boundary polygon
    for( unsigned int it_bp=0; it_bp < p_mesh->GetNumBoundaryP(); it_bp++ )
    {
        unsigned int it_he( p_mesh->BP_FirstHEID(it_bp) );
        // For each vertex
        do
        {
            feature_index_type vid( p_mesh->HE_OriginVID(it_he) );
            Vec2 pos_local( p_mesh->V_Pos( vid, vec_sdof ) );
            Vec2 pos_global( mesh_tr *  pos_local );
            Real dist( mal::Dot( pos_global, plane_normal ) + plane_coeff_d );
            if( dist < 0 ) //\todo use margin?!
            {
                cd.AddCP( pos_global, //on mesh
                          pos_global - dist*plane_normal, //on plane
                          plane_normal, -dist,
                          // Radius = 0.5 * ( lengt(prev_he) + length(he) )
                          Real(0.5)*  ( mal::Norm( pos_local - p_mesh->V_Pos( p_mesh->HE_OriginVID( p_mesh->HE_Prev(it_he) ), vec_sdof ) )
                                        + mal::Norm( pos_local - p_mesh->V_Pos( p_mesh->HE_FinalVID( it_he ), vec_sdof ) ) ) );
                //cd.AddFeaturePair( feature_id(eFT_Vertex,vid), feature_id() );
                cd.AddPOF( PointOnFeature( feature_id(eFT_Vertex,vid), Vec4(1,0,0,0) ),
                           PointOnFeature() );
            }
            it_he = p_mesh->HE_Next(it_he);
        } while ( it_he != p_mesh->BP_FirstHEID(it_bp) );
    }
    if( cd.Size() > 0 ) cd.SetNumDisjointManifolds( 1 ); //\todo Count it properly, for a deformable object it may be > 1
    return cd.End();
}

/*! Mesh Vs Mesh
  \todo By now it's not a fully functional Contact_X_Y method, CP are
  wrong, but binary overlap test and and m_NumDisjointManifolds are
  right
*/
bool Contact_MeshSolid2_MeshSolid2_BruteForce( const MeshSolidShape2* p_mesh1, const Transform2& mesh_tr1, const Real* vec_dof1,
                                               const MeshSolidShape2* p_mesh2, const Transform2& mesh_tr2, const Real* vec_dof2,
                                               ContactData2& cd, ContactCache2* p_cc,
                                               const np::Context* p_context )
{
    //GEO_LOG("Contact_MeshSolid2_MeshSolid2_BruteForce()");
    cd.Begin();

    //\todo THIS IS dangerous... Should have a generic SRV type with safe casts...
    const Vec2* vec_sdof1( reinterpret_cast<const Vec2* >(vec_dof1) );
    const Vec2* vec_sdof2( reinterpret_cast<const Vec2* >(vec_dof2) );
    Transform2 tr12( mal::Inverse(mesh_tr2) * mesh_tr1 );

    /* BruteForce: Test crossing for each boundary polygon pair, for
     each edge pair, working in mesh_tr2 ref sys
    */
    // For each boundary polygon bp1
    unsigned int num_crossings(0);
    for( unsigned int it_bp1=0; it_bp1 < p_mesh1->GetNumBoundaryP(); it_bp1++ )
    {
        // For each boundary polygon bp2
        for( unsigned int it_bp2=0; it_bp2 < p_mesh2->GetNumBoundaryP(); it_bp2++ )
        {
            unsigned int it_he1( p_mesh1->BP_FirstHEID(it_bp1) );
            // for each segment (a1,b1)
            do
            {
                Vec2 a1_2( tr12 * p_mesh1->V_Pos( p_mesh1->HE_OriginVID(it_he1), vec_sdof1 ) );
                Vec2 b1_2( tr12 * p_mesh1->V_Pos( p_mesh1->HE_FinalVID(it_he1), vec_sdof1 ) );
                unsigned int it_he2( p_mesh2->BP_FirstHEID(it_bp2) );
                // for each segment (a2,b2)
                do
                {
                    Vec2 a2_2( p_mesh2->V_Pos( p_mesh2->HE_OriginVID(it_he2), vec_sdof2 ) );
                    Vec2 b2_2( p_mesh2->V_Pos( p_mesh2->HE_FinalVID(it_he2), vec_sdof2 ) );
                    Real lambda1, lambda2;
                    if( Intersection_Segment2_Segment2( a1_2, b1_2, a2_2, b2_2, lambda1, lambda2 ) )
                    {
                        //TEMP: wrong point...
                        cd.AddCP( mesh_tr2*( a1_2 + lambda1 * (b1_2-a1_2) ),
                                  mesh_tr2*( a2_2 + lambda2 * (b2_2-a2_2) ),
                                  Vec2(0,1), Real(0) ); //\todo N,d hacked...
                        cd.AddPOF( PointOnFeature( feature_id( eFT_Segment, it_he1 ), Vec4(Real(1)-lambda1,lambda1,0,0) ),
                                   PointOnFeature( feature_id( eFT_Segment, it_he2 ), Vec4(Real(1)-lambda2,lambda2,0,0) ) );
                        num_crossings++;
                    }
                    it_he2 = p_mesh2->HE_Next(it_he2);
                } while ( it_he2 != p_mesh2->BP_FirstHEID(it_bp2) );
                it_he1 = p_mesh1->HE_Next(it_he1);
            } while ( it_he1 != p_mesh1->BP_FirstHEID(it_bp1) );
        }
    }
    cd.SetNumDisjointManifolds( num_crossings / 2 );
    return cd.End();
}

bool Contact_DCRMS2_DCRMS2_BruteForce_DCR( const MeshSolidShape2* p_mesh1, const Transform2& mesh_tr1, const Real* vec_dof1,
                                           const MeshSolidShape2* p_mesh2, const Transform2& mesh_tr2, const Real* vec_dof2,
                                           ContactData2& cd, ContactCache2* p_cc,
                                           const np::Context* p_context )
{
    cd.Begin();

    const DCR_MeshSolidShape2* pDCR1( p_mesh1->GetDCR() );
    const DCR_MeshSolidShape2* pDCR2( p_mesh2->GetDCR() );
    if( 0 == pDCR1 || 0 == pDCR2 ) return false; //\todo Fallback to a different method?

    //\todo THIS IS dangerous... Should have a generic SRV type with safe casts...
    const Vec2* vec_sdof1( reinterpret_cast<const Vec2*>(vec_dof1) );
    const Vec2* vec_sdof2( reinterpret_cast<const Vec2*>(vec_dof2) );
    const Vec2* default_sdof1( p_mesh1->GetVecDefaultSDOF() );
    const Vec2* default_sdof2( p_mesh2->GetVecDefaultSDOF() );
    Transform2 tr12( mal::Inverse(mesh_tr2) * mesh_tr1 );

    unsigned int num_crossings(0);
    {
#ifdef __GEO_ENABLE_NP_SCRATCHPAD
        // Eval deformed DCR vertices into p_context->m_ScratchPad, perform BF Contact determination and exit
        util::ScopedAllocator scoped_allocator( p_context->m_ScratchPad, "Contact_DCRMS2_DCRMS2_BruteForce_DCR" );
        Vec2* vec_v1_2 = scoped_allocator.NewArrayPOD<Vec2>(pDCR1->m_NumVertices);
        Vec2* vec_v2_2 = scoped_allocator.NewArrayPOD<Vec2>(pDCR2->m_NumVertices);
#else
        // Faked scope-allocation, original test, slower than true ScopedAllocator
        std::vector<Vec2> alloc_v1(pDCR1->m_NumVertices);
        std::vector<Vec2> alloc_v2(pDCR2->m_NumVertices);
        Vec2* vec_v1_2 = &alloc_v1[0];
        Vec2* vec_v2_2 = &alloc_v2[0];
        /* TEMP This works...
        std::vector<Vec2> alloc(pDCR1->m_NumVertices + pDCR2->m_NumVertices);
        Vec2* vec_v1_2 = &alloc[0];
        Vec2* vec_v2_2 = &alloc[pDCR1->m_NumVertices];
        */
#endif
        // Apply baricentric transform and tr12 to all V1
        for( unsigned int it_e=0; it_e<pDCR1->m_NumElements; it_e++ )
        {
            const DCR_MeshSolidShape2::Element& ed( pDCR1->m_vecE[it_e] );
            uint32 vec_nid[3];
            p_mesh1->P_VecVID(it_e,vec_nid,3);
            Mat2x2 B( mal::GMat2x2_From_Columns(vec_sdof1[vec_nid[1]]-vec_sdof1[vec_nid[0]],
                                                vec_sdof1[vec_nid[2]]-vec_sdof1[vec_nid[0]]) );
            Mat2x2 Bm( mal::GMat2x2_From_Columns(default_sdof1[vec_nid[1]]-default_sdof1[vec_nid[0]],
                                                 default_sdof1[vec_nid[2]]-default_sdof1[vec_nid[0]]) ); //\todo Could be precomputed in DCR::ED
            Mat2x2 B_invBm( B * mal::Inverse(Bm) );
            Transform2 tr12_B_invBm( tr12.m_Pos, tr12.m_Rot * B_invBm );
            Vec2 p0_2( tr12.m_Rot * vec_sdof1[vec_nid[0]] );
            for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
                vec_v1_2[ed.m_FirstVID + it_vie] = tr12_B_invBm * (pDCR1->m_vecV[ed.m_FirstVID + it_vie] - default_sdof1[vec_nid[0]]) + p0_2;
            GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed.m_NumVertices );
        }
        // Apply baricentric transform to all V2
        for( unsigned int it_e=0; it_e<pDCR2->m_NumElements; it_e++ )
        {
            const DCR_MeshSolidShape2::Element& ed( pDCR2->m_vecE[it_e] );
            uint32 vec_nid[3];
            p_mesh2->P_VecVID(it_e,vec_nid,3);
            Mat2x2 B( mal::GMat2x2_From_Columns(vec_sdof2[vec_nid[1]]-vec_sdof2[vec_nid[0]],
                                                vec_sdof2[vec_nid[2]]-vec_sdof2[vec_nid[0]]) );
            Mat2x2 Bm( mal::GMat2x2_From_Columns(default_sdof2[vec_nid[1]]-default_sdof2[vec_nid[0]],
                                                 default_sdof2[vec_nid[2]]-default_sdof2[vec_nid[0]]) ); //\todo Could be precomputed in DCR::ED
            Mat2x2 B_invBm( B * mal::Inverse(Bm) );
            Vec2 p0_2( vec_sdof2[vec_nid[0]] );
            for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
                vec_v2_2[ed.m_FirstVID + it_vie] = B_invBm * (pDCR2->m_vecV[ed.m_FirstVID + it_vie] - default_sdof2[vec_nid[0]]) + p0_2;
            GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed.m_NumVertices );
        }
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
        // //if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_IP ) )
        // {
        //     for( unsigned int it_v=0; it_v<pDCR1->m_NumVertices; it_v++ )
        //         cd.m_VD.m_vecCP.push_back( std::make_pair(vec_v1_2[it_v],mesh_tr2 * vec_v1_2[it_v]) );
        //     for( unsigned int it_v=0; it_v<pDCR2->m_NumVertices; it_v++ )
        //         cd.m_VD.m_vecCP.push_back( std::make_pair(vec_v2_2[it_v],mesh_tr2 * vec_v2_2[it_v]) );
        // }
#endif
        //BruteForce: for each patch pair, for each segment pair...
        // For each patch1
        for( unsigned int it_patch1=0; it_patch1 < pDCR1->m_NumPatches; it_patch1++ )
        {
            const geo::DCR_MeshSolidShape2::Patch& pd1( pDCR1->m_vecP[it_patch1] );
            // For each patch2
            for( unsigned int it_patch2=0; it_patch2 < pDCR2->m_NumPatches; it_patch2++ )
            {
                const geo::DCR_MeshSolidShape2::Patch& pd2( pDCR2->m_vecP[it_patch2] );
                // for each segment (a1,b1)
                for( unsigned int it_sid1=0; it_sid1 < pd1.m_NumSegments; it_sid1++ )
                {
                    unsigned int sid1( pd1.m_FirstSID + it_sid1 );
                    Vec2 a1_2( vec_v1_2[ pDCR1->m_vecS[sid1].GetVID(0) ] );
                    Vec2 b1_2( vec_v1_2[ pDCR1->m_vecS[sid1].GetVID(1) ] );
                    //\todo IMPORTANT! Early-out if (a1,b1) does NOT overlap bv2!! \todo use bslab2 here, it's FAST to test and provides the smallest volume among all BV types
                    // for each segment (a2,b2)
                    for( unsigned int it_sid2=0; it_sid2 < pd2.m_NumSegments; it_sid2++ )
                    {
                        unsigned int sid2( pd2.m_FirstSID + it_sid2 );
                        Vec2 a2_2( vec_v2_2[ pDCR2->m_vecS[sid2].GetVID(0) ] );
                        Vec2 b2_2( vec_v2_2[ pDCR2->m_vecS[sid2].GetVID(1) ] );
                        Real lambda1, lambda2;
                        if( Intersection_Segment2_Segment2( a1_2, b1_2, a2_2, b2_2, lambda1, lambda2 ) )
                        {
                            //TEMP: wrong point...
                            cd.AddCP( mesh_tr2*( a1_2 + lambda1 * (b1_2-a1_2) ),
                                      mesh_tr2*( a2_2 + lambda2 * (b2_2-a2_2) ),
                                      Vec2(0,1), Real(0) ); //\todo N,d hacked...
                            //IMPORTANT: POF SHOULD BE be relative to the MSS, NOT the DCR features
                            cd.AddPOF( PointOnFeature( feature_id( eFT_Segment, sid1 ), Vec4(Real(1)-lambda1,lambda1,0,0) ),
                                       PointOnFeature( feature_id( eFT_Segment, sid2 ), Vec4(Real(1)-lambda2,lambda2,0,0) ) );
                            num_crossings++;
                        }
                    }
                }
            }
        }
    }
    cd.SetNumDisjointManifolds( num_crossings / 2 );
    return cd.End();
}

/* Test all DCR2DCR element pairs and, for overlapping (e1,e2), test
   all segment pairs (s1,s2) in all patch pairs (patch1,patch2)
   \note Notice that NOT all MSS elements are considered, only the
         subset used by their respective DCR, which ideally
         corresponding to MSS layer[0]
*/
bool Contact_DCRMS2_DCRMS2_BruteForce_MSS_Lazy_DCR( const MeshSolidShape2* p_mesh1, const Transform2& mesh_tr1, const Real* vec_dof1,
                                                    const MeshSolidShape2* p_mesh2, const Transform2& mesh_tr2, const Real* vec_dof2,
                                                    ContactData2& cd, ContactCache2* p_cc,
                                                    const np::Context* p_context )
{
    cd.Begin();

    typedef bv::AABB2 bv_type;
    //typedef bv::DOP2_K4 bv_type;
    //typedef bv::DOP2_K8 bv_type;

    const DCR_MeshSolidShape2* pDCR1( p_mesh1->GetDCR() );
    const DCR_MeshSolidShape2* pDCR2( p_mesh2->GetDCR() );
    if( 0 == pDCR1 || 0 == pDCR2 ) return false; //\todo Fallback to a different method?

    //\todo THIS IS dangerous... Should have a generic SRV type with safe casts...
    const Vec2* vec_sdof1( reinterpret_cast<const Vec2*>(vec_dof1) );
    const Vec2* vec_sdof2( reinterpret_cast<const Vec2*>(vec_dof2) );
    const Vec2* default_sdof1( p_mesh1->GetVecDefaultSDOF() );
    const Vec2* default_sdof2( p_mesh2->GetVecDefaultSDOF() );
    Transform2 tr12( mal::Inverse(mesh_tr2) * mesh_tr1 );

    /* TEMP: compute mesh2 BV (recomputation, it was ALREADY
             AVAILABLE in broad-phase/mid-phase, consider propagating
             it to narrow-phase (as a ContactCache memeber, for
             example)
             \todo Alternatively, compute and cache the array of e2 BV
             and the global mesh2 BV in a single pass, and reuse the
             cache in (e1,e2) overlap tests.
    */
    bv_type bv2;
    p_mesh2->ComputeBVD( bv2, mesh_tr2, vec_sdof2 );

    //\todo Eval deformed DCR vertices into p_context->m_ScratchPad, perform BF Contact determination and exit
    //\todo p_context->m_ScratchPad.BeginScope( "Contact_DCRMS2_DCRMS2_BruteForce_MSS_Lazy_DCR" );
    unsigned int num_crossings(0);
    {
        // for each overlapping element pair test all patch pairs for intersection
        for( unsigned int it_e1=0; it_e1<pDCR1->m_NumElements; it_e1++ )
        {
            const DCR_MeshSolidShape2::Element& ed1( pDCR1->m_vecE[it_e1] );
            uint32 vec_nid1[3];
            p_mesh1->P_VecVID( it_e1, vec_nid1, 3 ); //\todo NID could be stored in DCR::ED if not directly available in MSS
            // Test e1 vs p_mesh2
            bv_type bv_e1( bv_type( mesh_tr1*vec_sdof1[vec_nid1[0]] )
                            .Merge( mesh_tr1*vec_sdof1[vec_nid1[1]] )
                            .Merge( mesh_tr1*vec_sdof1[vec_nid1[2]] ) );
            if( TestOverlap(bv_e1,bv2) )//\todo Ideally: Overlap( ed1, mesh2 )
            {
#ifdef __GEO_ENABLE_NP_SCRATCHPAD
                // Transform ALL vertices e1 and cache
                util::ScopedAllocator scoped_allocator_e1( p_context->m_ScratchPad, "Contact_DCRMS2_DCRMS2_BruteForce_MSS_Lazy_DCR_e1" );
                Vec2* vec_v1_2 = scoped_allocator_e1.NewArrayPOD<Vec2>(ed1.m_NumVertices);
#else
                // Transform ALL vertices e1 and cache \todo If frame/scope allocated, this would be FREE, now it's EXPENSIVE
                std::vector<Vec2> alloc_v1(ed1.m_NumVertices);
                Vec2* vec_v1_2 = &alloc_v1[0];
#endif
                {
                    Mat2x2 B( mal::GMat2x2_From_Columns(vec_sdof1[vec_nid1[1]]-vec_sdof1[vec_nid1[0]],
                                                        vec_sdof1[vec_nid1[2]]-vec_sdof1[vec_nid1[0]]) );
                    Mat2x2 Bm( mal::GMat2x2_From_Columns(default_sdof1[vec_nid1[1]]-default_sdof1[vec_nid1[0]],
                                                         default_sdof1[vec_nid1[2]]-default_sdof1[vec_nid1[0]]) ); //\todo Could be precomputed in DCR::ED
                    Mat2x2 B_invBm( B * mal::Inverse(Bm) );
                    Transform2 tr12_B_invBm( tr12.m_Pos, tr12.m_Rot * B_invBm );
                    Vec2 p0_2( tr12.m_Rot * vec_sdof1[vec_nid1[0]] );
                    // Transform ed1 m_NumVertices vertices from m_First and store them in 0-based indices inside vec_v1_2
                    for( unsigned int it_vie=0; it_vie<ed1.m_NumVertices; it_vie++ )
                        vec_v1_2[it_vie] = tr12_B_invBm * (pDCR1->m_vecV[ed1.m_FirstVID + it_vie] - default_sdof1[vec_nid1[0]]) + p0_2;
                    GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed1.m_NumVertices );
                }
                // For each e2
                for( unsigned int it_e2=0; it_e2<pDCR2->m_NumElements; it_e2++ )
                {
                    const DCR_MeshSolidShape2::Element& ed2( pDCR2->m_vecE[it_e2] );
                    uint32 vec_nid2[3];
                    p_mesh2->P_VecVID( it_e2, vec_nid2, 3 ); //\todo NID could be stored in DCR::ED if not directly available in MSS
                    // Test e1 vs e2 //TEMP: Using element bv, not exact geom
                    //\todo Ideally, TestOverlap( e1, e2 ) using precomputed per-element or per-patch BV in barycentric coords (eg: BDOP)
                    bv_type bv_e2( bv_type( mesh_tr2*vec_sdof2[vec_nid2[0]] )
                                   .Merge( mesh_tr2*vec_sdof2[vec_nid2[1]] )
                                   .Merge( mesh_tr2*vec_sdof2[vec_nid2[2]] ) );
                    if( TestOverlap(bv_e1,bv_e2) )
                    {
#ifdef __GEO_ENABLE_NP_SCRATCHPAD
                        // Transform ALL vertices e1 and cache
                        util::ScopedAllocator scoped_allocator_e2( p_context->m_ScratchPad, "Contact_DCRMS2_DCRMS2_BruteForce_MSS_Lazy_DCR_e2" );
                        Vec2* vec_v2_2 = scoped_allocator_e2.NewArrayPOD<Vec2>(ed2.m_NumVertices);
#else
                        // Transform ALL vertices e2 and cache \todo If frame/scope allocated, this would be FREE, now it's EXPENSIVE
                        std::vector<Vec2> alloc_v2(ed2.m_NumVertices);
                        Vec2* vec_v2_2 = &alloc_v2[0];
#endif
                        // e2
                        /*\note Here we're re-transforming e2 for EACH
                          overlapping e1, consider caching SOME
                          transformed e2 (LRU?) but avoid storing ALL
                          vertices, as #V >> #E
                        */
                        {
                            Mat2x2 B( mal::GMat2x2_From_Columns(vec_sdof2[vec_nid2[1]]-vec_sdof2[vec_nid2[0]],
                                                                vec_sdof2[vec_nid2[2]]-vec_sdof2[vec_nid2[0]]) );
                            Mat2x2 Bm( mal::GMat2x2_From_Columns(default_sdof2[vec_nid2[1]]-default_sdof2[vec_nid2[0]],
                                                                 default_sdof2[vec_nid2[2]]-default_sdof2[vec_nid2[0]]) ); //\todo Could be precomputed in DCR::ED
                            Mat2x2 B_invBm( B * mal::Inverse(Bm) );
                            Vec2 p0_2( vec_sdof2[vec_nid2[0]] );
                            for( unsigned int it_vie=0; it_vie<ed2.m_NumVertices; it_vie++ )
                                vec_v2_2[it_vie] = B_invBm * (pDCR2->m_vecV[ed2.m_FirstVID + it_vie] - default_sdof2[vec_nid2[0]]) + p0_2;
                            GEO_NP_STAT_ADD( contact.m_Num_Vertices_Transformed_Barycentric, ed2.m_NumVertices );
                        }
                        //-- Test all patch pairs using cached transformed vertices
                        // For each patch1 in e1
                        for( unsigned int it_patch1=0; it_patch1 < ed1.m_NumPatches; it_patch1++ )
                        {
                            const geo::DCR_MeshSolidShape2::Patch& pd1( pDCR1->m_vecP[ed1.m_FirstPID + it_patch1] );
                            // For each patch2 in e2
                            for( unsigned int it_patch2=0; it_patch2 < ed2.m_NumPatches; it_patch2++ )
                            {
                                const geo::DCR_MeshSolidShape2::Patch& pd2( pDCR2->m_vecP[ed2.m_FirstPID + it_patch2] );
                                // for each segment (a1,b1)
                                for( unsigned int it_sid1=0; it_sid1 < pd1.m_NumSegments; it_sid1++ )
                                {
                                    unsigned int sid1( pd1.m_FirstSID + it_sid1 );
                                    Vec2 a1_2( vec_v1_2[ pDCR1->m_vecS[sid1].GetVID(0) - ed1.m_FirstVID ] );
                                    Vec2 b1_2( vec_v1_2[ pDCR1->m_vecS[sid1].GetVID(1) - ed1.m_FirstVID ] );
                                    //\todo IMPORTANT! Early-out if (a1,b1) does NOT overlap bv2!! \todo use bslab2 here, it's FAST to test and provides the smallest volume among all BV types
                                    // for each segment (a2,b2)
                                    for( unsigned int it_sid2=0; it_sid2 < pd2.m_NumSegments; it_sid2++ )
                                    {
                                        unsigned int sid2( pd2.m_FirstSID + it_sid2 );
                                        Vec2 a2_2( vec_v2_2[ pDCR2->m_vecS[sid2].GetVID(0) - ed2.m_FirstVID ] );
                                        Vec2 b2_2( vec_v2_2[ pDCR2->m_vecS[sid2].GetVID(1) - ed2.m_FirstVID ] );
                                        Real lambda1, lambda2;
                                        if( Intersection_Segment2_Segment2( a1_2, b1_2, a2_2, b2_2, lambda1, lambda2 ) )
                                        {
                                            //TEMP: wrong point...
                                            cd.AddCP( mesh_tr2*( a1_2 + lambda1 * (b1_2-a1_2) ),
                                                      mesh_tr2*( a2_2 + lambda2 * (b2_2-a2_2) ),
                                                      Vec2(0,1), Real(0) ); //\todo N,d hacked...
                                            //IMPORTANT: POF SHOULD BE be relative to the MSS, NOT the DCR
                                            cd.AddPOF( PointOnFeature( feature_id( eFT_Segment, sid1 ), Vec4(Real(1)-lambda1,lambda1,0,0) ),
                                                       PointOnFeature( feature_id( eFT_Segment, sid2 ), Vec4(Real(1)-lambda2,lambda2,0,0) ) );
                                            num_crossings++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    cd.SetNumDisjointManifolds( num_crossings / 2 );
    return cd.End();
}

//--------------------------------------------------------------------------------------------------------------------------------
// TetSolid3 Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_TetSolid3_Plane3_BruteForce( const TetSolidShape3* p_solid, const Transform3& mesh_tr, const Real* vec_dof,
                                          const Vec3& plane_normal, Real plane_coeff_d,
                                          ContactData3& cd, ContactCache3* p_cc,
                                          const np::Context* p_context )
{
    cd.Begin();
    //GEO_LOG_WARNING("Contact_MeshSolid2_Plane2_BruteForce %llx", (machine_uint_type)p_mesh_dof);

    //\todo THIS IS dangerous... Should have a generic SRV type with safe casts...
    const Vec3* vec_sdof( reinterpret_cast<const Vec3* >(vec_dof) );

    // For each boundary polygon
    for( unsigned int it_v=0; it_v < p_solid->GetNumV(); it_v++ )
    {
        Vec3 pos_local( p_solid->V_Pos( it_v, vec_sdof ) );
        Vec3 pos_global( mesh_tr *  pos_local );
        Real dist( mal::Dot( pos_global, plane_normal ) + plane_coeff_d );
        if( dist < 0 ) //\todo use margin?!
        {
            cd.AddCP( pos_global, //on mesh
                      pos_global - dist*plane_normal, //on plane
                      plane_normal, -dist,
                      1.0f ); //\todo FAKE RADIUS
            cd.AddPOF( PointOnFeature( feature_id(eFT_Vertex,it_v), Vec4(1,0,0,0) ),
                       PointOnFeature() );
            GEO_NP_STAT_INC( mig2015.m_Num_CP );
        }
    }
    if( cd.Size() > 0 ) cd.SetNumDisjointManifolds( 1 ); //\todo Count it properly, for a deformable object it may be > 1
    return cd.End();
}

}} //namespace geo::np
