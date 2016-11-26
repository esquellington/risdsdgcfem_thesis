#include "Vertices3_In_Tetrahedrons3_Barycentric_Embedding.h"
#include <Geo/shape/shape.h> //Include all shapes
#include <map>
#include <vector>

#define __USE_BVH_IN_EMBEDDING_V3INT3BE
#ifdef __USE_BVH_IN_EMBEDDING_V3INT3BE
#  include <boost/bind.hpp>
#endif

namespace geo {

void Vertices3_In_Tetrahedrons3_Barycentric_Embedding::SetBakedData( bool b_shared,
                                                                     const Transform3& transform_e2d,
                                                                     unsigned int num_tevs, unsigned num_embedded_vid,
                                                                     const BakedTetrahedronEmbeddedVertexSubset* vec_tevs, const embedded_feature_index_type* vec_embedded_vid )
{
    GEO_ASSERT( b_shared ); //\todo CAN ONLY BE SHARED BY NOW, NO m_pBuffer single allocation yet
    GEO_ASSERT( num_tevs < cDeformer_InvalidFeatureIndex );
    GEO_ASSERT( num_embedded_vid < cEmbedded_InvalidFeatureIndex );
    m_TransformE2D = transform_e2d;
    m_NumTEVS = num_tevs;
    m_NumEmbeddedVID = num_embedded_vid;
    m_vecTEVS = vec_tevs;
    m_vecEmbeddedVID = vec_embedded_vid;
}

void Vertices3_In_Tetrahedrons3_Barycentric_Embedding::ClearBakedData()
{
    m_TransformE2D = Transform3::Identity();
    m_NumTEVS = 0;
    m_NumEmbeddedVID = 0;
    m_vecTEVS = 0;
    m_vecEmbeddedVID = 0;
}

void Vertices3_In_Tetrahedrons3_Barycentric_Embedding::Apply( const IObject* p_deformer, IObject* p_embedded ) const
{
    //-- Apply embedding transform
    static_cast<IObject3*>(p_embedded)->SetTransform_Embedded( static_cast<const IObject3*>(p_deformer)->GetTransform() * m_TransformE2D );

    //-- Apply embedding DOF
    const Vec3* vec_deformer_sdof( reinterpret_cast<const Vec3*>( p_deformer->GetVecDOF() ) ); //GetVecDOF() updates deformer DOF if it's also embedded and out of date
    const Vec3* vec_embedded_default_sdof( reinterpret_cast<const Vec3*>( const_cast<const IObject*>(p_embedded)->GetShapeInterface()->GetVecDefaultDOF() ) );
    Vec3* vec_embedded_sdof( reinterpret_cast<Vec3*>( p_embedded->GetVecDOF_WriteOnly() ) ); //\note GetVecDOF_WriteOnly() just returns pointer, otherwise cyclic SyncEmbedding()->Apply()->SyncEmbedding()!
    //foreach TEVS: foreach embedded vid: compute deformed pos and write to vec_dof_embedded
    for( unsigned int it_tevs=0; it_tevs < m_NumTEVS; it_tevs++ )
    {
        const BakedTetrahedronEmbeddedVertexSubset& btevs( m_vecTEVS[it_tevs] );
        // Compute deformer matrix
        Vec3 vec_points[4];
        vec_points[0] = vec_deformer_sdof[ btevs.m_vecTetrahedronVID[0] ];
        vec_points[1] = vec_deformer_sdof[ btevs.m_vecTetrahedronVID[1] ];
        vec_points[2] = vec_deformer_sdof[ btevs.m_vecTetrahedronVID[2] ];
        vec_points[3] = vec_deformer_sdof[ btevs.m_vecTetrahedronVID[3] ];
        Mat4x4 B( 1, 1, 1, 1,
                  vec_points[0].x(), vec_points[1].x(), vec_points[2].x(), vec_points[3].x(),
                  vec_points[0].y(), vec_points[1].y(), vec_points[2].y(), vec_points[3].y(),
                  vec_points[0].z(), vec_points[1].z(), vec_points[2].z(), vec_points[3].z() );
        /*\todo THIS IS WRONG, because m_InvTetrahedronBCM was
         * computed in global coords, but B is in deformer
         * local coords!!... seems to work if both Tr are Id
         * at embedding-time, but I'm afraid it WOULD FAIL if
         * any object was rotated at embedding time... we
         * should use m_TransformE2D to fix it, or, better,
         * store baked m_InvTetrahedronBCM in the properly
         * transformed ref sys. This problem may also apply to
         * B, write it down.
         */
        Mat4x4 M( B * btevs.m_InvTetrahedronBCM );
        mal::GMat<Real,3,4> m( mal::GRange<1,0,3,3>( M ) );
        for( unsigned int it_ev=0; it_ev < btevs.m_NumEmbeddedVertices; it_ev++ )
        {
            embedded_feature_index_type vid( m_vecEmbeddedVID[ btevs.m_FirstEmbeddedVertexIndex + it_ev ] );
            /*TEMP: Slow way
              Vec4 v4 = M * Vec4( 1, vec_embedded_default_sdof[vid].x(), vec_embedded_default_sdof[vid].y(), vec_embedded_default_sdof[vid].z() );
              vec_embedded_sdof[ vid ] = Vec3( v4[1], v4[2], v[3] );
            */
            vec_embedded_sdof[ vid ] = mal::GColumn<0>(m) + mal::GRange<0,1,2,3>(m) * vec_embedded_default_sdof[vid]; // == M * [1 x y z]^T, faster
        }
    }
}

struct TetrahedronEmbeddedVertexSubset
{
    deformer_feature_index_type m_TetrahedronID;
    deformer_feature_index_type m_vecTetrahedronVID[4];
    Mat4x4 m_InvTetrahedronBCM;
    std::vector<embedded_feature_index_type> m_vecEmbeddedVID;
};

/* foreach vertex v in polygonal
   find closest or embedding element e (use distance from barycentric to clamped barycentric coords as measure)
   add to e vertex subset
   if v outside e, report warning
   \todo
   bake "crust" TEVS, with a single global vertex array
   consider reordering vertices
   \todo Consider using Object DOF instead of default Shape DOF
   for the embedding. This would allow user translation/rotation
   to match E and S better.
*/
void ComputeEmbedding_TriSurface3_In_TetSolid3_As_Vertices3_In_Tetrahedrons3_Barycentric( const TetSolidShape3& solid, const Transform3& tr_deformer,
                                                                                          const TriSurfaceShape3& surface, const Transform3& tr_embedded,
                                                                                          std::map< deformer_feature_index_type, TetrahedronEmbeddedVertexSubset* >& map_EID2TEVS )
{
#ifdef __USE_BVH_IN_EMBEDDING_V3INT3BE
    GEO_ASSERT( map_EID2TEVS.empty() );
    /* IMPORTANT Ideally, we'd retrieve the BVH from solid, but this may or may not be up to date/in sync with tr_deformer... so we'll just rebuild it here.
    const BVH_TetSolidShape3* pSolidBVH( solid.GetBVH() );
    GEO_ASSERT( pSolidBVH );
    */
    //TEMP: build a local BVH with solid's given transform and default DOF
    BVH_TetSolidShape3 bvh;
    bvh.Rebuild_TopDown( solid.GetNumT(),
                         boost::bind<void>( &GEBV_TetSolidShape3_E<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
                                            &solid, tr_deformer, solid.GetVecDefaultSDOF(),
                                            _1, _2) );
    unsigned int num_unmatched(0);
    for( unsigned int it_v=0; it_v < surface.GetNumV(); it_v++ )
    {
        Vec3 vertex_pos( tr_embedded * surface.V_Pos_0(it_v) );
        // Find best embedding element e
        deformer_feature_index_type best_TetrahedronID( cDeformer_InvalidFeatureIndex );
        deformer_feature_index_type best_vecTetrahedronVID[4] = { cDeformer_InvalidFeatureIndex,
                                                                  cDeformer_InvalidFeatureIndex,
                                                                  cDeformer_InvalidFeatureIndex,
                                                                  cDeformer_InvalidFeatureIndex };
        Real best_barycentric_distance_sq( mal::Infinity<Real>() );
        Mat4x4 best_inv_BCM( Mat4x4::Identity() );
        std::vector< BVH_TetSolidShape3::entry_index_type > vec_overlaps;
        auto tetBV = BVH_TetSolidShape3::bv_type(vertex_pos);
        tetBV.Extend(0.01f);
        if( bvh.Test( tetBV, vec_overlaps ) )
        {
            for( unsigned int it_o=0;
                 it_o < vec_overlaps.size() && best_barycentric_distance_sq > 0; //stop if vertex is exactly inside an element
                 it_o++ )
            {
                uint32 tid( vec_overlaps[it_o] );
                // Gather element nodes
                unsigned int vid0( solid.T_VID(tid,0) );
                unsigned int vid1( solid.T_VID(tid,1) );
                unsigned int vid2( solid.T_VID(tid,2) );
                unsigned int vid3( solid.T_VID(tid,3) );
                Vec3 vec_element_nodes[4];
                vec_element_nodes[0] = tr_deformer * solid.V_Pos_0( vid0 );
                vec_element_nodes[1] = tr_deformer * solid.V_Pos_0( vid1 );
                vec_element_nodes[2] = tr_deformer * solid.V_Pos_0( vid2 );
                vec_element_nodes[3] = tr_deformer * solid.V_Pos_0( vid3 );
                // Compute element barycentric coordinates matrix in global coords
                Mat4x4 BCM( 1, 1, 1, 1,
                            vec_element_nodes[0].x(), vec_element_nodes[1].x(), vec_element_nodes[2].x(), vec_element_nodes[3].x(),
                            vec_element_nodes[0].y(), vec_element_nodes[1].y(), vec_element_nodes[2].y(), vec_element_nodes[3].y(),
                            vec_element_nodes[0].z(), vec_element_nodes[1].z(), vec_element_nodes[2].z(), vec_element_nodes[3].z() );
                Mat4x4 inv_BCM( mal::Inverse( BCM ) );
                // Compute vertex barycentric coords (may be out of range)
                Vec4 barycentric_coords( inv_BCM * mal::Concat(1,vertex_pos) );
                Vec4 clamped_barycentric_coords = mal::Clamp( barycentric_coords,
                                                              Vec4(0,0,0,0),
                                                              Vec4(1,1,1,1) );
                Real barycentric_distance_sq = mal::NormSq( barycentric_coords - clamped_barycentric_coords );
                if( best_TetrahedronID == cDeformer_InvalidFeatureIndex
                    || barycentric_distance_sq < best_barycentric_distance_sq )
                {
                    best_TetrahedronID = tid;
                    best_vecTetrahedronVID[0] = vid0; best_vecTetrahedronVID[1] = vid1; best_vecTetrahedronVID[2] = vid2; best_vecTetrahedronVID[3] = vid3;
                    best_barycentric_distance_sq = barycentric_distance_sq;
                    best_inv_BCM = inv_BCM;
                }
            }
            // Create best e TEVS if required
            if( map_EID2TEVS.find( best_TetrahedronID ) == map_EID2TEVS.end() )
            {
                TetrahedronEmbeddedVertexSubset *pTEVS = new TetrahedronEmbeddedVertexSubset();
                pTEVS->m_TetrahedronID = best_TetrahedronID;
                for( int i=0; i<4; i++ ) pTEVS->m_vecTetrahedronVID[i] = best_vecTetrahedronVID[i];
                pTEVS->m_InvTetrahedronBCM = best_inv_BCM;
                map_EID2TEVS.insert( std::make_pair( best_TetrahedronID, pTEVS ) );
            }
            // Add v to e
            map_EID2TEVS[best_TetrahedronID]->m_vecEmbeddedVID.push_back( it_v );
        }
        else
        {
            // Create TEVS 0 if required
            if( map_EID2TEVS.find( 0 ) == map_EID2TEVS.end() )
            {
                TetrahedronEmbeddedVertexSubset *pTEVS = new TetrahedronEmbeddedVertexSubset();
                pTEVS->m_TetrahedronID = 0;
                for( int i=0; i<4; i++ ) pTEVS->m_vecTetrahedronVID[i] = solid.T_VID(0,i);
                pTEVS->m_InvTetrahedronBCM = Mat4x4::Identity();
                map_EID2TEVS.insert( std::make_pair( 0, pTEVS ) );
            }
            // Add v to TEVS 0
            map_EID2TEVS[0]->m_vecEmbeddedVID.push_back( it_v );
            num_unmatched++;
        }
    }
    if( num_unmatched > 0 )
        GEO_LOG_ERROR("ComputeEmbedding_TriSurface3_In_TetSolid3_As_Vertices3_In_Tetrahedrons3_Barycentric() %u unmatched V associated to Tet[0]", num_unmatched );
#else //__USE_BVH_IN_EMBEDDING_V3INT3BE
    GEO_ASSERT( map_EID2TEVS.empty() );
    for( unsigned int it_v=0; it_v < surface.GetNumV(); it_v++ )
    {
        Vec3 vertex_pos( tr_embedded * surface.V_Pos_0(it_v) );
        // Find best embedding element e
        deformer_feature_index_type best_TetrahedronID( cDeformer_InvalidFeatureIndex );
        deformer_feature_index_type best_vecTetrahedronVID[4] = { cDeformer_InvalidFeatureIndex,
                                                                  cDeformer_InvalidFeatureIndex,
                                                                  cDeformer_InvalidFeatureIndex,
                                                                  cDeformer_InvalidFeatureIndex };
        Real best_barycentric_distance_sq( mal::Infinity<Real>() );
        Mat4x4 best_inv_BCM( Mat4x4::Identity() );
        for( unsigned int it_tet=0;
             it_tet < solid.GetNumT() && best_barycentric_distance_sq > 0; //stop if vertex is exactly inside an element
             it_tet++ )
        {
            // Gather element nodes
            unsigned int vid0( solid.T_VID(it_tet,0) );
            unsigned int vid1( solid.T_VID(it_tet,1) );
            unsigned int vid2( solid.T_VID(it_tet,2) );
            unsigned int vid3( solid.T_VID(it_tet,3) );
            Vec3 vec_element_nodes[4];
            vec_element_nodes[0] = tr_deformer * solid.V_Pos_0( vid0 );
            vec_element_nodes[1] = tr_deformer * solid.V_Pos_0( vid1 );
            vec_element_nodes[2] = tr_deformer * solid.V_Pos_0( vid2 );
            vec_element_nodes[3] = tr_deformer * solid.V_Pos_0( vid3 );
            // Compute element barycentric coordinates matrix in global coords
            Mat4x4 BCM( 1, 1, 1, 1,
                        vec_element_nodes[0].x(), vec_element_nodes[1].x(), vec_element_nodes[2].x(), vec_element_nodes[3].x(),
                        vec_element_nodes[0].y(), vec_element_nodes[1].y(), vec_element_nodes[2].y(), vec_element_nodes[3].y(),
                        vec_element_nodes[0].z(), vec_element_nodes[1].z(), vec_element_nodes[2].z(), vec_element_nodes[3].z() );
            Mat4x4 inv_BCM( mal::Inverse( BCM ) );
            // Compute vertex barycentric coords (may be out of range)
            Vec4 barycentric_coords( inv_BCM * mal::Concat(1,vertex_pos) );
            Vec4 clamped_barycentric_coords = mal::Clamp( barycentric_coords,
                                                          Vec4(0,0,0,0),
                                                          Vec4(1,1,1,1) );
            Real barycentric_distance_sq = mal::NormSq( barycentric_coords - clamped_barycentric_coords );
            if( best_TetrahedronID == cDeformer_InvalidFeatureIndex
                || barycentric_distance_sq < best_barycentric_distance_sq )
            {
                best_TetrahedronID = it_tet;
                best_vecTetrahedronVID[0] = vid0; best_vecTetrahedronVID[1] = vid1; best_vecTetrahedronVID[2] = vid2; best_vecTetrahedronVID[3] = vid3;
                best_barycentric_distance_sq = barycentric_distance_sq;
                best_inv_BCM = inv_BCM;
            }
        }
        // Create best e TEVS if required
        if( map_EID2TEVS.find( best_TetrahedronID ) == map_EID2TEVS.end() )
        {
            TetrahedronEmbeddedVertexSubset *pTEVS = new TetrahedronEmbeddedVertexSubset();
            pTEVS->m_TetrahedronID = best_TetrahedronID;
            for( int i=0; i<4; i++ ) pTEVS->m_vecTetrahedronVID[i] = best_vecTetrahedronVID[i];
            pTEVS->m_InvTetrahedronBCM = best_inv_BCM;
            map_EID2TEVS.insert( std::make_pair( best_TetrahedronID, pTEVS ) );
        }
        // Add v to e
        map_EID2TEVS[best_TetrahedronID]->m_vecEmbeddedVID.push_back( it_v );
    }
#endif
}

void Editable_Vertices3_In_Tetrahedrons3_Barycentric_Embedding::Init( const IObject3* p_deformer, const IObject3* p_embedded )
{
    Clear();

    Transform3 transform_e2d = mal::Inverse( p_deformer->GetTransform() ) * p_embedded->GetTransform();

    // Compute embedding \todo SHAPE-SPECIFIC, could use a double switch, or a DD table, whatever...
    unsigned int num_embedded_vid = static_cast<const TriSurfaceShape3*>( p_embedded->GetShapeInterface() )->GetNumV();
    std::map< deformer_feature_index_type, TetrahedronEmbeddedVertexSubset* > map_EID2TEVS;
    ComputeEmbedding_TriSurface3_In_TetSolid3_As_Vertices3_In_Tetrahedrons3_Barycentric( *static_cast<const TetSolidShape3*>( p_deformer->GetShapeInterface() ), p_deformer->GetTransform(),
                                                                                         *static_cast<const TriSurfaceShape3*>( p_embedded->GetShapeInterface() ), p_embedded->GetTransform(),
                                                                                         map_EID2TEVS );

    // Bake embedding \todo NOT SHAPE-SPECIFIC, but embedding-type specific (V2inT2)
    m_NumTEVS = map_EID2TEVS.size();
    m_NumEmbeddedVID = num_embedded_vid;
    m_allocTEVS = new BakedTetrahedronEmbeddedVertexSubset[ m_NumTEVS ];
    m_allocEmbeddedVID = new embedded_feature_index_type[ m_NumEmbeddedVID ];
    unsigned int last_embedded_vertex_index(0);
    unsigned int last_btevs(0);
    for( std::map< deformer_feature_index_type, TetrahedronEmbeddedVertexSubset* >::iterator it_map = map_EID2TEVS.begin();
         it_map != map_EID2TEVS.end();
         ++it_map )
    {
        const TetrahedronEmbeddedVertexSubset& tevs = *it_map->second;
        BakedTetrahedronEmbeddedVertexSubset& btevs = m_allocTEVS[last_btevs];
        btevs.m_TetrahedronID = tevs.m_TetrahedronID;
        for( int i=0; i<4; i++ ) btevs.m_vecTetrahedronVID[i] = tevs.m_vecTetrahedronVID[i];
        btevs.m_InvTetrahedronBCM = tevs.m_InvTetrahedronBCM;
        btevs.m_FirstEmbeddedVertexIndex = last_embedded_vertex_index;
        btevs.m_NumEmbeddedVertices = tevs.m_vecEmbeddedVID.size();
        for( unsigned int it_ev=0; it_ev < btevs.m_NumEmbeddedVertices; it_ev++ )
            m_allocEmbeddedVID[ btevs.m_FirstEmbeddedVertexIndex + it_ev ] = tevs.m_vecEmbeddedVID[it_ev];
        last_embedded_vertex_index += btevs.m_NumEmbeddedVertices;
        last_btevs++;
    }
    GEO_ASSERT( m_NumTEVS == last_btevs );
    GEO_ASSERT( m_NumEmbeddedVID == last_embedded_vertex_index );
    SetBakedData( true, transform_e2d, m_NumTEVS, m_NumEmbeddedVID, m_allocTEVS, m_allocEmbeddedVID );

    // Free temporary data
    for( std::map< deformer_feature_index_type, TetrahedronEmbeddedVertexSubset* >::iterator it_map = map_EID2TEVS.begin();
         it_map != map_EID2TEVS.end();
         ++it_map )
        delete it_map->second;

}
void Editable_Vertices3_In_Tetrahedrons3_Barycentric_Embedding::Clear()
{
    ClearEditData();
    ClearBakedData();
}

void Editable_Vertices3_In_Tetrahedrons3_Barycentric_Embedding::ClearEditData()
{
    if(m_allocTEVS) delete[] m_allocTEVS;
    if(m_allocEmbeddedVID) delete[] m_allocEmbeddedVID;
    m_allocTEVS = 0;
    m_allocEmbeddedVID = 0;
}

} //namespace geo
