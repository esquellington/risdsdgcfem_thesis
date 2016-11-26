#include "Vertices2_In_Triangles2_Barycentric_Embedding.h"
#include <Geo/shape/shape.h> //Include all shapes
#include <map>
#include <vector>

namespace geo {

void Vertices2_In_Triangles2_Barycentric_Embedding::SetBakedData( bool b_shared,
                                                                  const Transform2& transform_e2d,
                                                                  unsigned int num_tevs, unsigned num_embedded_vid,
                                                                  const BakedTriangleEmbeddedVertexSubset* vec_tevs, const feature_index_type* vec_embedded_vid )
{
    GEO_ASSERT( b_shared ); //\todo CAN ONLY BE SHARED BY NOW, NO m_pBuffer single allocation yet
    m_TransformE2D = transform_e2d;
    m_NumTEVS = num_tevs;
    m_NumEmbeddedVID = num_embedded_vid;
    m_vecTEVS = vec_tevs;
    m_vecEmbeddedVID = vec_embedded_vid;
}

void Vertices2_In_Triangles2_Barycentric_Embedding::ClearBakedData()
{
    m_TransformE2D = Transform2::Identity();
    m_NumTEVS = 0;
    m_NumEmbeddedVID = 0;
    m_vecTEVS = 0;
    m_vecEmbeddedVID = 0;
}

void Vertices2_In_Triangles2_Barycentric_Embedding::Apply( const IObject* p_deformer, IObject* p_embedded ) const
{
    //-- Apply embedding transform
    static_cast<IObject2*>(p_embedded)->SetTransform_Embedded( static_cast<const IObject2*>(p_deformer)->GetTransform() * m_TransformE2D );

    //-- Apply embedding DOF
    const Vec2* vec_deformer_sdof( reinterpret_cast<const Vec2*>( p_deformer->GetVecDOF() ) ); //GetVecDOF() updates deformer DOF if it's also embedded and out of date
    const Vec2* vec_embedded_default_sdof( reinterpret_cast<const Vec2*>( const_cast<const IObject*>(p_embedded)->GetShapeInterface()->GetVecDefaultDOF() ) );
    Vec2* vec_embedded_sdof( reinterpret_cast<Vec2*>( p_embedded->GetVecDOF_WriteOnly() ) ); //\note GetVecDOF_WriteOnly() just returns pointer, otherwise cyclic SyncEmbedding()->Apply()->SyncEmbedding()!
    //foreach TEVS: foreach embedded vid: compute deformed pos and write to vec_dof_embedded
    for( unsigned int it_tevs=0; it_tevs < m_NumTEVS; it_tevs++ )
    {
        const BakedTriangleEmbeddedVertexSubset& btevs( m_vecTEVS[it_tevs] );
        // Compute deformer matrix
        Vec2 vec_points[3];
        vec_points[0] = vec_deformer_sdof[ btevs.m_vecTriangleVID[0] ];
        vec_points[1] = vec_deformer_sdof[ btevs.m_vecTriangleVID[1] ];
        vec_points[2] = vec_deformer_sdof[ btevs.m_vecTriangleVID[2] ];
        Mat3x3 B( 1, 1, 1,
                  vec_points[0].x(), vec_points[1].x(), vec_points[2].x(),
                  vec_points[0].y(), vec_points[1].y(), vec_points[2].y() );
        /*\todo THIS IS WRONG, because m_InvTriangleBCM was
         * computed in global coords, but B is in deformer
         * local coords!!... seems to work if both Tr are Id
         * at embedding-time, but I'm afraid it WOULD FAIL if
         * any object was rotated at embedding time... we
         * should use m_TransformE2D to fix it, or, better,
         * store baked m_InvTriangleBCM in the properly
         * transformed ref sys. This problem may also apply to
         * B, write it down.
         */
        Mat3x3 M( B * btevs.m_InvTriangleBCM );
        mal::GMat<Real,2,3> m( mal::GRange<1,0,2,2>( M ) );
        for( unsigned int it_ev=0; it_ev < btevs.m_NumEmbeddedVertices; it_ev++ )
        {
            feature_index_type vid( m_vecEmbeddedVID[ btevs.m_FirstEmbeddedVertexIndex + it_ev ] );
            /*TEMP: Slow way
              Vec3 v3 = M * Vec3( 1, vec_embedded_default_sdof[vid].x(), vec_embedded_default_sdof[vid].y() );
              vec_embedded_sdof[ vid ] = Vec2( v3[1], v3[2] );
            */
            vec_embedded_sdof[ vid ] = mal::GColumn<0>(m) + mal::GRange<0,1,1,2>(m) * vec_embedded_default_sdof[vid]; // == M * [1 x y]^T, faster
        }
    }
}

struct TriangleEmbeddedVertexSubset
{
    feature_index_type m_TriangleID;
    feature_index_type m_vecTriangleVID[3];
    Mat3x3 m_InvTriangleBCM;
    std::vector<feature_index_type> m_vecEmbeddedVID;
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
void ComputeEmbedding_Polygonal2_In_MeshSolid2_As_Vertices2_In_Triangles2_Barycentric( const MeshSolidShape2& mss, const Transform2& tr_deformer,
                                                                                       const PolygonalShape2& pss, const Transform2& tr_embedded,
                                                                                       std::map< feature_index_type, TriangleEmbeddedVertexSubset* >& map_EID2TEVS )
{
    GEO_ASSERT( map_EID2TEVS.empty() );

    for( unsigned int it_v=0; it_v < pss.GetNumVertices(); it_v++ )
    {
        Vec2 vertex_pos( tr_embedded * pss.V_Pos_0(it_v) );
        // Find best embedding element e
        feature_index_type best_TriangleID( cInvalidFeatureIndex );
        feature_index_type best_vecTriangleVID[3] = { cInvalidFeatureIndex, cInvalidFeatureIndex, cInvalidFeatureIndex };
        Real best_barycentric_distance_sq( mal::Infinity<Real>() );
        Mat3x3 best_inv_BCM( Mat3x3::Identity() );
        for( unsigned int it_p=0;
             it_p < mss.GetNumP() && best_barycentric_distance_sq > 0; //stop if vertex is exactly inside an element
             it_p++ )
        {
            // Gather element nodes
            unsigned int it_he( mss.P_FirstHEID(it_p) );
            unsigned int vid0 = mss.HE_OriginVID(it_he); it_he = mss.HE_Next(it_he);
            unsigned int vid1 = mss.HE_OriginVID(it_he); it_he = mss.HE_Next(it_he);
            unsigned int vid2 = mss.HE_OriginVID(it_he); it_he = mss.HE_Next(it_he);
            GEO_ASSERT( it_he == mss.P_FirstHEID(it_p) ); //Check that it's a triangle
            Vec2 vec_element_nodes[3];
            vec_element_nodes[0] = tr_deformer * mss.V_Pos_0( vid0 );
            vec_element_nodes[1] = tr_deformer * mss.V_Pos_0( vid1 );
            vec_element_nodes[2] = tr_deformer * mss.V_Pos_0( vid2 );
            // Compute element barycentric coordinates matrix in global coords
            Mat3x3 BCM( 1, 1, 1,
                        vec_element_nodes[0].x(), vec_element_nodes[1].x(), vec_element_nodes[2].x(),
                        vec_element_nodes[0].y(), vec_element_nodes[1].y(), vec_element_nodes[2].y() );
            Mat3x3 inv_BCM( mal::Inverse( BCM ) );
            // Compute vertex barycentric coords (may be out of range)
            Vec3 barycentric_coords( inv_BCM * mal::Concat(1,vertex_pos) );
            Vec3 clamped_barycentric_coords = mal::Clamp( barycentric_coords,
                                                          Vec3(0,0,0),
                                                          Vec3(1,1,1) );
            Real barycentric_distance_sq = mal::NormSq( barycentric_coords - clamped_barycentric_coords );
            if( best_TriangleID == cInvalidFeatureIndex
                || barycentric_distance_sq < best_barycentric_distance_sq )
            {
                best_TriangleID = it_p;
                best_vecTriangleVID[0] = vid0; best_vecTriangleVID[1] = vid1; best_vecTriangleVID[2] = vid2;
                best_barycentric_distance_sq = barycentric_distance_sq;
                best_inv_BCM = inv_BCM;
            }
        }
        // Create best e TEVS if required
        if( map_EID2TEVS.find( best_TriangleID ) == map_EID2TEVS.end() )
        {
            TriangleEmbeddedVertexSubset *pTEVS = new TriangleEmbeddedVertexSubset();
            pTEVS->m_TriangleID = best_TriangleID;
            for( int i=0; i<3; i++ ) pTEVS->m_vecTriangleVID[i] = best_vecTriangleVID[i];
            pTEVS->m_InvTriangleBCM = best_inv_BCM;
            map_EID2TEVS.insert( std::make_pair( best_TriangleID, pTEVS ) );
        }
        // Add v to e
        map_EID2TEVS[best_TriangleID]->m_vecEmbeddedVID.push_back( it_v );
    }
    /*TEMP: Debug output
      std::cout << "#TEVS = " << map_EID2TEVS.size() << std::endl;
      for( std::map< feature_index_type, TriangleEmbeddedVertexSubset* >::iterator it_map = map_EID2TEVS.begin();
      it_map != map_EID2TEVS.end();
      ++it_map )
      {
      std::cout << "TEVS[" << it_map->second->m_TriangleID << "] with " << it_map->second->m_vecEmbeddedVID.size() << " vid " << std::endl;
      for( std::vector<feature_index_type>::const_iterator it_vid=it_map->second->m_vecEmbeddedVID.begin();
      it_vid != it_map->second->m_vecEmbeddedVID.end();
      ++it_vid )
      std::cout << "\tVID " << *it_vid << std::endl;
      }
    */
}

void Editable_Vertices2_In_Triangles2_Barycentric_Embedding::Init( const IObject2* p_deformer, const IObject2* p_embedded )
{
    Clear();

    Transform2 transform_e2d = mal::Inverse( p_deformer->GetTransform() ) * p_embedded->GetTransform();

    // Compute embedding \todo SHAPE-SPECIFIC, could use a double switch, or a DD table, whatever...
    unsigned int num_embedded_vid = static_cast<const PolygonalShape2*>( p_embedded->GetShapeInterface() )->GetNumVertices();
    std::map< feature_index_type, TriangleEmbeddedVertexSubset* > map_EID2TEVS;
    ComputeEmbedding_Polygonal2_In_MeshSolid2_As_Vertices2_In_Triangles2_Barycentric( *static_cast<const MeshSolidShape2*>( p_deformer->GetShapeInterface() ), p_deformer->GetTransform(),
                                                                                      *static_cast<const PolygonalShape2*>( p_embedded->GetShapeInterface() ), p_embedded->GetTransform(),
                                                                                      map_EID2TEVS );

    // Bake embedding \todo NOT SHAPE-SPECIFIC, but embedding-type specific (V2inT2)
    m_NumTEVS = map_EID2TEVS.size();
    m_NumEmbeddedVID = num_embedded_vid;
    m_allocTEVS = new BakedTriangleEmbeddedVertexSubset[ m_NumTEVS ];
    m_allocEmbeddedVID = new feature_index_type[ m_NumEmbeddedVID ];
    unsigned int last_embedded_vertex_index(0);
    unsigned int last_btevs(0);
    for( std::map< feature_index_type, TriangleEmbeddedVertexSubset* >::iterator it_map = map_EID2TEVS.begin();
         it_map != map_EID2TEVS.end();
         ++it_map )
    {
        const TriangleEmbeddedVertexSubset& tevs = *it_map->second;
        BakedTriangleEmbeddedVertexSubset& btevs = m_allocTEVS[last_btevs];
        btevs.m_TriangleID = tevs.m_TriangleID;
        for( int i=0; i<3; i++ ) btevs.m_vecTriangleVID[i] = tevs.m_vecTriangleVID[i];
        btevs.m_InvTriangleBCM = tevs.m_InvTriangleBCM;
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
    for( std::map< feature_index_type, TriangleEmbeddedVertexSubset* >::iterator it_map = map_EID2TEVS.begin();
         it_map != map_EID2TEVS.end();
         ++it_map )
        delete it_map->second;

}
void Editable_Vertices2_In_Triangles2_Barycentric_Embedding::Clear()
{
    ClearEditData();
    ClearBakedData();
}

void Editable_Vertices2_In_Triangles2_Barycentric_Embedding::ClearEditData()
{
    if(m_allocTEVS) delete[] m_allocTEVS;
    if(m_allocEmbeddedVID) delete[] m_allocEmbeddedVID;
    m_allocTEVS = 0;
    m_allocEmbeddedVID = 0;
}

} //namespace geo
