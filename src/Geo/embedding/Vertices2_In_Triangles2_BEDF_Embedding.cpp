#include "Vertices2_In_Triangles2_BEDF_Embedding.h"
#include <Geo/shape/shape.h> //Include all shapes
#include <map>
#include <vector>

namespace geo {

void Vertices2_In_Triangles2_BEDF_Embedding::SetBakedData( bool b_shared,
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

void Vertices2_In_Triangles2_BEDF_Embedding::ClearBakedData()
{
    m_TransformE2D = Transform2::Identity();
    m_NumTEVS = 0;
    m_NumEmbeddedVID = 0;
    m_vecTEVS = 0;
    m_vecEmbeddedVID = 0;
}

/*!
  \todo OPTIMIZATION:
    - Should manage to compute deformation field ONLY on crust nodes
      that will be actually used when deforming the embedded v_i
    - Per-node rest area (or its inverse) could be precomputed and
      stored (ONLY on crust nodes)
*/
void Vertices2_In_Triangles2_BEDF_Embedding::Apply( const IObject* p_deformer, IObject* p_embedded ) const
{
    //-- Apply embedding transform
    static_cast<IObject2*>(p_embedded)->SetTransform_Embedded( static_cast<const IObject2*>(p_deformer)->GetTransform() * m_TransformE2D );

    const Vec2* vec_deformer_sdof( reinterpret_cast<const Vec2*>( p_deformer->GetVecDOF() ) ); //GetVecDOF() updates deformer DOF if it's also embedded and out of date
    GEO_ASSERT( p_deformer->GetShapeInterface()->GetType() == eShape_MeshSolid2 );
    const MeshSolidShape2* pMSS( static_cast<const MeshSolidShape2*>(p_deformer->GetShapeInterface()) );
    const Vec2* vec_deformer_default_sdof( pMSS->GetVecDefaultSDOF() );

    //---- Compute deformation field FS(T_i) on deformer Triangles (Simplices)
    //==> USE p_embedded as MeshSolidShape2 explicitly, easier...
    Mat2x2* vec_simplex_F = new Mat2x2[ pMSS->GetNumP() ]; //\todo should alloc persistently for efficiency or use scratchpad
    for( unsigned int it_simplex=0; it_simplex<pMSS->GetNumP(); it_simplex++ )
    {
        feature_index_type vec_nid[3];
        unsigned int it_he( pMSS->P_FirstHEID(it_simplex) );
        vec_nid[0] = pMSS->HE_OriginVID(it_he); it_he = pMSS->HE_Next(it_he);
        vec_nid[1] = pMSS->HE_OriginVID(it_he); it_he = pMSS->HE_Next(it_he);
        vec_nid[2] = pMSS->HE_OriginVID(it_he); it_he = pMSS->HE_Next(it_he);
        GEO_ASSERT( it_he == pMSS->P_FirstHEID(it_simplex) ); //Closed triangle
        Mat2x2 Dm = mal::GMat2x2_From_Columns( vec_deformer_default_sdof[ vec_nid[1] ] - vec_deformer_default_sdof[ vec_nid[0] ],
                                               vec_deformer_default_sdof[ vec_nid[2] ] - vec_deformer_default_sdof[ vec_nid[0] ] );
        Mat2x2 Ds = mal::GMat2x2_From_Columns( vec_deformer_sdof[ vec_nid[1] ] - vec_deformer_sdof[ vec_nid[0] ],
                                               vec_deformer_sdof[ vec_nid[2] ] - vec_deformer_sdof[ vec_nid[0] ] );
        vec_simplex_F[it_simplex] = Ds*mal::Inverse(Dm);
    }
    //---- Apply FN(x_i) to embedded v_i using barycentric weights b_i of v_i inside containing triangle T
    //==> b_i Could be precomputed, but let's wait until final expression to decide what should we precompute/store and what should be computed at runtime
    //==> w_i Could be precomputed too, probably BETTER than precomputing b_i
    //-- Apply embedding DOF
    const Vec2* vec_embedded_default_sdof( reinterpret_cast<const Vec2*>( const_cast<const IObject*>(p_embedded)->GetShapeInterface()->GetVecDefaultDOF() ) );
    Vec2* vec_embedded_sdof( reinterpret_cast<Vec2*>( p_embedded->GetVecDOF_WriteOnly() ) ); //\note GetVecDOF_WriteOnly() just returns pointer, otherwise cyclic SyncEmbedding()->Apply()->SyncEmbedding()!
    //foreach TEVS: foreach embedded vid: compute deformed pos and write to vec_dof_embedded
    for( unsigned int it_tevs=0; it_tevs < m_NumTEVS; it_tevs++ )
    {
        const BakedTriangleEmbeddedVertexSubset& btevs( m_vecTEVS[it_tevs] );
        const feature_index_type* vec_nid( &btevs.m_vecTriangleVID[0] );
        for( unsigned int it_ev=0; it_ev < btevs.m_NumEmbeddedVertices; it_ev++ )
        {
            feature_index_type vid( m_vecEmbeddedVID[ btevs.m_FirstEmbeddedVertexIndex + it_ev ] );
            Vec3 b = btevs.m_InvTriangleBCM * mal::Concat( 1, vec_embedded_default_sdof[vid] );
            Real vec_w[4]; //w[1],w[2],w[3] correspond to neighbour simplices 0,1,2 with opposite nodes 2,0,1 respectively
            Vec2 vec_pos[4];
            Mat2x2 vec_F[4];

            // Self w0 \in [0.5..1]
            vec_w[0] = (3*3*3*3*3*3)*mal::Sq(b[0])*mal::Sq(b[1])*mal::Sq(b[2]); //3^6 * b_0^2 * b_1^2 * b_2^2 ==> [0..1]
            //vec_w[0] = (3*3*3*3*3*3*3*3*3*3*3*3)*mal::Sq(b[0]*b[0])*mal::Sq(b[1]*b[1])*mal::Sq(b[2]*b[2]); //3^9 * b_0^3 * b_1^3 * b_2^3 ==> [0..1]
            vec_w[0] = 0.5f + 0.5f*vec_w[0]; //==> [0.5..1]
            // Neighbour weights normalized to 1-w0 to guarantee PoU
            const Real cExponent(1);
            {
                // w(b_i) = 1/(1-b_i) - 1 is 0 at the boundary where b_i=1, but has dw/db(0) != 0 and dw/db(1) != 0
                vec_w[1] = mal::Rcp( mal::Pow(1-b[2],cExponent) - 1 );
                vec_w[2] = mal::Rcp( mal::Pow(1-b[0],cExponent) - 1 );
                vec_w[3] = mal::Rcp( mal::Pow(1-b[1],cExponent) - 1 );

                /* Tried to use cubic spline to guarantee
                   0-derivatives at endpoints, but non-infinity at b=0
                   breaks everything

                vec_w[1] = 2*mal::Sq(b[2])*b[2] - 3*mal::Sq(b[2]) + 1;
                vec_w[2] = 2*mal::Sq(b[2])*b[0] - 3*mal::Sq(b[0]) + 1;
                vec_w[3] = 2*mal::Sq(b[2])*b[1] - 3*mal::Sq(b[1]) + 1;
                vec_w[1] *= 1000.0f;
                vec_w[1] *= 1000.0f;
                vec_w[1] *= 1000.0f;
                vec_w[1] = mal::Sq( mal::Sq( vec_w[1] ) );
                vec_w[2] = mal::Sq( mal::Sq( vec_w[2] ) );
                vec_w[3] = mal::Sq( mal::Sq( vec_w[3] ) );
                */

                Real neighbour_weight_scaling( (1-vec_w[0]) / (vec_w[1] + vec_w[2] + vec_w[3]) ); //\todo THIS MAY BE +Inf near deformer nodes... consider clamping 1/(1-w_i) to some LARGE value instead of +Inf
                vec_w[1] *= neighbour_weight_scaling;
                vec_w[2] *= neighbour_weight_scaling;
                vec_w[3] *= neighbour_weight_scaling;

                /* Neighbour influences may NOT exist on boundary, we
                   transfer their influence to Self instead, so that
                   at a boundary edge the only non-zero influence is
                   vec_pos[0], but split into 2 parts, w[0] and
                   w[i+1]=1-w[0] for boundary edge i, other w[j] being
                   0 due to edge i infinite partial influence and edge
                   weights partition of unity.
                */
                //
                vec_pos[0] = ( vec_deformer_sdof[ vec_nid[0] ] + vec_simplex_F[ btevs.m_TriangleID ] * ( vec_embedded_default_sdof[vid] - vec_deformer_default_sdof[ vec_nid[0] ] ) );
                vec_pos[1] = ( btevs.m_vecTriangleNeighbourPID[0] != cInvalidFeatureIndex )
                             ? ( vec_deformer_sdof[ vec_nid[0] ] + vec_simplex_F[ btevs.m_vecTriangleNeighbourPID[0] ] * ( vec_embedded_default_sdof[vid] - vec_deformer_default_sdof[ vec_nid[0] ] ) )
                             : vec_pos[0];
                vec_pos[2] = ( btevs.m_vecTriangleNeighbourPID[1] != cInvalidFeatureIndex )
                             ? ( vec_deformer_sdof[ vec_nid[1] ] + vec_simplex_F[ btevs.m_vecTriangleNeighbourPID[1] ] * ( vec_embedded_default_sdof[vid] - vec_deformer_default_sdof[ vec_nid[1] ] ) )
                             : vec_pos[0];
                vec_pos[3] = ( btevs.m_vecTriangleNeighbourPID[2] != cInvalidFeatureIndex )
                             ? ( vec_deformer_sdof[ vec_nid[2] ] + vec_simplex_F[ btevs.m_vecTriangleNeighbourPID[2] ] * ( vec_embedded_default_sdof[vid] - vec_deformer_default_sdof[ vec_nid[2] ] ) )
                             : vec_pos[0];
                //
                /*\note Using barycenters has the same effect
                vec_pos[0] = ( pMSS->P_Barycenter(btevs.m_TriangleID,vec_deformer_sdof) + vec_simplex_F[ btevs.m_TriangleID ] * ( vec_embedded_default_sdof[vid] - pMSS->P_Barycenter_0(btevs.m_TriangleID) ) );
                vec_pos[1] = ( btevs.m_vecTriangleNeighbourPID[0] != cInvalidFeatureIndex )
                             ? ( pMSS->P_Barycenter(btevs.m_vecTriangleNeighbourPID[0],vec_deformer_sdof)
                                 + vec_simplex_F[ btevs.m_vecTriangleNeighbourPID[0] ] * ( vec_embedded_default_sdof[vid] - pMSS->P_Barycenter_0(btevs.m_vecTriangleNeighbourPID[0]) ) )
                             : vec_pos[0];
                vec_pos[2] = ( btevs.m_vecTriangleNeighbourPID[1] != cInvalidFeatureIndex )
                             ? ( pMSS->P_Barycenter(btevs.m_vecTriangleNeighbourPID[1],vec_deformer_sdof)
                                 + vec_simplex_F[ btevs.m_vecTriangleNeighbourPID[1] ] * ( vec_embedded_default_sdof[vid] - pMSS->P_Barycenter_0(btevs.m_vecTriangleNeighbourPID[1]) ) )
                             : vec_pos[0];
                vec_pos[3] = ( btevs.m_vecTriangleNeighbourPID[2] != cInvalidFeatureIndex )
                             ? ( pMSS->P_Barycenter(btevs.m_vecTriangleNeighbourPID[2],vec_deformer_sdof)
                                 + vec_simplex_F[ btevs.m_vecTriangleNeighbourPID[2] ] * ( vec_embedded_default_sdof[vid] - pMSS->P_Barycenter_0(btevs.m_vecTriangleNeighbourPID[2]) ) )
                             : vec_pos[0];
                */

                vec_F[0] = vec_simplex_F[ btevs.m_TriangleID ];
                /*
                vec_F[1] = vec_simplex_F[ btevs.m_TriangleID ];
                vec_F[2] = vec_simplex_F[ btevs.m_TriangleID ];
                vec_F[3] = vec_simplex_F[ btevs.m_TriangleID ];
                */
                vec_F[1] = ( btevs.m_vecTriangleNeighbourPID[0] != cInvalidFeatureIndex )
                             ? vec_simplex_F[ btevs.m_vecTriangleNeighbourPID[0] ]
                             : vec_F[0];
                vec_F[2] = ( btevs.m_vecTriangleNeighbourPID[1] != cInvalidFeatureIndex )
                             ? vec_simplex_F[ btevs.m_vecTriangleNeighbourPID[1] ]
                             : vec_F[0];
                vec_F[3] = ( btevs.m_vecTriangleNeighbourPID[2] != cInvalidFeatureIndex )
                             ? vec_simplex_F[ btevs.m_vecTriangleNeighbourPID[2] ]
                             : vec_F[0];
            }
            // Apply LBS
            vec_embedded_sdof[ vid ] = vec_w[0] * vec_pos[0]
                                       + vec_w[1] * vec_pos[1]
                                       + vec_w[2] * vec_pos[2]
                                       + vec_w[3] * vec_pos[3];
            //

            /* This fails
            Mat2x2 blended_F( vec_w[0] * vec_F[0]
                              + vec_w[1] * vec_F[1]
                              + vec_w[2] * vec_F[2]
                              + vec_w[3] * vec_F[3] );
            vec_embedded_sdof[ vid ] = vec_deformer_sdof[ vec_nid[0] ]
                                       + blended_F * ( vec_embedded_default_sdof[vid] - vec_deformer_default_sdof[ vec_nid[0] ] );
            */
        }
    }
    delete[] vec_simplex_F;
}

struct TriangleEmbeddedVertexSubset
{
    feature_index_type m_TriangleID;
    feature_index_type m_vecTriangleVID[3];
    feature_index_type m_vecTriangleNeighbourPID[3]; //\todo This is implicit in TriangleId, but requires accessing topology, it's cheaper storing ids here and processing only the vec_cage_dof/geometry
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
void ComputeEmbedding_Polygonal2_In_MeshSolid2_As_Vertices2_In_Triangles2_BEDF( const MeshSolidShape2& mss, const Transform2& tr_deformer,
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
        feature_index_type best_vecTriangleNeighbourPID[3] = { cInvalidFeatureIndex, cInvalidFeatureIndex, cInvalidFeatureIndex };
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
            Vec3 barycentric_coords( inv_BCM * mal::Concat( 1, vertex_pos ) );
            Vec3 clamped_barycentric_coords = mal::Clamp( barycentric_coords,
                                                          Vec3(0,0,0),
                                                          Vec3(1,1,1) );
            Real barycentric_distance_sq = mal::NormSq( barycentric_coords - clamped_barycentric_coords );
            if( best_TriangleID == cInvalidFeatureIndex
                || barycentric_distance_sq < best_barycentric_distance_sq )
            {
                best_TriangleID = it_p;
                best_vecTriangleVID[0] = vid0; best_vecTriangleVID[1] = vid1; best_vecTriangleVID[2] = vid2;
                unsigned int it_he( mss.P_FirstHEID(it_p) );
                best_vecTriangleNeighbourPID[0] = mss.HE_RightPID(it_he); it_he = mss.HE_Next(it_he);
                best_vecTriangleNeighbourPID[1] = mss.HE_RightPID(it_he); it_he = mss.HE_Next(it_he);
                best_vecTriangleNeighbourPID[2] = mss.HE_RightPID(it_he); it_he = mss.HE_Next(it_he);
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
            for( int i=0; i<3; i++ ) pTEVS->m_vecTriangleNeighbourPID[i] = best_vecTriangleNeighbourPID[i];
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

void Editable_Vertices2_In_Triangles2_BEDF_Embedding::Init( const IObject2* p_deformer, const IObject2* p_embedded )
{
    Clear();

    Transform2 transform_e2d = mal::Inverse( p_deformer->GetTransform() ) * p_embedded->GetTransform();

    // Compute embedding \todo SHAPE-SPECIFIC, could use a double switch, or a DD table, whatever...
    unsigned int num_embedded_vid = static_cast<const PolygonalShape2*>( p_embedded->GetShapeInterface() )->GetNumVertices();
    std::map< feature_index_type, TriangleEmbeddedVertexSubset* > map_EID2TEVS;
    ComputeEmbedding_Polygonal2_In_MeshSolid2_As_Vertices2_In_Triangles2_BEDF( *static_cast<const MeshSolidShape2*>( p_deformer->GetShapeInterface() ), p_deformer->GetTransform(),
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
        for( int i=0; i<3; i++ ) btevs.m_vecTriangleNeighbourPID[i] = tevs.m_vecTriangleNeighbourPID[i];
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
void Editable_Vertices2_In_Triangles2_BEDF_Embedding::Clear()
{
    ClearEditData();
    ClearBakedData();
}

void Editable_Vertices2_In_Triangles2_BEDF_Embedding::ClearEditData()
{
    if(m_allocTEVS) delete[] m_allocTEVS;
    if(m_allocEmbeddedVID) delete[] m_allocEmbeddedVID;
    m_allocTEVS = 0;
    m_allocEmbeddedVID = 0;
}

} //namespace geo
