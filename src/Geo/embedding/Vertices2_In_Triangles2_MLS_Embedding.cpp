#include "Vertices2_In_Triangles2_MLS_Embedding.h"
#include <Geo/shape/shape.h> //Include all shapes
#include <map>
#include <vector>
#include <algorithm>

namespace geo {

void Vertices2_In_Triangles2_MLS_Embedding::SetBakedData( bool b_shared,
                                                          const Transform2& transform_e2d,
                                                          unsigned int num_ev, const BakedEmbeddedVertex* vec_ev )
{
    GEO_ASSERT( b_shared ); //\todo CAN ONLY BE SHARED BY NOW, NO m_pBuffer single allocation yet
    m_TransformE2D = transform_e2d;
    m_NumEV = num_ev;
    m_vecEV = vec_ev;
}

void Vertices2_In_Triangles2_MLS_Embedding::ClearBakedData()
{
    m_TransformE2D = Transform2::Identity();
    m_NumEV = 0;
    m_vecEV = 0;
}

/*!
*/
void Vertices2_In_Triangles2_MLS_Embedding::Apply( const IObject* p_deformer, IObject* p_embedded ) const
{
    //-- Apply embedding transform
    static_cast<IObject2*>(p_embedded)->SetTransform_Embedded( static_cast<const IObject2*>(p_deformer)->GetTransform() * m_TransformE2D );
    //-- Apply embedding DOF
    const Vec2* vec_deformer_sdof( reinterpret_cast<const Vec2*>( p_deformer->GetVecDOF() ) ); //GetVecDOF() updates deformer DOF if it's also embedded and out of date
    Vec2* vec_embedded_sdof( reinterpret_cast<Vec2*>( p_embedded->GetVecDOF_WriteOnly() ) ); //\note GetVecDOF_WriteOnly() just returns pointer, otherwise cyclic SyncEmbedding()->Apply()->SyncEmbedding()!
    //foreach TEVS: foreach embedded vid: compute deformed pos and write to vec_dof_embedded
    for( unsigned int it_ev=0; it_ev < m_NumEV; it_ev++ )
    {
        const BakedEmbeddedVertex& bev( m_vecEV[it_ev] );
        Vec2 vertex_pos( Vec2::Zero() );
        for( unsigned int it_influence=0; it_influence < bev.m_NumInfluences; it_influence++ )
            vertex_pos += bev.m_vecWeight[it_influence] * vec_deformer_sdof[ bev.m_vecNID[it_influence] ];
        vec_embedded_sdof[it_ev] = vertex_pos;
    }
}

//----------------------------------------------------------------
// MLS baking
//  We use doubles for precomputation to avoid precision issues in
// Inverse(M) and elsewhere (det==0,NaN...), but storing geo::Real
// in baked structures for runtime performance
//----------------------------------------------------------------
typedef double mls_real_t;
typedef mal::GVec<mls_real_t,2> mls_vec2_t;
typedef mal::GVec<mls_real_t,3> mls_vec3_t;
typedef mal::GMat<mls_real_t,3,3> mls_mat3x3_t;

//! Editing embedded vertex, influences ordered in weight decreasing order
struct EmbeddedVertex
{
    struct Influence {
        feature_index_type m_NID;
        mls_real_t m_Weight;
        inline Influence() : m_NID(cInvalidFeatureIndex), m_Weight(0) {}
        inline Influence( feature_index_type nid, mls_real_t weight ) : m_NID(nid), m_Weight(weight) {}
        //\todo No longer required, using lambda in C++11: inline bool operator<( const Influence& other ) const { return m_Weight > other.m_Weight; } //\note DECREASING order in sort()
    };
    std::vector< Influence > m_vecInfluences;
};

//! Weight Function for MLS \todo consider other?
template <typename T>
inline T WF_Spline4( T q )
{
    if( q<T(1) )
    {
        T q2( q*q );
        T q3( q2*q );
        T q4( q2*q2 );
        return T(1) - T(6)*q2 + T(8)*q3 - T(3)*q4;
    }
    else
        return T(0);
}

/* Foreach deformer node n_i
   - Compute support radius h_i
   foreach embedded Vertex v
   - Find influences n_i such that v is inside h_i around n_i
   - Precompute MLS weights
   \todo
   bake "crust" Nodes, with a single global vertex array
   consider reordering vertices
   \todo Consider using Object DOF instead of default Shape DOF
   for the embedding. This would allow user translation/rotation
   to match E and S better.
*/
void ComputeEmbedding_Polygonal2_In_MeshSolid2_As_Vertices2_In_Triangles2_MLS( const MeshSolidShape2& mss, const Transform2& tr_deformer,
                                                                               const PolygonalShape2& pss, const Transform2& tr_embedded,
                                                                               std::vector< EmbeddedVertex >& vecEV )
{
    //---- Compute support radius h_i per deformer node
    //\todo By now, use largest 1-ring distance as h_i
    const mls_real_t cFactorOneRingRadius( 1.05f );
    std::vector<mls_real_t> vec_node_h( mss.GetNumV(), 0 );
    for( unsigned int it_node=0; it_node < mss.GetNumV(); it_node++ )
    {
        mls_real_t max_dist_sq(0);
        for( geo::MeshSolidShape2::iterator_polygons_around_vertex_ccw it_ccw( mss.GetIterator_PolygonsAroundVertexCCW(it_node) );
             it_ccw.IsValid();
             ++it_ccw )
        {
            mls_real_t dist_sq( mal::NormSq( mss.V_Pos_0(it_node) - mss.V_Pos_0( mss.HE_FinalVID( it_ccw.HEID() ) ) ) );
            if( dist_sq > max_dist_sq ) max_dist_sq = dist_sq;
        }
        vec_node_h[it_node] = cFactorOneRingRadius * mal::Sqrt( max_dist_sq );
    }
    //---- For each embedded vertex find influences and compute MLS weights
    vecEV.resize( pss.GetNumVertices() );
    for( unsigned int it_v=0; it_v < pss.GetNumVertices(); it_v++ )
    {
        // Gather influences with non-zero weighting function "kernel"
        EmbeddedVertex& ev( vecEV[it_v] );
        mls_vec2_t vertex_pos( tr_embedded * pss.V_Pos_0(it_v) );
        for( unsigned int it_node=0; it_node < mss.GetNumV(); it_node++ )
        {
            mls_vec2_t node_pos( tr_deformer * mss.V_Pos_0( it_node ) );
            // Compute weighting function w_i( v )
            mls_real_t w = WF_Spline4( mal::Norm( vertex_pos - node_pos ) / vec_node_h[it_node] );
            if( w > mls_real_t(0) ) ev.m_vecInfluences.push_back( EmbeddedVertex::Influence( it_node, w ) );
        }
        // Check if there's enough influences, otherwise, M will be singular!
        //\todo To avoid this potential error we should detect it and increase some h_i radius and re-gather all influences
        GEO_ASSERT( ev.m_vecInfluences.size() > 2 );
        // Compute moment matrix M from all incluences
        // M(x) = \sum_i_1^N w( \|x-x_i\| / h_i ) p(x_i) p(x_i)^T
        mls_vec3_t P_vertex( 1, vertex_pos.x(), vertex_pos.y() );
        mls_mat3x3_t M( mls_mat3x3_t::Zero() );
        for( unsigned int it_influence=0; it_influence < ev.m_vecInfluences.size(); it_influence++ )
        {
            const EmbeddedVertex::Influence& evi( ev.m_vecInfluences[it_influence] );
            mls_vec2_t node_pos( tr_deformer * mss.V_Pos_0( evi.m_NID ) );
            mls_vec3_t P_node( 1, node_pos.x(), node_pos.y() );
            for( int i=0; i<3; i++ )
                for( int j=0; j<3; j++ )
                    M(i,j) += evi.m_Weight * P_node[i] * P_node[j];
        }
        // Compute shape-function weights for all influences from their "kernel" weights
        // \phy_i(x) = w( \|x-x_i\| / h_i ) p(x)^T M(x)^-1 p(x_i)
        mls_mat3x3_t invM( mal::Inverse(M) );
        for( unsigned int it_influence=0; it_influence < ev.m_vecInfluences.size(); it_influence++ )
        {
            EmbeddedVertex::Influence& evi( ev.m_vecInfluences[it_influence] );
            mls_vec2_t node_pos( tr_deformer * mss.V_Pos_0( evi.m_NID ) );
            mls_vec3_t P_node( 1, node_pos.x(), node_pos.y() );
            evi.m_Weight = evi.m_Weight * mal::Dot( P_vertex, invM * P_node );
        }
        //TEMP: Check Partition-Of-Unity (PoU)
        mls_real_t total_weight(0);
        for( unsigned int it_influence=0; it_influence < ev.m_vecInfluences.size(); it_influence++ ) total_weight += ev.m_vecInfluences[it_influence].m_Weight;
        if( !mal::ApproxEq( total_weight, mls_real_t(1) ) )
            GEO_LOG_ERROR( "Embedded Vertex %d with %d influences has total weight = %f", it_v, (int)ev.m_vecInfluences.size(), total_weight );
    }
}

void Editable_Vertices2_In_Triangles2_MLS_Embedding::Init( const IObject2* p_deformer, const IObject2* p_embedded )
{
    Clear();

    Transform2 transform_e2d = mal::Inverse( p_deformer->GetTransform() ) * p_embedded->GetTransform();

    // Compute embedding \todo SHAPE-SPECIFIC, could use a double switch, or a DD table, whatever...
    unsigned int num_ev = static_cast<const PolygonalShape2*>( p_embedded->GetShapeInterface() )->GetNumVertices();
    std::vector< EmbeddedVertex > vec_EV;
    ComputeEmbedding_Polygonal2_In_MeshSolid2_As_Vertices2_In_Triangles2_MLS( *static_cast<const MeshSolidShape2*>( p_deformer->GetShapeInterface() ), p_deformer->GetTransform(),
                                                                              *static_cast<const PolygonalShape2*>( p_embedded->GetShapeInterface() ), p_embedded->GetTransform(),
                                                                              vec_EV );
    // Bake embedding \todo NOT SHAPE-SPECIFIC, but embedding-type specific (V2inT2)
    m_NumEV = vec_EV.size();
    m_allocEV = new BakedEmbeddedVertex[ m_NumEV ];
    unsigned int max_influences(0);
    for( unsigned int it_ev=0; it_ev < m_NumEV; it_ev++ )
    {
        EmbeddedVertex& ev( vec_EV[it_ev] );
        BakedEmbeddedVertex& bev( m_allocEV[it_ev] );
        bev.m_NumInfluences = mal::Min<uint32>( ev.m_vecInfluences.size(), cMaxInfluencesPerVertex );
        max_influences = mal::Max<uint32>( ev.m_vecInfluences.size(), max_influences );

        // Sort and Crop influences up to a given cMaxInfluencesPerVertex
        std::sort( ev.m_vecInfluences.begin(), ev.m_vecInfluences.end(),
                   []
                   ( const EmbeddedVertex::Influence& a, const EmbeddedVertex::Influence& b )
                   { return a.m_Weight > b.m_Weight; } ); //\todo Sorts in INCREASING order
        if( bev.m_NumInfluences < ev.m_vecInfluences.size() )
            GEO_LOG_WARNING("Editable_Vertices2_In_Triangles2_MLS_Embedding::Init() Vertex %d has %d > %d influences, %d smallest influences cropped",
                            it_ev, (int)ev.m_vecInfluences.size(), Vertices2_In_Triangles2_MLS_Embedding::cMaxInfluencesPerVertex,
                            (int)ev.m_vecInfluences.size() - Vertices2_In_Triangles2_MLS_Embedding::cMaxInfluencesPerVertex );

        // Gather
        mls_real_t total_weight(0);
        for( unsigned int it_influence=0; it_influence < bev.m_NumInfluences; it_influence++ )
        {
            bev.m_vecNID[it_influence] = ev.m_vecInfluences[it_influence].m_NID;
            bev.m_vecWeight[it_influence] = ev.m_vecInfluences[it_influence].m_Weight;
            total_weight += bev.m_vecWeight[it_influence];
        }
        // Normalize
        for( unsigned int it_influence=0; it_influence < bev.m_NumInfluences; it_influence++ )
            bev.m_vecWeight[it_influence] /= total_weight;
    }
    GEO_LOG("Editable_Vertices2_In_Triangles2_MLS_Embedding::Init() max_influences = %d", max_influences );
    SetBakedData( true, transform_e2d, m_NumEV, m_allocEV );
}
void Editable_Vertices2_In_Triangles2_MLS_Embedding::Clear()
{
    ClearEditData();
    ClearBakedData();
}

void Editable_Vertices2_In_Triangles2_MLS_Embedding::ClearEditData()
{
    if(m_allocEV) delete[] m_allocEV;
    m_allocEV = 0;
}

} //namespace geo
