#include "TriSurfaceShape3.h"
#include <memory.h> //TEMPORAL req by memcpy()
#include "../util/GSimpleSpatialHash.h"

#ifdef __GEO_TRISS_ENABLE_BVH
#  include <boost/bind.hpp> //for ILSS
#endif

#define MIN_EDGES_TO_USE_SH 100
#define MAX_DIMENSION_SH 128
#define GEO_TRISS_EPSILON_LENGTH 1e-3 //\todo 1e-3 leaves dragon|happy 100K|HD open on Load due to excessive merging... and but armadillo100K CLP does NOT CLOSE if eps<1e-3  \todo CONSIDER using relative epsilon, an optional parameter OR even the Context epsilon...
#define GEO_TRISS_EPSILON_LENGTH_SQ mal::Sq(GEO_TRISS_EPSILON_LENGTH)

#define __ENABLE_ADDITIONAL_FIX_PASSES
#define __NUM_ADDITIONAL_FIX_PASSES 10

// float g_TTS_SliverThreshold(1e-3);
float g_TTS_SliverThreshold(25);

namespace geo
{

//-----------------------------------------------------------------------------
//---- TriSurfaceShape3 Implementation
//-----------------------------------------------------------------------------

TriSurfaceShape3::TriSurfaceShape3()
: m_Flags(eFlag_None)
, m_NumV(0), m_NumT(0)
, m_vecPoints(0), m_vecV(0), m_vecT(0)
, m_pBuffer(0)
#ifdef __GEO_TRISS_ENABLE_BVH
, m_pBVH(0)
#endif
{}

TriSurfaceShape3::~TriSurfaceShape3()
{
    ClearBakedData();
}

void TriSurfaceShape3::ClearBakedData()
{
    if( m_pBuffer ) delete [] m_pBuffer;
    m_pBuffer = 0;
    m_NumV = 0; m_vecPoints = 0; m_vecV = 0;
    m_NumT = 0; m_vecT = 0;
#ifdef __GEO_TRISS_ENABLE_BVH
    // BVH is strictly nonshared, by now
    if( m_pBVH ) delete m_pBVH;
    m_pBVH = 0;
#endif
}

void TriSurfaceShape3::ComputeBVD( bv::BoundingVolume3 &bv, const transform_type &transform, const sdof_type *vec_sdof ) const
{
    const sdof_type *actual_sdof( ( 0 != vec_sdof ) ? vec_sdof : m_vecPoints );
    switch( bv.GetType() )
    {
    case bv::eBV_Sphere3:
        GEO_ASSERT( false ); //Not yet implemented
        //bv.As<bv::Sphere2>().SetPosRadius( transform.m_Pos, m_Radius );
        break;
    case bv::eBV_AABB3:
        {
            bv::AABB3 aabb( transform*actual_sdof[0], vec_type::Zero() );
            for( unsigned int i=1; i < m_NumV; i++ )
                aabb.Merge( transform*actual_sdof[i] );
            bv.As<bv::AABB3>() = aabb;
        }
        break;
    case bv::eBV_LSS3:
        GEO_ASSERT( false ); //Not yet implemented
        //bv.As<bv::LSS2>().SetPosRadius( transform.m_Pos, transform.m_Pos, m_Radius );
        break;
    case bv::eBV_Void: break;
    case bv::eBV_Infinite: break;
    default:
        GEO_ASSERT( false ); //wrong type or dimension
        break;
    }
}

/* Init from external arrays

   If shared, arrays are NOT allocated and copied, only
   referenced. Otherwise all data allocated in m_pBuffer and copied,
   and the shape is self-contained.
*/
void TriSurfaceShape3::SetBakedData( bool b_shared,
                                     Flags32 flags,
                                     uint32 num_vertices, uint32 num_triangles,
                                     const Vec3 *vec_points, const vertex_type *vec_vertices, const triangle_type *vec_triangles )
{
    GEO_ASSERT( m_NumV == 0 && m_vecPoints == 0 && m_vecV == 0 && num_vertices > 0
                && m_NumT == 0 && m_vecT == 0 && num_triangles > 0
                && m_pBuffer == 0 );
    m_Flags = flags;
    m_NumV = num_vertices;
    m_NumT = num_triangles;
    if( b_shared )
    {
        m_vecPoints = vec_points;
        m_vecV = vec_vertices;
        m_vecT = vec_triangles;
    }
    else
    {
        // Alloc and fill single buffer
        size_t size_points_4aligned = 4*((sizeof(Vec3)*m_NumV+3)/4);
        size_t size_v_4aligned = 4*((sizeof(vertex_type)*m_NumV+3)/4);
        size_t size_t_4aligned = 4*((sizeof(triangle_type)*m_NumT+3)/4);
        size_t total_size_4aligned = size_points_4aligned + size_v_4aligned + size_t_4aligned;
        m_pBuffer = new uint32[ total_size_4aligned ];
        Vec3 *p_buffer_points = reinterpret_cast<Vec3*>( &m_pBuffer[0] );
        vertex_type *p_buffer_v = reinterpret_cast<vertex_type*>( &m_pBuffer[ size_points_4aligned ] );
        triangle_type *p_buffer_t = reinterpret_cast<triangle_type*>( &m_pBuffer[ size_points_4aligned + size_v_4aligned ] );
        memcpy( p_buffer_points, vec_points, sizeof(Vec3)*m_NumV );
        memcpy( p_buffer_v, vec_vertices, sizeof(vertex_type)*m_NumV );
        memcpy( p_buffer_t, vec_triangles, sizeof(triangle_type)*m_NumT );
        // Save const pointers
        m_vecPoints = p_buffer_points;
        m_vecV = p_buffer_v;
        m_vecT = p_buffer_t;
    }
}

Vec3 TriSurfaceShape3::T_Barycenter_0( uint32 tid ) const
{
    return mal::Rcp<Real>(3) * ( V_Pos_0( m_vecT[tid].m_vecVID[0] ) +
                                 V_Pos_0( m_vecT[tid].m_vecVID[1] ) +
                                 V_Pos_0( m_vecT[tid].m_vecVID[2] ) );
}

Vec3 TriSurfaceShape3::T_Barycenter( uint32 tid, const sdof_type *vec_sdof ) const
{
    return mal::Rcp<Real>(3) * ( vec_sdof[ m_vecT[tid].m_vecVID[0] ] +
                                 vec_sdof[ m_vecT[tid].m_vecVID[1] ] +
                                 vec_sdof[ m_vecT[tid].m_vecVID[2] ] );
}

void TriSurfaceShape3::T_Edge_0( uint32 tid, uint32 eit, Vec3& p0, Vec3& p1 ) const
{
    GEO_LOG_ASSERT( tid < m_NumT && eit < 3, "TriSurfaceShape3(%d,%d) out of range!", tid, eit );
    p0 = V_Pos_0( T_VID(tid,eit) );
    p1 = V_Pos_0( T_VID(tid,(eit+1)%3) );
}

//-----------------------------------------------------------------------------
//---- EditableTriSurfaceShape3 Implementation
//-----------------------------------------------------------------------------
EditableTriSurfaceShape3::EditableTriSurfaceShape3()
: m_IsBeingEdited(false)
, m_IsClosed(false), m_IsConnected(false)
{
}

EditableTriSurfaceShape3::~EditableTriSurfaceShape3()
{
    Clear();
}

void EditableTriSurfaceShape3::Clear()
{
    GEO_ASSERT( !IsBeingEdited() );
    ClearEditData();
    TriSurfaceShape3::ClearBakedData();
}

void EditableTriSurfaceShape3::Set( const TriSurfaceShape3 &tss3 )
{
    GEO_ASSERT( !IsBeingEdited() );
    Clear();
    BeginEdition();
    // Copy other mesh baked data
    for( unsigned int i=0; i<tss3.GetNumV(); i++ )
        m_addV.push_back( editable_vertex_type( tss3.GetVecPoints()[i] ) );
    for( unsigned int i=0; i<tss3.GetNumT(); i++ )
        m_addT.push_back( editable_triangle_type( tss3.GetVecT()[i].m_vecVID[0],
                                                  tss3.GetVecT()[i].m_vecVID[1],
                                                  tss3.GetVecT()[i].m_vecVID[2] ) );
    EndEdition();
}

void EditableTriSurfaceShape3::BeginEdition()
{
    GEO_ASSERT( !IsBeingEdited() );
    // Clear anything added outside Begin/End
    ClearEditData();
    // Copy baked data
    for( unsigned int i=0; i<m_NumV; i++ )
        m_addV.push_back( editable_vertex_type( m_vecPoints[i] ) );
    for( unsigned int i=0; i<m_NumT; i++ )
        m_addT.push_back( editable_triangle_type( m_vecT[i].m_vecVID[0],
                                                  m_vecT[i].m_vecVID[1],
                                                  m_vecT[i].m_vecVID[2] ) );
    // Delete baked data
    TriSurfaceShape3::ClearBakedData();
    m_IsBeingEdited = true;
}

bool EditableTriSurfaceShape3::EndEdition()
{
    if( IsBeingEdited() )
    {
        RebuildAABB();
#ifdef __ENABLE_ADDITIONAL_FIX_PASSES
        bool bRepeat(false);
        unsigned int num_passes(0);
        do
        {
            bool bFD = FixDegeneracies();
            bool bRT = RebuildTopology();
            bRepeat = (bFD || bRT);
        } while( bRepeat && num_passes++ < __NUM_ADDITIONAL_FIX_PASSES );
        // //TEMPORAL:
        // if( RemoveIsolatedTriangles() )
        // {
        //     // This should fix bug in Clip_TriSS_TetSS()...sometimes,
        //     // but MAY break other stuff if done in early Fix/Rebuild
        //     // iterations... we do a SINGLE post-pass
        //     bool bFD = FixDegeneracies();
        //     bool bRT = RebuildTopology();
        //     bRepeat = (bFD || bRT);
        // }
        //TEMPORAL
        if( bRepeat ) { GEO_LOG_WARNING("EditableTriSurfaceShape3::EndEdition() FixDegeneracies/RebuildTopology additional passes required but not performed!"); }
        else if( num_passes > 1 ) { GEO_LOG("EditableTriSurfaceShape3::EndEdition() FixDegeneracies/RebuildTopology with %d additional passes", num_passes-1 ); }
#else
        FixDegeneracies();
        RebuildTopology();
#endif
        if( m_IsClosed && m_addV.size() > 3 && m_addT.size() > 3 ) RebuildClosed();
        RebuildBakedData();
        ClearEditData();
        m_IsBeingEdited = false;
        return true;
    }
    else
        return false;
}

trisurface3_feature_index_type EditableTriSurfaceShape3::AddVertex( const Vec3 &point )
{
    GEO_ASSERT( IsBeingEdited() );
    GEO_ASSERT( m_addV.size() < cTriSurface3_InvalidFeatureIndex-1 );
    m_addV.push_back( editable_vertex_type( point ) );
    return trisurface3_feature_index_type(m_addV.size()-1);
}

trisurface3_feature_index_type EditableTriSurfaceShape3::AddTriangle( trisurface3_feature_index_type vid0, trisurface3_feature_index_type vid1, trisurface3_feature_index_type vid2 )
{
    GEO_ASSERT( IsBeingEdited() );
    //GEO_LOG_WARNING("EditableTriSurfaceShape3::AddTriangle(%d,%d,%d) start...", vid0, vid1, vid2 );
    //\todo Consider looking for existing T
    GEO_ASSERT( m_addT.size() < cTriSurface3_InvalidFeatureIndex-1 );
    GEO_ASSERT( vid0 < m_addV.size() && vid1 < m_addV.size() && vid2 < m_addV.size()  );
    // Alloc new TID
    trisurface3_feature_index_type tid( m_addT.size() );
    // Create Triangle
    m_addT.push_back( editable_triangle_type( vid0, vid1, vid2 ) );
    return tid;
}

void EditableTriSurfaceShape3::Transform( const Transform3& tr )
{
    GEO_ASSERT( IsBeingEdited() );
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
        m_addV[it_v].m_Pos = tr * m_addV[it_v].m_Pos;
}

void EditableTriSurfaceShape3::Scale( Real scale )
{
    GEO_ASSERT( IsBeingEdited() );
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
        m_addV[it_v].m_Pos *= scale;
}

/* Given a non-closed but connected TSS, try to fix holes by:
   - If __ENABLE_TRYTOCLOSEHOLES_MULTI Snap all edges below
     epsilon_length threshold
   - ELSE: Collapse the single shortest open edge, expecting it to
     "fix" a single hole.
   - \todo MAYBE we could close each hole explicitly, finding "good"
     triplets of open edges and generating a new tri between them, but
     real holes may be much more complex than a single missing
     triangle, therefore we just collapse open edges and expect
     EndEdition->FixDegeneracies to do the rest.  only fix simple
     holes...
*/
bool EditableTriSurfaceShape3::TryToCloseHoles( Real epsilon_length )
{
    GEO_ASSERT( IsBeingEdited() );
    //IMPORTANT: We NEED up to date topology here to identify open edges correctly
    RebuildTopology();
    if( m_IsClosed ) return false;
#define __ENABLE_TRYTOCLOSEHOLES_NEW
#ifdef __ENABLE_TRYTOCLOSEHOLES_NEW
    // Gather all open edges below epsilon_length, and also the shortest
    Real epsilon_length_sq( mal::Sq(epsilon_length) );
    Real shortest_open_edge_length_sq(mal::Infinity<Real>());
    trisurface3_edge_id_type shortest_open_eid(cTriSurface3_InvalidFeatureIndex,-1);
    std::vector<trisurface3_edge_id_type> vec_open_edges;
    for( unsigned int it_tri=0; it_tri<m_addT.size(); it_tri++ )
        for( unsigned int it_vit=0; it_vit<3; it_vit++ )
            if( m_addT[it_tri].m_vecNTID[it_vit] == cTriSurface3_InvalidFeatureIndex )
            {
                Real length_sq( mal::NormSq( m_addV[m_addT[it_tri].m_vecVID[it_vit]].m_Pos - m_addV[m_addT[it_tri].m_vecVID[(it_vit+1)%3]].m_Pos ) );
                if( length_sq < shortest_open_edge_length_sq )
                {
                    shortest_open_edge_length_sq = length_sq;
                    shortest_open_eid = trisurface3_edge_id_type(it_tri,it_vit);
                }
                if( length_sq < epsilon_length_sq )
                    vec_open_edges.push_back( trisurface3_edge_id_type(it_tri,it_vit) );
            }
    // Collapse all below epsilon_length
    if( !vec_open_edges.empty() )
    {
        for( auto oeid : vec_open_edges )
        {
            uint32 vid0( m_addT[oeid.first].m_vecVID[oeid.second] );
            uint32 vid1( m_addT[oeid.first].m_vecVID[(oeid.second+1)%3] );
            Vec3& p0( m_addV[ vid0 ].m_Pos );
            Vec3& p1( m_addV[ vid1 ].m_Pos );
            Vec3 snapped_p0p1 = Real(0.5)*(p0+p1);
            p0 = p1 = snapped_p0p1;
        }
        return true;
    }
    else if( shortest_open_eid.first != cTriSurface3_InvalidFeatureIndex )
    {
        uint32 vid0( m_addT[shortest_open_eid.first].m_vecVID[shortest_open_eid.second] );
        uint32 vid1( m_addT[shortest_open_eid.first].m_vecVID[(shortest_open_eid.second+1)%3] );
        Vec3& p0( m_addV[ vid0 ].m_Pos );
        Vec3& p1( m_addV[ vid1 ].m_Pos );
        Vec3 snapped_p0p1 = Real(0.5)*(p0+p1);
        p0 = p1 = snapped_p0p1;
        //GEO_LOG_WARNING("SNAPPING edge (%d,%d) with VID (%d,%d)", shortest_open_eid.first, shortest_open_eid.second, vid0, vid1 );
        return true;
    }
    else
        return false;
#else
    if( epsilon_length > 0 ) // Remove all open edges below epsilon_length
    {
        // Gather all open edges
        std::vector<trisurface3_edge_id_type> vec_open_edges;
        for( unsigned int it_tri=0; it_tri<m_addT.size(); it_tri++ )
            for( unsigned int it_neighbour=0; it_neighbour<3; it_neighbour++ )
                if( m_addT[it_tri].m_vecNTID[it_neighbour] == cTriSurface3_InvalidFeatureIndex )
                    vec_open_edges.push_back( trisurface3_edge_id_type(it_tri,it_neighbour) );
        // Collapse
        Real epsilon_length_sq( mal::Sq(epsilon_length) );
        unsigned int num_collapsed(0);
        for( auto oeid : vec_open_edges )
        {
            uint32 vid0( m_addT[oeid.first].m_vecVID[oeid.second] );
            uint32 vid1( m_addT[oeid.first].m_vecVID[(oeid.second+1)%3] );
            Vec3& p0( m_addV[ vid0 ].m_Pos );
            Vec3& p1( m_addV[ vid1 ].m_Pos );
            if( mal::NormSq(p0-p1) < epsilon_length_sq )
            {
                Vec3 snapped_p0p1 = Real(0.5)*(p0+p1);
                p0 = p1 = snapped_p0p1;
                num_collapsed++;
                //GEO_LOG_WARNING("SNAPPING edge (%d,%d) with VID (%d,%d)", shortest_open_eid.first, shortest_open_eid.second, vid0, vid1 );
            }
        }
        return num_collapsed > 0;
    }
    else // Remove shortest open edge
    {
        // Find shortest open edge
        Real shortest_open_edge_length_sq(mal::Infinity<Real>());
        trisurface3_edge_id_type shortest_open_eid(cTriSurface3_InvalidFeatureIndex,-1);
        for( unsigned int it_tri=0; it_tri<m_addT.size(); it_tri++ )
            for( unsigned int it_neighbour=0; it_neighbour<3; it_neighbour++ )
                if( m_addT[it_tri].m_vecNTID[it_neighbour] == cTriSurface3_InvalidFeatureIndex )
                {
                    Real length_sq( mal::NormSq( m_addV[m_addT[it_tri].m_vecVID[it_neighbour]].m_Pos - m_addV[m_addT[it_tri].m_vecVID[(it_neighbour+1)%3]].m_Pos ) );
                    if( length_sq < shortest_open_edge_length_sq )
                    {
                        shortest_open_edge_length_sq = length_sq;
                        shortest_open_eid = std::make_pair(it_tri,it_neighbour);
                    }
                }
        // Collapse
        if( shortest_open_eid.first != cTriSurface3_InvalidFeatureIndex )
        {
            uint32 vid0( m_addT[shortest_open_eid.first].m_vecVID[shortest_open_eid.second] );
            uint32 vid1( m_addT[shortest_open_eid.first].m_vecVID[(shortest_open_eid.second+1)%3] );
            Vec3& p0( m_addV[ vid0 ].m_Pos );
            Vec3& p1( m_addV[ vid1 ].m_Pos );
            Vec3 snapped_p0p1 = Real(0.5)*(p0+p1);
            p0 = p1 = snapped_p0p1;
            //GEO_LOG_WARNING("SNAPPING edge (%d,%d) with VID (%d,%d)", shortest_open_eid.first, shortest_open_eid.second, vid0, vid1 );
            return true;
        }
        else
            return false;
    }
#endif
}

/* Dumb simplification, collapse all edges below epsilon_length, no
   error criterion, just a quick'n'dirty test...
*/
unsigned int EditableTriSurfaceShape3::Simplify_EdgeCollapse( Real epsilon_length )
{
    GEO_ASSERT( IsBeingEdited() );
    Real epsilon_length_sq(mal::Sq(epsilon_length));
    //IMPORTANT: We NEED up to date topology here to identify open edges correctly
    RebuildTopology();
    unsigned int num_collapsed(0);
    for( unsigned int it_tri=0; it_tri<m_addT.size(); it_tri++ )
        for( unsigned int it_neighbour=0; it_neighbour<3; it_neighbour++ )
        {
            Vec3& p0( m_addV[ m_addT[it_tri].m_vecVID[it_neighbour] ].m_Pos );
            Vec3& p1( m_addV[m_addT[it_tri].m_vecVID[(it_neighbour+1)%3]].m_Pos );
            if( mal::NormSq( p0 - p1 ) < epsilon_length_sq )
            {
                p0 = p1 = Real(0.5)*(p0+p1);
                num_collapsed++;
            }
        }
    return num_collapsed;
}

#ifdef __GEO_TRISS_ENABLE_BVH
bool EditableTriSurfaceShape3::AddBVH()
{
    if( m_pBVH ) delete m_pBVH;
    m_pBVH = new BVH_TriSurfaceShape3;
    uint32 num_entries( m_NumT ); //\todo uint32 num_entries( m_pDCR ? m_pDCR->m_NumElements : m_NumP ); //TEMP: If there's a DCR, only use those P in the BVH
    m_pBVH->Rebuild_TopDown( num_entries,
                             boost::bind<void>( &GEBV_TriSurfaceShape3_E<BVH_TriSurfaceShape3::entry_index_type,BVH_TriSurfaceShape3::bv_type>,
                                                this, transform_type::Identity(), GetVecDefaultSDOF(),
                                                _1, _2) );
    /*
    m_pBVH->Rebuild( num_entries,
                     boost::bind<void>( &GEBV_TriSurfaceShape3_E<BVH_TriSurfaceShape3::entry_index_type,BVH_TriSurfaceShape3::bv_type>,
                                        this, transform_type::Identity(), GetVecDefaultSDOF(),
                                        _1, _2) );
    */
    return true;
}
#endif

//---- Internal methods
void EditableTriSurfaceShape3::ClearEditData()
{
    m_addV.clear();
    m_addT.clear();
}

void EditableTriSurfaceShape3::RebuildBakedData()
{
    // Bake Vertices
    Vec3 *vecPoints( new Vec3[m_addV.size()] );
    vertex_type *vecV( new vertex_type[m_addV.size()] );
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
    {
        vecPoints[it_v] = m_addV[it_v].m_Pos;
        vecV[it_v] = m_addV[it_v];
    }
    // Bake Triangles
    triangle_type *vecT( new triangle_type[m_addT.size()] );
    for( unsigned int it_t=0; it_t<m_addT.size(); it_t++ )
        vecT[it_t] = m_addT[it_t];
    // Set baked data as non-shared
    SetBakedData( false,
                  (m_IsClosed?eFlag_Closed:0) | (m_IsConnected?eFlag_Connected:0),
                  m_addV.size(), m_addT.size(),
                  vecPoints, vecV, vecT );

    // Clear temporal stuff
    delete [] vecPoints;
    delete [] vecV;
    delete [] vecT;
}

void EditableTriSurfaceShape3::RebuildAABB()
{
    m_AABB = bv::AABB3( m_addV[0].m_Pos );
    for( unsigned int it_v=1; it_v<m_addV.size(); it_v++ ) m_AABB.Merge( m_addV[it_v].m_Pos );
}

bool EditableTriSurfaceShape3::FixDegeneracies()
{
    bool bFixDV = FixDegenerateVertices();
    bool bFixDT = FixDegenerateTriangles();
    return bFixDV || bFixDT;
}

bool EditableTriSurfaceShape3::RebuildTopology()
{
    if( m_addV.size() == 0 || m_addT.size() == 0 ) return false;

    // Rebuild triangle adjacency
    typedef std::pair<trisurface3_feature_index_type,int> trisurface3_edge_id_type;

    uint32 num_edges = 3*m_addT.size();
    std::vector< std::pair<trisurface3_edge_id_type,trisurface3_edge_id_type> > vec_topology_unmatched_edge_pairs;
    // Use SH only for non-trivial meshes, otherwise the O(n) complexity is not worth the grid overhead
    //TEMP: Do it always, by now, shorter code... if( num_edges > MIN_EDGES_TO_USE_SH )
    {
        // Create SH (AABB assumed up-to-date from EndFaces()->Rebuild())
        uint32 dimensions[3] = { MAX_DIMENSION_SH, MAX_DIMENSION_SH, MAX_DIMENSION_SH };
        // Compute SH resolution for the largest extent
        Vec3 sizes( m_AABB.GetMax()-m_AABB.GetMin() );
        Real resolution = mal::Max(sizes) / MAX_DIMENSION_SH;
        for( int i=0; i<3; i++ )
            dimensions[i] = mal::Clamp<uint32>( mal::IntPart<uint32>( mal::Ceil( sizes[i] / resolution ) ), 1, MAX_DIMENSION_SH );
        //IMPORTANT, We force uint32 entry_index_type in GSimpleSpatialHash because even if trisurface3_feature_index_type = uint16, it must allocate 3*T edge-entries
        GSimpleSpatialHash<uint32, trisurface3_edge_id_type> ssh( num_edges,
                                                                  dimensions,
                                                                  m_AABB.GetMin(), m_AABB.GetMax() );
        // GEO_LOG( "EditableTriSurfaceShape3::RebuildTopology() SH dim = (%d,%d,%d)", dimensions[0], dimensions[1], dimensions[2] );

        // First pass: classify edge midpoints AND CLEAR TOPOLOGY
        ssh.BeginClassify();
        {
            for( unsigned int it_tri=0; it_tri<m_addT.size(); it_tri++ )
            {
                editable_triangle_type& tri( m_addT[it_tri] );
                for( unsigned int it_neighbour=0; it_neighbour<3; it_neighbour++ )
                    ssh.Classify( std::make_pair( it_tri, it_neighbour ),
                                  Real(0.5)*( m_addV[ tri.m_vecVID[it_neighbour] ].m_Pos
                                              + m_addV[ tri.m_vecVID[(it_neighbour+1)%3] ].m_Pos ) );
                tri.ResetNeigbours();
            }
        }
        ssh.EndClassify();

        // Second pass, test overlaps
        /* For each edge, we test its midpoint and test all potential overlaps for the opposite edge.
           - Opposite half-edges have EXACTLY the same midpoint, which ensures classification in the SAME cell
           - However, we DO NOT compare midpoints, we explicitly compare edge VID to identify opposite half-edges.
        */
        uint32 num_tests_sh(0);
        ssh.BeginTestAdd();
        {
            std::vector<trisurface3_edge_id_type> vec_potential_overlaps;
            for( unsigned int it_tri=0; it_tri<m_addT.size(); it_tri++ )
            {
                editable_triangle_type& tri( m_addT[it_tri] );
                for( unsigned int it_neighbour=0; it_neighbour<3; it_neighbour++ )
                {
                    // if open, find potential matching
                    if( tri.m_vecNTID[it_neighbour] == cTriSurface3_InvalidFeatureIndex //\todo THIS requires ResetNeigbours(), otherwise previous passe topology is STILL available!
                        &&
                        ssh.TestAdd( std::make_pair( it_tri, it_neighbour ),
                                     Real(0.5)*( m_addV[ tri.m_vecVID[it_neighbour] ].m_Pos
                                                 + m_addV[ tri.m_vecVID[(it_neighbour+1)%3] ].m_Pos ),
                                     GEO_TRISS_EPSILON_LENGTH, //GEO_TRISS_EPSILON_LENGTH_SQ, //\todo CRITICAL: this SHOULD BE GEO_TRISS_EPSILON_LENGTH, not squared... MANY topology problems could derive from this shit!!
                                     vec_potential_overlaps ) )
                    {
                        // Given edge (it_tri,it_neighbour) with VID
                        // (vid0,vid1), check if any potential overlap
                        // is the opposite edge (vid1,vid0)
                        //\note That (it_tri,it_neighbour) is *NEVER* in vec_potential_overlaps (as it's incrementally added in Test())
                        for( uint32 it_pov=0; it_pov < vec_potential_overlaps.size(); it_pov++ )
                        {
                            trisurface3_edge_id_type neid( vec_potential_overlaps[it_pov] );
                            // Check not already closed candidate (otherwise, non-2-manifold!)
                            bool bOpen( m_addT[neid.first].m_vecNTID[neid.second] == cTriSurface3_InvalidFeatureIndex );
                            // Check not on same triangle (otherwise, 2 edges on a degenerate triangle!)
                            bool bDifferentTriangle( it_tri != neid.first );
                            // Topology-matching is conclusive but requires coincident vertices being already merged (same VID)
                            bool bMatchingTopology( tri.m_vecVID[it_neighbour] == m_addT[neid.first].m_vecVID[ (neid.second+1)%3 ]
                                                    && tri.m_vecVID[(it_neighbour+1)%3] == m_addT[neid.first].m_vecVID[ neid.second ] );
                            // Geometry-matching is nonconclusive
                            bool bMatchingGeometry( mal::NormSq( Real(0.5)*( m_addV[ tri.m_vecVID[it_neighbour] ].m_Pos + m_addV[ tri.m_vecVID[(it_neighbour+1)%3] ].m_Pos )
                                                                 - Real(0.5)*( m_addV[ m_addT[neid.first].m_vecVID[neid.second] ].m_Pos + m_addV[ m_addT[neid.first].m_vecVID[(neid.second+1)%3] ].m_Pos ) )
                                                    < GEO_TRISS_EPSILON_LENGTH_SQ );
                            if( bMatchingTopology )
                            {
                                if( bOpen && bDifferentTriangle )
                                {
                                    tri.m_vecNTID[it_neighbour] = neid.first;
                                    m_addT[neid.first].m_vecNTID[neid.second] = it_tri;
                                }
                                /*TEMP: VERBOSE, enable to debug
                                else if( !bOpen )
                                {
                                    GEO_LOG_WARNING("EditableTriSurfaceShape3::RebuildTopology() open edge (%d,%d) topology-matches closed edge (%d,%d) connected to triangle %d",
                                                    it_tri, it_neighbour, neid.first, neid.second, m_addT[neid.first].m_vecNTID[neid.second] );
                                }
                                else //!bDifferentTriangle
                                {
                                    GEO_LOG_WARNING("EditableTriSurfaceShape3::RebuildTopology() open edge (%d,%d) topology-matches open edge (%d,%d) on the same triangle %d",
                                                    it_tri, it_neighbour, neid.first, neid.second, it_tri );
                                }
                                */
                            }
                            else if( bMatchingGeometry )
                            {
                                /* Geometry-based matching is NOT
                                   CONCLUSIVE, but may represent edge
                                   pairs with one/both endpoints VERY
                                   close but not necesarily merged
                                   under GEO_TRISS_EPSILON_LENGTH_SQ
                                   criterion. we do NOT have enough
                                   information here to accept/discard
                                   the matching. We store the edge
                                   pair as unmatched and, if NO OTHER
                                   matches are found for any of them,
                                   we'll accept their matching by
                                   means of vertex snapping followed
                                   an additional
                                   FixDegenerateVertices() and
                                   RebuildTopology() pass.
                                */
                                if( bOpen && bDifferentTriangle )
                                {
                                    vec_topology_unmatched_edge_pairs.push_back( std::make_pair( std::make_pair(it_tri,it_neighbour), neid ) );
                                }
                                /*\note This is NOISY... and already summarized in vec_topology_unmatched_edge_pairs post-process
                                else if( !bOpen )
                                {
                                    GEO_LOG_WARNING("EditableTriSurfaceShape3::RebuildTopology() open edge (%d,%d) geometry-matches closed edge (%d,%d) connected to triangle %d",
                                                    it_tri, it_neighbour, neid.first, neid.second, m_addT[neid.first].m_vecNTID[neid.second] );
                                }
                                else //!bDifferentTriangle
                                {
                                    GEO_LOG_WARNING("EditableTriSurfaceShape3::RebuildTopology() open edge (%d,%d) geometry-matches open edge (%d,%d) on the same triangle %d",
                                                    it_tri, it_neighbour, neid.first, neid.second, it_tri );
                                }
                                */
                            }
                            // else Neither topology nor geometry matching, discard
                            num_tests_sh++;
                        }
                    }
                }
            }
        }
        ssh.EndTestAdd();
        // GEO_LOG( "EditableTriSurfaceShape3::RebuildTopology() SH tests = %d << %d", num_tests_sh, (num_edges * (num_edges-1)) / 2 );
    }

    // For bMatchingGeometry && !bMatchingTopology edge pairs, if BOTH still open, SNAL their vertices geometrically and AND require a new FixDegeneracies()->RebuildTopology() iteration
    unsigned int num_open_topology_unmatched_edge_pairs(0);
    if( !vec_topology_unmatched_edge_pairs.empty() )
    {
        // GEO_LOG_WARNING("EditableTriSurfaceShape3::RebuildTopology() found %d topology-unmatched edge pairs", (int)vec_topology_unmatched_edge_pairs.size() );
        for( auto it_edge_pair : vec_topology_unmatched_edge_pairs )
        {
            if( m_addT[it_edge_pair.first.first].m_vecNTID[it_edge_pair.first.second] == cTriSurface3_InvalidFeatureIndex
                && m_addT[it_edge_pair.second.first].m_vecNTID[it_edge_pair.second.second] == cTriSurface3_InvalidFeatureIndex )
            {
                Vec3& p0( m_addV[ m_addT[it_edge_pair.first.first].m_vecVID[it_edge_pair.first.second] ].m_Pos );
                Vec3& p1( m_addV[ m_addT[it_edge_pair.first.first].m_vecVID[(it_edge_pair.first.second+1)%3] ].m_Pos );
                Vec3& q0( m_addV[ m_addT[it_edge_pair.second.first].m_vecVID[it_edge_pair.second.second] ].m_Pos );
                Vec3& q1( m_addV[ m_addT[it_edge_pair.second.first].m_vecVID[(it_edge_pair.second.second+1)%3] ].m_Pos );
                Vec3 snapped_p0q1 = Real(0.5)*(p0+q1);
                Vec3 snapped_p1q0 = Real(0.5)*(p1+q0);
                p0 = q1 = snapped_p0q1;
                p1 = q0 = snapped_p1q0;
                num_open_topology_unmatched_edge_pairs++;
            }
        }
    }
    if( num_open_topology_unmatched_edge_pairs > 0 )
        GEO_LOG_WARNING("EditableTriSurfaceShape3::RebuildTopology() additional pass required for %u (snapped) open topology-unmatched edge pairs",
                        num_open_topology_unmatched_edge_pairs );

    // Compute topology stuff (Closed, \todo Connected,\todo SimpleConnected...)
    unsigned int num_open_edges(0);
    unsigned int num_nonmanifold_edges(0);
    for( unsigned int it_tri=0; it_tri<m_addT.size(); it_tri++ )
    {
        for( unsigned int it_neighbour=0; it_neighbour<3; it_neighbour++ )
        {
            uint32 ntid( m_addT[it_tri].m_vecNTID[it_neighbour] );
            bool bOpen( ntid == cTriSurface3_InvalidFeatureIndex );
            bool bNonManifold( !bOpen &&
                               (m_addT[ntid].m_vecNTID[0] != it_tri
                                && m_addT[ntid].m_vecNTID[1] != it_tri
                                && m_addT[ntid].m_vecNTID[2] != it_tri) );
            if( bOpen )
            {
                num_open_edges++;
                /*TEMP: Debugging open edges...
                uint32 vid0( m_addT[it_tri].m_vecVID[it_neighbour] );
                uint32 vid1( m_addT[it_tri].m_vecVID[(it_neighbour+1)%3] );
                Vec3 p0( m_addV[ vid0 ].m_Pos );
                Vec3 p1( m_addV[ vid1 ].m_Pos );
                GEO_LOG_WARNING( "Open edge (%d,%d) VID (%d,%d) POS (%f,%f,%f)->(%f,%f,%f)",
                                 it_tri, it_neighbour,
                                 vid0, vid1,
                                 p0[0], p0[1], p0[2],
                                 p1[0], p1[1], p1[2] );
                */
            }
            else if( bNonManifold )
            {
                num_nonmanifold_edges++;
                /*TEMP: Debugging nonmanifold edges...
                uint32 vid0( m_addT[it_tri].m_vecVID[it_neighbour] );
                uint32 vid1( m_addT[it_tri].m_vecVID[(it_neighbour+1)%3] );
                GEO_LOG_WARNING( "Non-2-manifold edge (%d,%d) VID (%d,%d) has opposite NTID (%d,%d,%d)",
                                 it_tri, it_neighbour,
                                 vid0, vid1,
                                 m_addT[ntid].m_vecNTID[0], m_addT[ntid].m_vecNTID[1], m_addT[ntid].m_vecNTID[2] );
                */
            }
        }
    }
    //\todo We could try to REPAIR open/nonmanifold edges here... sometimes one causes the other, and a nonmanifold edge could become closed by stealing its adjacency from another one...
    m_IsClosed = num_open_edges == 0;
    // GEO_LOG("EditableTriSurfaceShape3::IsClosed() = %s (%u open, %u non-manifold)", m_IsClosed?"true":"false", num_open_edges, num_nonmanifold_edges );

    // Compute IsConnected
    uint32 num_flooded(0);
    std::vector<bool> vec_tri_flooded( m_addT.size(), false );
    std::vector<uint32> stackTID;
    stackTID.push_back(0);
    while( !stackTID.empty() )
    {
        uint32 tid( stackTID.back() );
        stackTID.pop_back();
        const editable_triangle_type& tri( m_addT[tid] );
        for( unsigned int it_neighbour=0; it_neighbour<3; it_neighbour++ )
        {
            uint32 ntid( tri.m_vecNTID[it_neighbour] );
            if( ntid != cTriSurface3_InvalidFeatureIndex )
            {
                if( !vec_tri_flooded[ntid] )
                {
                    num_flooded++;
                    vec_tri_flooded[ntid] = true;
                    stackTID.push_back( ntid );
                }
            }
        }
    }
    m_IsConnected = num_flooded == m_addT.size();
    // GEO_LOG("EditableTriSurfaceShape3::IsConnected() = %s (%d unflooded from tid 0)", m_IsConnected?"true":"false", (int)m_addT.size()-num_flooded );

    //\todo Rebuild vertex topology

    return num_open_topology_unmatched_edge_pairs > 0;
}

/*Closed-specific rebuild
  - Enforce CCW face orientation
    - Test ray from potentially-inside point ray_pos to
      guaranteed-outside point ray_pos+ray_length*ray_dir, and, if p0
      crosses surface n=2k+1 times (odd), it's ACTUALLY inside and the
      surface is correctly oriented, if not, we must reverse all
      triangle orientations.
      - \todo We only chech Tri[0] orientation and assume all T have the
        same CW/CCW now. We could perform this test per-triangle, if
        required (however, if tri orientation was not coherent, the
        trisurface would not have been classified as Closed...)
*/
void EditableTriSurfaceShape3::RebuildClosed()
{
    // Initialize ray
    uint32 vid(0);
    Vec3 vec_pos[3] = { m_addV[ m_addT[vid].m_vecVID[0] ].m_Pos, m_addV[ m_addT[vid].m_vecVID[1] ].m_Pos, m_addV[ m_addT[vid].m_vecVID[2] ].m_Pos };
    Vec3 normal = mal::Normalized( mal::Cross( vec_pos[1]-vec_pos[0], vec_pos[2]-vec_pos[0] ) );
    Vec3 ray_pos = mal::Rcp<Real>(3.0f)*(vec_pos[0]+vec_pos[1]+vec_pos[2]) - 0.002f*normal; //potentially-inside point if orientation is CCW \todo 0.001 ABSOLUTE EPSILON, ARBITRARY 1mm, SHOULD BE RELATIVE?!?!?!
    Vec3 ray_dir( normal );
    Real ray_length( 3*mal::Norm( m_AABB.GetHalfSizes() ) );
    // Count ray-tri hits
    uint32 num_hits(0);
    geo::np::RayHit3 rh;
    for( unsigned int it_tri=0; it_tri<m_addT.size(); it_tri++ )
        if( geo::np::RayCast_Triangle3_DoubleSided( ray_pos, ray_dir, Interval(0,ray_length),
                                                    m_addV[ m_addT[it_tri].m_vecVID[0] ].m_Pos,
                                                    m_addV[ m_addT[it_tri].m_vecVID[1] ].m_Pos,
                                                    m_addV[ m_addT[it_tri].m_vecVID[2] ].m_Pos,
                                                    rh ) )
            num_hits++;
    // If there's an even num_hits, ray_pos was OUTSIDE the closed surface, and thus normals were CW
    if( num_hits % 2 == 0 )
    {
        GEO_LOG_WARNING("EditableTriSurfaceShape3::RebuildClosed() enforcing CCW triangles");
        for( unsigned int it_tri=0; it_tri<m_addT.size(); it_tri++ )
        {
            //Swap VID 1,2 and NTID 0,2 (corresponding to e0 and e2)
            std::swap( m_addT[it_tri].m_vecVID[1], m_addT[it_tri].m_vecVID[2] );
            std::swap( m_addT[it_tri].m_vecNTID[0], m_addT[it_tri].m_vecNTID[2] );
        }
    }
}

// \todo CONSIDER keeping map Old2New for viz purposes (coincident vertices may have different texcoords/tangentspaces)
bool EditableTriSurfaceShape3::FixDegenerateVertices()
{
    //---- Merge coincident vertices
    unsigned int num_merged(0);
#define __ENABLE_MERGE
#ifdef __ENABLE_MERGE
    {
        // (1) For each vtx, store to which "parent" vertex does it merge (a la MF-set)
        std::vector< trisurface3_feature_index_type > vecParentVID( m_addV.size() );
        // Initially, each vertex is unique
        for( uint32 it_v=0; it_v < m_addV.size(); it_v++ ) vecParentVID[it_v] = it_v;
        // Use SH only for non-trivial meshes, otherwise the O(n) complexity is not worth the grid overhead
        //TEMP: Do it always, by now, shorter code... if( m_addV.size() > MIN_VERTICES_TO_USE_SH )
        {
            // Create SH (AABB assumed up-to-date from EndFaces()->Rebuild())
            uint32 dimensions[3] = { MAX_DIMENSION_SH, MAX_DIMENSION_SH, MAX_DIMENSION_SH };
            // Compute SH resolution for the largest extent
            Vec3 sizes( m_AABB.GetMax()-m_AABB.GetMin() );
            Real resolution = mal::Max(sizes) / MAX_DIMENSION_SH;
            for( int i=0; i<3; i++ )
                dimensions[i] = mal::Clamp<uint32>( mal::IntPart<uint32>( mal::Ceil( sizes[i] / resolution ) ), 1, MAX_DIMENSION_SH );
            //IMPORTANT, We force uint32 entry_index_type in GSimpleSpatialHash because even if trisurface3_feature_index_type = uint16, it must allocate 3*T edge-entries
            GSimpleSpatialHash<uint32, trisurface3_feature_index_type> ssh( m_addV.size(),
                                                                            dimensions,
                                                                            m_AABB.GetMin(), m_AABB.GetMax() );
            // GEO_LOG( "EditableTriSurfaceShape3::FixDegenerateVertices() SH dim = (%d,%d,%d)", dimensions[0], dimensions[1], dimensions[2] );

            // First pass: classify vertices
            ssh.BeginClassify();
            {
                for( uint32 it_v=0; it_v < m_addV.size(); it_v++ )
                    ssh.Classify( it_v, m_addV[it_v].m_Pos );
            }
            ssh.EndClassify();

            // Second pass, test overlaps
            /* For each vertex, we try to try to merge with all previous ones:
               - If none, the vertex remains unique.
               - Otherwise, we merge it into existing one.
               \todo No need for path-compression when merging into an already merged vtx, because we NEVER merge 2 already added vtx
            */
            uint32 num_tests_sh(0);
            ssh.BeginTestAdd();
            {
                std::vector<trisurface3_feature_index_type> vec_potential_overlaps;
                for( uint32 it_v1=0; it_v1 < m_addV.size(); it_v1++ )
                    if( ssh.TestAdd( it_v1,
                                     m_addV[it_v1].m_Pos,
                                     GEO_TRISS_EPSILON_LENGTH, //GEO_TRISS_EPSILON_LENGTH_SQ, //\todo CRITICAL: this SHOULD BE GEO_TRISS_EPSILON_LENGTH, not squared... MANY topology problems could derive from this shit!!
                                     vec_potential_overlaps ) )
                        for( uint32 it_pov=0; it_pov < vec_potential_overlaps.size(); it_pov++ )
                        {
                            uint32 it_v2 = vec_potential_overlaps[it_pov];
                            if( mal::NormSq(m_addV[it_v1].m_Pos-m_addV[it_v2].m_Pos) < GEO_TRISS_EPSILON_LENGTH_SQ ) //\todo Could use some complex merging condition here...
                                vecParentVID[it_v1] = vecParentVID[it_v2];
                            num_tests_sh++;
                        }
            }
            ssh.EndTestAdd();
            // GEO_LOG( "EditableTriSurfaceShape3::FixDegenerateVertices() SH tests = %d << %d", num_tests_sh, uint32(m_addV.size() * (m_addV.size()-1)) / 2 );
        }
        // (2) Compactify vertices using merge-to relationship
        for( uint32 it_v=0; it_v < m_addV.size(); it_v++ )
        {
            if( vecParentVID[it_v] == it_v )
            {
                // Unique vtx, compactify moving it num_merged slots to the left
                vecParentVID[it_v] = it_v - num_merged;
                m_addV[ vecParentVID[it_v] ] = m_addV[it_v];
            }
            else
            {
                // Non-unique vtx, Map to parent's ALREADY compactified index
                vecParentVID[it_v] = vecParentVID[ vecParentVID[it_v] ];
                num_merged++;
            }
        }
        m_addV.resize( m_addV.size() - num_merged ); //crop last num_merged V
        // (3) Remap T vid
        for( uint32 it_tri=0; it_tri < m_addT.size(); it_tri++ )
            for( int it_vit=0; it_vit < 3; it_vit++ )
                m_addT[it_tri].m_vecVID[it_vit] = vecParentVID[ m_addT[it_tri].m_vecVID[it_vit] ];
    }
#endif //__ENABLE_MERGE
    //---- Remove 0-degree vertices
    unsigned int num_isolated(0);
    {
        // Compute V degree (topology is NOT available, only m_addT data is up to date)
        std::vector< int > vecVertexDegree( m_addV.size(), 0 );
        for( unsigned int it_tri=0; it_tri<m_addT.size(); it_tri++ )
            for( int it_vit=0; it_vit<3; it_vit++ )
                vecVertexDegree[ m_addT[it_tri].m_vecVID[it_vit] ]++;
        // Remove V with 0-degree and remap to their new index
        std::vector< trisurface3_feature_index_type > vecOld2NewVID( m_addV.size() );
        for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
        {
            if( vecVertexDegree[it_v] > 0 )
            {
                vecOld2NewVID[it_v] = it_v - num_isolated;
                m_addV[ vecOld2NewVID[it_v] ] = m_addV[it_v];
            }
            else
            {
                vecOld2NewVID[it_v] = cTriSurface3_InvalidFeatureIndex;
                num_isolated++;
            }
        }
        m_addV.resize( m_addV.size() - num_isolated ); //crop last num_isolated V
        // Remap T vid
        for( unsigned int it_tri=0; it_tri<m_addT.size(); it_tri++ )
            for( int it_vit=0; it_vit<3; it_vit++ )
                m_addT[it_tri].m_vecVID[it_vit] = vecOld2NewVID[ m_addT[it_tri].m_vecVID[it_vit] ];
        // GEO_LOG_WARNING("%d 0-degree vertices removed",num_isolated);
    }
    if( num_merged > 0 || num_isolated > 0 )
        GEO_LOG( "EditableTriSurfaceShape3::FixDegenerateVertices() merged %u, isolated %u", num_merged, num_isolated );
    return num_merged > 0 || num_isolated > 0;
}

/*\note This breaks triangle topology/neighbours, must
        RebuildTopology() afterwards as done in EndEdition()
  \todo Removing degenerate tri with repeated VID may produce isolated
        0-degree vertices on the border, if the TriSurfaceShape3 is not
        closed

  \todo Naive sliver removal opens a hole in the mesh, so in addition,
  we collapse the shortest edge E to force vertex weld in the next
  iteration. We DO NOT modify the neighbour triangle across E, which
  becomes collapsed and will be hopefully removed in a later
  iteration.
*/
bool EditableTriSurfaceShape3::FixDegenerateTriangles()
{
    GEO_LOG("g_TTS_SliverThreshold = %f",g_TTS_SliverThreshold);

    // Remove triangles that have repeated VID, after FixDegenerateVertices()
    unsigned int num_collapsed(0);
    unsigned int num_slivers(0);
    {
        // std::vector< trisurface3_feature_index_type > vecOld2NewTID( m_addT.size() ); \todo nothing to remap, by now
        for( unsigned int it_tri=0; it_tri<m_addT.size(); it_tri++ )
        {
            uint32 vid0( m_addT[it_tri].m_vecVID[0] );
            uint32 vid1( m_addT[it_tri].m_vecVID[1] );
            uint32 vid2( m_addT[it_tri].m_vecVID[2] );
            bool bCollapsed( vid0 == vid1 || vid1 == vid2 || vid2 == vid0 );
            Vec3 d01( m_addV[vid1].m_Pos-m_addV[vid0].m_Pos );
            Vec3 d02( m_addV[vid2].m_Pos-m_addV[vid0].m_Pos );
            Vec3 d12( m_addV[vid2].m_Pos-m_addV[vid1].m_Pos );
            //\todo This sliver condition is WRONG and ignores d12!! it actually removed bad tri after Clip(), but complicates TryToCloseHoles() afterwards...
            // bool bSliver( mal::Abs( mal::Norm(mal::Cross(d01,d02)) ) < g_TTS_SliverThreshold*mal::Norm(d01)*mal::Norm(d02) );
            Real a( mal::Norm(d01) );
            Real b( mal::Norm(d02) );
            Real c( mal::Norm(d12) );
            Real circumradius( a*b*c / mal::Sqrt( (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c) ) ); //\see http://mathworld.wolfram.com/Circumradius.html
            bool bSliver( circumradius / mal::Min(a,mal::Min(b,c)) > g_TTS_SliverThreshold ); //larger ratio => worse sliver (min value = 0.5)
            if( !bCollapsed && !bSliver )
                m_addT[it_tri-num_collapsed-num_slivers] = m_addT[it_tri];
            else if( bCollapsed )
                num_collapsed++;
            else if( bSliver )
            {
                num_slivers++;
                // Find shortest edge
                int vid_a(0),vid_b(0);
                Real length01_sq(mal::NormSq(d01));
                Real length02_sq(mal::NormSq(d02));
                Real length12_sq(mal::NormSq(d12));
                if( length01_sq <= length02_sq && length01_sq <= length12_sq ) { vid_a = vid0; vid_b = vid1; }
                else if( length02_sq <= length01_sq && length02_sq <= length12_sq ) { vid_a = vid0; vid_b = vid2; }
                else { vid_a = vid1; vid_b = vid2; }
                // Snap shortest edge
                Vec3& a( m_addV[ vid_a ].m_Pos );
                Vec3& b( m_addV[ vid_b ].m_Pos );
                Vec3 snapped_ab = Real(0.5)*(a+b);
                a = b = snapped_ab;
            }
        }
        m_addT.resize( m_addT.size() - num_collapsed - num_slivers ); //crop
    }

    if( num_collapsed > 0 || num_slivers > 0 )
        GEO_LOG( "EditableTriSurfaceShape3::FixDegenerateTriangles() collapsed %u, slivers %u", num_collapsed, num_slivers );
    return num_collapsed > 0 || num_slivers > 0;
}

/*\note This DOES NOT AFFECT triangle topology/neighbours, can be
  safely done after RebuildTopology() but BEFORE computing if closed.
*/
bool EditableTriSurfaceShape3::RemoveIsolatedTriangles()
{
    // Remove triangles that have all edges open
    unsigned int num_isolated(0);
    {
        for( unsigned int it_tri=0; it_tri<m_addT.size(); it_tri++ )
        {
            unsigned int num_open_edges(0);
            for( unsigned int it_neighbour=0; it_neighbour<3; it_neighbour++ )
                if( m_addT[it_tri].m_vecNTID[it_neighbour] == cTriSurface3_InvalidFeatureIndex )
                    num_open_edges++;
            if( num_open_edges < 2 )
                m_addT[it_tri-num_isolated] = m_addT[it_tri];
            else
                num_isolated++;
        }
        m_addT.resize( m_addT.size() - num_isolated ); //crop last num_removed V
    }
    if( num_isolated > 0 )
        GEO_LOG( "EditableTriSurfaceShape3::RemoveIsolatedTriangles() isolated %u", num_isolated );
    return num_isolated > 0;
}


//-----------------------------------------------------------------------------
// Free functions
//-----------------------------------------------------------------------------
Real ComputeVolume( const TriSurfaceShape3& tss3, const Transform3& transform, const Vec3* vec_sdof )
{
    const Vec3* actual_sdof( ( 0 != vec_sdof ) ? vec_sdof : tss3.GetVecDefaultSDOF() );
    Real volume_6(0);
    for( unsigned int it_tri=0; it_tri<tss3.GetNumT(); it_tri++ )
        volume_6 += mal::Det(mal::GMat3x3_From_Columns( tss3.V_Pos( tss3.T_VID(it_tri,0), actual_sdof ),
                                                        tss3.V_Pos( tss3.T_VID(it_tri,1), actual_sdof ),
                                                        tss3.V_Pos( tss3.T_VID(it_tri,2), actual_sdof ) ) );
    return volume_6/6;
}

//-----------------------------------------------------------------------------
//---- Make_TriSurfaceShape3_XXX Implementation
//-----------------------------------------------------------------------------

void Make_TriSurfaceShape3_Box( EditableTriSurfaceShape3 &etss3, const Vec3 &half_sizes )
{
    etss3.Clear();
    etss3.BeginEdition();
    {
        trisurface3_feature_index_type vid0 = etss3.AddVertex( Vec3( -half_sizes.x(), -half_sizes.y(), -half_sizes.z() ) );
        trisurface3_feature_index_type vid1 = etss3.AddVertex( Vec3( -half_sizes.x(),  half_sizes.y(), -half_sizes.z() ) );
        trisurface3_feature_index_type vid2 = etss3.AddVertex( Vec3(  half_sizes.x(),  half_sizes.y(), -half_sizes.z() ) );
        trisurface3_feature_index_type vid3 = etss3.AddVertex( Vec3(  half_sizes.x(), -half_sizes.y(), -half_sizes.z() ) );
        trisurface3_feature_index_type vid4 = etss3.AddVertex( Vec3( -half_sizes.x(), -half_sizes.y(),  half_sizes.z() ) );
        trisurface3_feature_index_type vid5 = etss3.AddVertex( Vec3( -half_sizes.x(),  half_sizes.y(),  half_sizes.z() ) );
        trisurface3_feature_index_type vid6 = etss3.AddVertex( Vec3(  half_sizes.x(),  half_sizes.y(),  half_sizes.z() ) );
        trisurface3_feature_index_type vid7 = etss3.AddVertex( Vec3(  half_sizes.x(), -half_sizes.y(),  half_sizes.z() ) );
        etss3.AddTriangle( vid0, vid1, vid3 );
        etss3.AddTriangle( vid0, vid4, vid1 );
        etss3.AddTriangle( vid0, vid3, vid4 );
        etss3.AddTriangle( vid1, vid2, vid3 );
        etss3.AddTriangle( vid1, vid4, vid5 );
        etss3.AddTriangle( vid1, vid5, vid6 );
        etss3.AddTriangle( vid1, vid6, vid2 );
        etss3.AddTriangle( vid2, vid6, vid3 );
        etss3.AddTriangle( vid3, vid6, vid7 );
        etss3.AddTriangle( vid3, vid7, vid4 );
        etss3.AddTriangle( vid4, vid6, vid5 );
        etss3.AddTriangle( vid4, vid7, vid6 );
    }
    etss3.EndEdition();
}

} //namespace geo
