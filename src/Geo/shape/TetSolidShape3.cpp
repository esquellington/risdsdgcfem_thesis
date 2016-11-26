#include "TetSolidShape3.h"
#include <memory.h> //TEMPORAL req by memcpy()
#include "../util/GSimpleSpatialHash.h"

// DCR \todo move elsewhere
#include <Geo/shape/TriSurfaceShape3.h>
#include <Geo/np/Overlap.h>

#ifdef __GEO_TETSS_ENABLE_BVH
#  include <boost/bind.hpp> //for ILSS
#endif

#define MIN_EDGES_TO_USE_SH 100
#define MAX_DIMENSION_SH 32
#define GEO_TETSS_EPSILON_LENGTH 0.001f //\todo CONSIDER using relative epsilon, an optional parameter OR even the Context epsilon...
#define GEO_TETSS_EPSILON_LENGTH_SQ mal::Sq(GEO_TETSS_EPSILON_LENGTH)

namespace geo
{

//-----------------------------------------------------------------------------
//---- TetSolidShape3 Implementation
//-----------------------------------------------------------------------------

TetSolidShape3::TetSolidShape3()
: m_NumV(0), m_NumT(0), m_NumBF(0), m_NumL(0)
, m_vecPoints(0), m_vecV(0), m_vecT(0), m_vecBF(0), m_vecL(0)
, m_pBuffer(0)
, m_pDCR(0)
#ifdef __GEO_TETSS_ENABLE_BVH
, m_pBVH(0)
#endif
{}

TetSolidShape3::~TetSolidShape3()
{
    ClearBakedData();
}

void TetSolidShape3::ClearBakedData()
{
    if( m_pBuffer ) delete [] m_pBuffer;
    m_pBuffer = 0;
    m_NumV = 0; m_vecPoints = 0; m_vecV = 0;
    m_NumT = 0; m_vecT = 0;
    m_NumBF = 0; m_vecBF = 0;
    m_NumL = 0; m_vecL = 0;
    // DCR is strictly nonshared, by now
    if( m_pDCR ) delete m_pDCR;
    m_pDCR = 0;
#ifdef __GEO_TETSS_ENABLE_BVH
    // BVH is strictly nonshared, by now
    if( m_pBVH ) delete m_pBVH;
    m_pBVH = 0;
#endif
}

void TetSolidShape3::ComputeBVD( bv::BoundingVolume3 &bv, const transform_type &transform, const sdof_type *vec_sdof ) const
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
void TetSolidShape3::SetBakedData( bool b_shared,
                                   uint32 num_v, uint32 num_t, uint32 num_bf, uint32 num_l,
                                   const Vec3 *vec_points, const vertex_type *vec_v, const tetrahedron_type *vec_t, const boundary_face_type *vec_bf, const tetrahedron_layer_type *vec_l )
{
    GEO_ASSERT( m_NumV == 0 && m_vecPoints == 0 && m_vecV == 0 && num_v > 0
                && m_NumT == 0 && m_vecT == 0 && num_t > 0
                && m_NumBF == 0 && m_vecBF == 0 && num_bf > 0
                && m_NumL == 0 && num_l > 0
                && m_pBuffer == 0 );
    m_NumV = num_v;
    m_NumT = num_t;
    m_NumBF = num_bf;
    m_NumL = num_l;
    if( b_shared )
    {
        m_vecPoints = vec_points;
        m_vecV = vec_v;
        m_vecT = vec_t;
        m_vecBF = vec_bf;
        m_vecL = vec_l;
    }
    else
    {
        // Alloc and fill single buffer
        size_t size_points_4aligned = 4*((sizeof(Vec3)*m_NumV+3)/4);
        size_t size_v_4aligned = 4*((sizeof(vertex_type)*m_NumV+3)/4);
        size_t size_t_4aligned = 4*((sizeof(tetrahedron_type)*m_NumT+3)/4);
        size_t size_bf_4aligned = 4*((sizeof(boundary_face_type)*m_NumBF+3)/4);
        size_t size_l_4aligned = 4*((sizeof(tetrahedron_layer_type)*m_NumL+3)/4);
        size_t total_size_4aligned = size_points_4aligned + size_v_4aligned + size_t_4aligned + size_bf_4aligned + size_l_4aligned;
        m_pBuffer = new uint32[ total_size_4aligned ];
        Vec3 *p_buffer_points = reinterpret_cast<Vec3*>( &m_pBuffer[0] );
        vertex_type *p_buffer_v = reinterpret_cast<vertex_type*>( &m_pBuffer[ size_points_4aligned ] );
        tetrahedron_type *p_buffer_t = reinterpret_cast<tetrahedron_type*>( &m_pBuffer[ size_points_4aligned + size_v_4aligned ] );
        boundary_face_type *p_buffer_bf = reinterpret_cast<boundary_face_type*>( &m_pBuffer[ size_points_4aligned + size_v_4aligned + size_t_4aligned ] );
        tetrahedron_layer_type *p_buffer_l = reinterpret_cast<tetrahedron_layer_type*>( &m_pBuffer[ size_points_4aligned + size_v_4aligned + size_t_4aligned + size_bf_4aligned ] );
        memcpy( p_buffer_points, vec_points, sizeof(Vec3)*m_NumV );
        memcpy( p_buffer_v, vec_v, sizeof(vertex_type)*m_NumV );
        memcpy( p_buffer_t, vec_t, sizeof(tetrahedron_type)*m_NumT );
        memcpy( p_buffer_bf, vec_bf, sizeof(boundary_face_type)*m_NumBF );
        memcpy( p_buffer_l, vec_l, sizeof(tetrahedron_layer_type)*m_NumL );
        // Save const pointers
        m_vecPoints = p_buffer_points;
        m_vecV = p_buffer_v;
        m_vecT = p_buffer_t;
        m_vecBF = p_buffer_bf;
        m_vecL = p_buffer_l;
    }

    /*TEMPORAL
    GEO_LOG_WARNING( "SetBakedData( %d %d %d )", num_vertices, num_tetrahedrons, num_boundary_faces );
    for( unsigned int it_tet=0; it_tet < GetNumT(); it_tet++ )
    {

        GEO_LOG_WARNING( "Tet[%d] = (%d,%d,%d,%d)",
                         it_tet,
                         m_vecT[it_tet].m_vecVID[0],
                         m_vecT[it_tet].m_vecVID[1],
                         m_vecT[it_tet].m_vecVID[2],
                         m_vecT[it_tet].m_vecVID[3] );
    }
    */
}

Vec3 TetSolidShape3::T_Barycenter_0( uint32 tid ) const
{
    return Real(0.25) * ( V_Pos_0( m_vecT[tid].m_vecVID[0] ) +
                          V_Pos_0( m_vecT[tid].m_vecVID[1] ) +
                          V_Pos_0( m_vecT[tid].m_vecVID[2] ) +
                          V_Pos_0( m_vecT[tid].m_vecVID[3] ) );
}

Vec3 TetSolidShape3::T_Barycenter( uint32 tid, const sdof_type *vec_sdof ) const
{
    return Real(0.25) * ( vec_sdof[ m_vecT[tid].m_vecVID[0] ] +
                          vec_sdof[ m_vecT[tid].m_vecVID[1] ] +
                          vec_sdof[ m_vecT[tid].m_vecVID[2] ] +
                          vec_sdof[ m_vecT[tid].m_vecVID[3] ] );
}

/* Enumerate 6 edges 0..1 in a given tet:
   -3 edges in "frame" (0->1,0->2,0->3)
   -3 edges in triangle (1->2->3->1)
*/
void TetSolidShape3::T_Edge_0( uint32 tid, uint32 eit, Vec3& p0, Vec3& p1 ) const
{
    GEO_LOG_ASSERT( tid < m_NumT && eit < 6, "TetSolidShape3(%d,%d) out of range!", tid, eit );
    switch(eit)
    {
        // Frame 0->1,2,3
    case 0: p0 = V_Pos_0( T_VID(tid,0) ); p1 = V_Pos_0( T_VID(tid,1) ); break;
    case 1: p0 = V_Pos_0( T_VID(tid,0) ); p1 = V_Pos_0( T_VID(tid,2) ); break;
    case 2: p0 = V_Pos_0( T_VID(tid,0) ); p1 = V_Pos_0( T_VID(tid,3) ); break;
        // Triangle 1->2->3->1
    case 3: p0 = V_Pos_0( T_VID(tid,1) ); p1 = V_Pos_0( T_VID(tid,2) ); break;
    case 4: p0 = V_Pos_0( T_VID(tid,2) ); p1 = V_Pos_0( T_VID(tid,3) ); break;
    case 5: p0 = V_Pos_0( T_VID(tid,3) ); p1 = V_Pos_0( T_VID(tid,1) ); break;
    default: break;
    }
}

Vec3 TetSolidShape3::BF_Barycenter( uint32 bfid, const sdof_type *vec_sdof ) const
{
    return Real(1.0/3.0) * ( vec_sdof[ m_vecBF[bfid].m_vecVID[0] ] +
                             vec_sdof[ m_vecBF[bfid].m_vecVID[1] ] +
                             vec_sdof[ m_vecBF[bfid].m_vecVID[2] ] );
}

//-----------------------------------------------------------------------------
//---- EditableTetSolidShape3 Implementation
//-----------------------------------------------------------------------------
EditableTetSolidShape3::EditableTetSolidShape3()
: m_IsBeingEdited(false)
{
}

EditableTetSolidShape3::~EditableTetSolidShape3()
{
    Clear();
}

void EditableTetSolidShape3::Clear()
{
    GEO_ASSERT( !IsBeingEdited() );
    ClearEditData();
    TetSolidShape3::ClearBakedData();
}

void EditableTetSolidShape3::Set( const TetSolidShape3 &tss3 )
{
    GEO_ASSERT( IsBeingEdited() );
    Clear();
    // Copy other mesh baked data
    for( unsigned int i=0; i<tss3.GetNumV(); i++ )
        m_addV.push_back( editable_vertex_type( tss3.GetVecPoints()[i] ) );
    for( unsigned int i=0; i<tss3.GetNumT(); i++ )
        m_addT.push_back( editable_tetrahedron_type( tss3.GetVecT()[i].m_vecVID[0],
                                                     tss3.GetVecT()[i].m_vecVID[1],
                                                     tss3.GetVecT()[i].m_vecVID[2],
                                                     tss3.GetVecT()[i].m_vecVID[3] ) );
    for( unsigned int i=0; i<tss3.GetNumBF(); i++ )
        m_addBF.push_back( editable_boundary_face_type( tss3.GetVecBF()[i].m_vecVID[0],
                                                        tss3.GetVecBF()[i].m_vecVID[1],
                                                        tss3.GetVecBF()[i].m_vecVID[2] ) );
    EndEdition();
}

void EditableTetSolidShape3::BeginEdition()
{
    GEO_ASSERT( !IsBeingEdited() );
    // Clear anything added outside Begin/End
    ClearEditData();
    // Copy baked data
    for( unsigned int i=0; i<m_NumV; i++ )
        m_addV.push_back( editable_vertex_type( m_vecPoints[i] ) );
    for( unsigned int i=0; i<m_NumT; i++ )
        m_addT.push_back( editable_tetrahedron_type( m_vecT[i].m_vecVID[0],
                                                     m_vecT[i].m_vecVID[1],
                                                     m_vecT[i].m_vecVID[2],
                                                     m_vecT[i].m_vecVID[3] ) );
    for( unsigned int i=0; i<m_NumBF; i++ )
        m_addBF.push_back( editable_boundary_face_type( m_vecBF[i].m_vecVID[0],
                                                        m_vecBF[i].m_vecVID[1],
                                                        m_vecBF[i].m_vecVID[2] ) );
    // Delete baked data
    TetSolidShape3::ClearBakedData();
    m_IsBeingEdited = true;
}

bool EditableTetSolidShape3::EndEdition()
{
    GEO_ASSERT( IsBeingEdited() );
    bool bOk = RebuildAllFromEditData();
    m_IsBeingEdited = false;
    return bOk;
}

tetsolid3_feature_index_type EditableTetSolidShape3::AddVertex( const Vec3 &point )
{
    GEO_ASSERT( IsBeingEdited() );
    GEO_ASSERT( m_addV.size() < cTetSolid3_InvalidFeatureIndex-1 );
    m_addV.push_back( editable_vertex_type( point ) );
    return tetsolid3_feature_index_type(m_addV.size()-1);
}

tetsolid3_feature_index_type EditableTetSolidShape3::AddTetrahedron( uint32 vid0, uint32 vid1, uint32 vid2, uint32 vid3 )
{
    GEO_ASSERT( IsBeingEdited() );
    //GEO_LOG_WARNING("EditableTetSolidShape3::AddTetrahedron(%d,%d,%d,%d) start...", vid0, vid1, vid2, vid3 );
    //\todo Consider looking for existing T
    GEO_ASSERT( m_addT.size() < cTetSolid3_InvalidFeatureIndex-1 );
    GEO_ASSERT( vid0 < m_addV.size() && vid1 < m_addV.size() && vid2 < m_addV.size() && vid3 < m_addV.size() );
    Real det( mal::Det( mal::GMat3x3_From_Columns( m_addV[vid1].m_Pos - m_addV[vid0].m_Pos,
                                                   m_addV[vid2].m_Pos - m_addV[vid0].m_Pos,
                                                   m_addV[vid3].m_Pos - m_addV[vid0].m_Pos ) ) );
    if( det <= 0 )
    {
        Vec3f v0( m_addV[vid0].m_Pos );
        Vec3f v1( m_addV[vid1].m_Pos );
        Vec3f v2( m_addV[vid2].m_Pos );
        Vec3f v3( m_addV[vid3].m_Pos );
        GEO_LOG_WARNING("EditableTetSolidShape3::AddTetrahedron(%d,%d,%d,%d) ignored!\n\tdet = %f <= 0\n\tv0 = (%f,%f,%f)\n\tv1 = (%f,%f,%f)\n\tv2 = (%f,%f,%f)\n\tv3 = (%f,%f,%f)",
                        vid0, vid1, vid2, vid3,
                        det,
                        v0[0], v0[1], v0[2],
                        v1[0], v1[1], v1[2],
                        v2[0], v2[1], v2[2],
                        v3[0], v3[1], v3[2] );
        return cTetSolid3_InvalidFeatureIndex;
    }

    //GEO_ASSERT( det > 0 );
    // Alloc new TID
    tetsolid3_feature_index_type tid( m_addT.size() );
    // Create Tetrahedron
    m_addT.push_back( editable_tetrahedron_type( vid0, vid1, vid2, vid3 ) );
    return tid;
}

/* DEPRECATED
tetsolid3_feature_index_type EditableTetSolidShape3::AddBoundaryFace( tetsolid3_feature_index_type vid0, tetsolid3_feature_index_type vid1, tetsolid3_feature_index_type vid2 )
{
    GEO_ASSERT( IsBeingEdited() );
    //\todo Consider looking for existing T
    GEO_ASSERT( m_addBF.size() < cTetSolid3_InvalidFeatureIndex-1 );
    GEO_ASSERT( vid0 < m_addV.size() && vid1 < m_addV.size() && vid2 < m_addV.size() );
    // Alloc new TID
    tetsolid3_feature_index_type bfid( m_addBF.size() );
    m_addBF.push_back( editable_boundary_face_type( vid0, vid1, vid2 ) );
    return bfid;
}
*/

/* Move an existing vertex.
   \note This will NOT change topology... unless vertices are made
   coincident! therefore we MUST fully-rebuild in EndEdition()
   afterwards
*/
void EditableTetSolidShape3::SetVertex( uint32 vid, const Vec3& point )
{
    GEO_ASSERT( IsBeingEdited() );
    GEO_ASSERT( vid < m_addV.size() );
    GEO_ASSERT( !IsNaN(point) );
    m_addV[vid].m_Pos = point;
}

void EditableTetSolidShape3::Transform( const Transform3& tr )
{
    GEO_ASSERT( IsBeingEdited() );
    GEO_ASSERT( !IsNaN(tr) );
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
        m_addV[it_v].m_Pos = tr * m_addV[it_v].m_Pos;
}

/* Slow but easy:
  - Rebuild the whole EditableTetSolidShape3 adding all existing
    tetrahedrons except the ones in vec_tid
  - Note that spurious vertices only adjacent to removed tetrahedrons may
    be added in the first pass, but will be removed in EndEdition()
    FixDegenerateVertices
  - \todo Optimizations:
    - Current cost is is O(#T*num_tid)
    - if vec_tid were sorted, the cost would become O(#T)
      - this is HIGLY PROBABLE if vec_tid comes from a "filtering
        predicate" run on all tet 0..#T-1 in order
*/
bool EditableTetSolidShape3::RemoveTetrahedrons( uint32 num_tid, const uint32* vec_tid )
{
    GEO_ASSERT( !IsBeingEdited() );
    GEO_ASSERT( m_NumT > 0 && num_tid < m_NumT && vec_tid != 0 );
    ClearEditData(); //we'll use baked data to re-add all edit data
    // Re-add all Vertices
    for( unsigned int it_v=0; it_v<m_NumV; it_v++ )
        m_addV.push_back( editable_vertex_type( m_vecPoints[it_v] ) );
    // Re-add non-removed tetrahedrons
    for( unsigned int it_tet=0; it_tet<m_NumT; it_tet++ )
    {
        bool bKeep(true);
        for( unsigned int it_remove_tid=0; bKeep && it_remove_tid < num_tid; it_remove_tid++ )
            bKeep = it_tet != vec_tid[it_remove_tid];
        if( bKeep )
            m_addT.push_back( editable_tetrahedron_type( m_vecT[it_tet].m_vecVID[0],
                                                         m_vecT[it_tet].m_vecVID[1],
                                                         m_vecT[it_tet].m_vecVID[2],
                                                         m_vecT[it_tet].m_vecVID[3] ) );
    }
    return RebuildAllFromEditData();
}

/* From EditData rebuild topology, create splitV with
   UNIQUE vertices per tet, and remap/merge them according to
   incident tet adjacency (merge splitV if edge-adjacent tets
   refer them). After all tets have been considered, the splitV
   will have been merged into disjoint sets, split them
   geometrically, reindex tet VID and rebuild all (to rebuild tet
   topology and BF). If required, add originalVID parent to each
   splitV or alternatively store list of splitV per originalV

   \note We use the pre-split topology during SV merge, but discard it
   and rebuild it completely afterwards to account for SV.

   \todo If we require face-adjacency instead of edge-adjacency to
   merge splitV, we will AUTOMATICALLY handle non-manifold edge
   vertices INDEPENDENTLY (so that 1 can be split and the other
   merged, which CAN HAPPEN)

   \todo Hopefully, this will ALSO help splitting unwanted
   face-adjacency (receive a list of faces to split explicitly,
   computed externally (eg, in
   Fit_TetSolidShape3_To_TriSurfaceShape3) and avoid using them to
   merge splitV)

   // Original idea
   //---- Split non-manifold V
   // Gather per-V 1-ring adjacency
   // Compute per-V 1-ring disjoint components
   // Split V into its disjoint components
   //---- Split non-manifold E
   // Gather per-E 1-ring adjacency (\note strictly E-adjacency, not V-adjacency
*/
bool EditableTetSolidShape3::FixNonManifoldFeatures( const std::vector<tetsolid3_face_id_type>& vec_split_fid )
{
    GEO_ASSERT( !IsBeingEdited() );
    BeginEdition();
    RebuildTopology();
    // Init all vertices as split
    std::vector< editable_vertex_type > vecSV;
    vecSV.reserve( 4*m_addT.size() );
    for( auto tet : m_addT )
        for( unsigned int it_vit=0; it_vit<4; it_vit++ )
            vecSV.push_back( m_addV[ tet.m_vecVID[it_vit] ] );
    // For each SV, store which "parent" svid does it merge to (a la MF-set)
    std::vector<tetsolid3_feature_index_type> vecParentSVID( vecSV.size() );
    for( unsigned int it_sv=0; it_sv < vecSV.size(); it_sv++ ) vecParentSVID[it_sv] = it_sv;
    //\todo multi-pass seems to be required to properly merge all SV
    //in and bunny and armadillo... otherwise some SV remain unmerged,
    //which signals a problem in the single-pass algorithm
    for( int i=0; i<5; i++ )
    {
        // For each tet, merge its svid with face-adjacent neighbour svid
        for( unsigned int it_tet=0; it_tet < m_addT.size(); it_tet++ )
        {
            // GEO_LOG("* Merging T %u", it_tet );
            const editable_tetrahedron_type& tet( m_addT[it_tet] );
            //merge all splitV shared through faces to the smallest splitVID
            for( unsigned int it_ntit=0; it_ntit<4; it_ntit++ )
            {
                uint32 ntid( tet.m_vecNTID[it_ntit] );
                // GEO_LOG("** Merging NTID %u", ntid );
                if( ntid != cTetSolid3_InvalidFeatureIndex
//                && it_tet < ntid //Avoid reprocessing the same face \todo Two-sided faces SEEM to be necessary for proper merge propagation... investigate
                    && std::find( vec_split_fid.begin(), vec_split_fid.end(), tetsolid3_face_id_type(it_tet,it_ntit) ) == vec_split_fid.end() //ignore split-F
                    ) // \todo && UniqueFace(it_tet,it_ntit) is NOT in list of split faces no...
                {
                    for( unsigned int it_vif=0; it_vif<3; it_vif++ )
                    {
                        // GEO_LOG("*** Merging VIF %u", it_vif );
                        uint32 vit( (it_ntit+it_vif+1)%4 ); //(+1,+2,+3) %4 are the vif of the opposite face to the vertex it_ntit
                        uint32 vid( tet.m_vecVID[vit] );
                        uint32 svid( 4*it_tet + vit ); //split vid
                        uint32 nsvid( cTetSolid3_InvalidFeatureIndex ); //matching neightbour split vid
                        for( unsigned int it_nvit=0; it_nvit<4; it_nvit++ )
                            if( m_addT[ntid].m_vecVID[it_nvit] == vid )
                                nsvid = 4*ntid + it_nvit;
                        GEO_ASSERT( nsvid != cTetSolid3_InvalidFeatureIndex );
                        /* merge-to-smallest guarantees sequential remapping will work later
                           if( svid < nsvid )
                           {
                           vecParentSVID[nsvid] = svid;
                           GEO_LOG("**** Merging V %u < %u", svid, nsvid );
                           }
                           else
                           {
                           vecParentSVID[svid] = nsvid;
                           GEO_LOG("**** Merging V %u < %u", nsvid, svid );
                           }
                        */
                        //IMPORTANT: merge-to-smallest must be done with
                        //current parent SVID, not with actual SVID, as
                        //any of the involved SVID may have already been
                        //merged and ignoring this fact causes missing
                        //merge chains (eg, SVID 13 being merged into SVID
                        //1 first and SVID 5 afterwards in Box 1x1x1 test
                        if( vecParentSVID[svid] < vecParentSVID[nsvid] ) vecParentSVID[nsvid] = vecParentSVID[svid];
                        else vecParentSVID[svid] = vecParentSVID[nsvid];
                    }
                }
            }
        }
    }
    // (2) Compactify svid using merge-to relationship
    uint32 num_merged(0);
    for( uint32 it_sv=0; it_sv < vecSV.size(); it_sv++ )
    {
        if( vecParentSVID[it_sv] == it_sv )
        {
            // Unique vtx, compactify moving it num_merged slots to the left
            vecParentSVID[it_sv] = it_sv - num_merged;
            vecSV[ vecParentSVID[it_sv] ] = vecSV[it_sv];
        }
        else
        {
            // Non-unique vtx, Map to parent's ALREADY compactified index
            vecParentSVID[it_sv] = vecParentSVID[ vecParentSVID[it_sv] ];
            num_merged++;
        }
    }
    vecSV.resize( vecSV.size() - num_merged ); //crop last num_merged SV
    // Swap old addV with vecSV
    std::swap(m_addV,vecSV);
    // (3) Remap T svid
    for( uint32 it_tet=0; it_tet < m_addT.size(); it_tet++ )
        for( int it_vit=0; it_vit < 4; it_vit++ )
            m_addT[it_tet].m_vecVID[it_vit] = vecParentSVID[ 4*it_tet + it_vit ];

    //\todo IMPORTANT: DISPLACE EACH SV TOWARDS ITS ADJACENT connected
    //component to avoid FixDegenerateVertices() to merge them when
    //vertex-weld is implemented!! => Can be done by "pulling" each SV
    //slightly towards its opposite faces, in a single per-tet global
    //iter (consider volume weighting displacement dir and enforcing
    //given displacement length = 2*epsilon_length.

    return EndEdition();
}

bool EditableTetSolidShape3::AddDCR( const IShape3* p_shape, const Transform3& tr_s2tss )
{
    if( p_shape->GetType() == eShape_TriSurface3 )
    {
        if( m_pDCR ) delete m_pDCR;
        m_pDCR = 0;
        const TriSurfaceShape3& surface( *static_cast<const TriSurfaceShape3*>(p_shape) );
        SetBakedDCR_StrictlyNonshared_UglyHack( Create_DCR_TetSolidShape3_From_TriSurfaceShape3( *this, Transform3::Identity(), GetVecDefaultSDOF(),
                                                                                                 surface, tr_s2tss, surface.GetVecDefaultSDOF() ) );
        return true;
    }
    else
        return false;
}

#ifdef __GEO_TETSS_ENABLE_BVH
bool EditableTetSolidShape3::AddBVH()
{
    if( m_pBVH ) delete m_pBVH;
    m_pBVH = new BVH_TetSolidShape3;
    uint32 num_entries( m_NumT );
    /*\todo NO! we use the BVH for full-volume queries too (in escalunya), so we CANNOT just consider DCR-covered elements...
      => Consider choosing behaviour with an enum-param... because for surface-only CD adding only DCR-covered elements should be more efficient...
      uint32 num_entries( m_pDCR ? m_pDCR->m_NumElements : m_NumP ); //TEMP: If there's a DCR, only use those P in the BVH
    */
    m_pBVH->Rebuild_TopDown( num_entries,
                             boost::bind<void>( &GEBV_TetSolidShape3_E<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
                                                this, transform_type::Identity(), GetVecDefaultSDOF(),
                                                _1, _2) );
    /*
    m_pBVH->Rebuild( num_entries,
                     boost::bind<void>( &GEBV_TetSolidShape3_E<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
                                        this, transform_type::Identity(), GetVecDefaultSDOF(),
                                        _1, _2) );
    */
    return true;
}
#endif

//---- Internal methods
bool EditableTetSolidShape3::RebuildAllFromEditData()
{
    FixDegeneracies();
    RebuildTopology();
    RebuildLayers();
    RebuildBakedData();
    ClearEditData();
    return true;
}

void EditableTetSolidShape3::ClearEditData()
{
    m_addV.clear();
    m_addT.clear();
    m_addBF.clear();
}

void EditableTetSolidShape3::RebuildBakedData()
{
    TetSolidShape3::ClearBakedData();
    // Bake Vertices
    Vec3* vecPoints( new Vec3[m_addV.size()] );
    vertex_type* vecV( new vertex_type[m_addV.size()] );
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
    {
        vecPoints[it_v] = m_addV[it_v].m_Pos;
        vecV[it_v] = m_addV[it_v];
    }
    // Bake Tetrahedrons
    tetrahedron_type* vecT( new tetrahedron_type[m_addT.size()] );
    for( unsigned int it_tet=0; it_tet<m_addT.size(); it_tet++ )
        vecT[it_tet] = m_addT[it_tet];
    // Bake boundary Faces
    boundary_face_type* vecBF( new boundary_face_type[m_addBF.size()] );
    for( unsigned int it_bf=0; it_bf<m_addBF.size(); it_bf++ )
        vecBF[it_bf] = m_addBF[it_bf];
    // Bake layers
    tetrahedron_layer_type* vecL( new tetrahedron_layer_type[m_addL.size()] );
    for( unsigned int it_l=0; it_l<m_addL.size(); it_l++ )
        vecL[it_l] = m_addL[it_l];

    // Set baked data as non-shared
    SetBakedData( false,
                  m_addV.size(), m_addT.size(), m_addBF.size(), m_addL.size(),
                  vecPoints, vecV, vecT, vecBF, vecL );

    // Clear temporal stuff
    delete [] vecL;
    delete [] vecBF;
    delete [] vecT;
    delete [] vecV;
    delete [] vecPoints;
}

bool EditableTetSolidShape3::FixDegeneracies()
{
    bool bFixDV = FixDegenerateVertices();
    bool bFixDT = FixDegenerateTetrahedrons();
    return bFixDV || bFixDT;
};

void EditableTetSolidShape3::RebuildTopology()
{
    if( m_addV.size() == 0 || m_addT.size() == 0 ) return;

    //---- Rebuild tetraedron adjacency
    // Compute AABB
    bv::AABB3 aabb( m_addV[0].m_Pos );
    for( unsigned int it_v=1; it_v<m_addV.size(); it_v++ ) aabb.Merge( m_addV[it_v].m_Pos );

    uint32 num_faces = 4*m_addT.size();
    // Use SH only for non-trivial meshes, otherwise the O(n) complexity is not worth the grid overhead
    //TEMP: Do it always, by now... if( num_faces > MIN_EDGES_TO_USE_SH )
    {
        // Create SH (AABB assumed up-to-date from EndFaces()->Rebuild())
        uint32 dimensions[3] = { MAX_DIMENSION_SH, MAX_DIMENSION_SH, MAX_DIMENSION_SH };
        // Compute SH resolution for the largest extent
        Vec3 sizes( aabb.GetMax()-aabb.GetMin() );
        Real resolution = mal::Max(sizes) / MAX_DIMENSION_SH;
        for( int i=0; i<3; i++ )
            dimensions[i] = mal::Clamp<uint32>( mal::IntPart<uint32>( mal::Ceil( sizes[i] / resolution ) ), 1, MAX_DIMENSION_SH );
        //IMPORTANT, We force uint32 entry_index_type in GSimpleSpatialHash because even if tetsolid3_feature_index_type = uint16, it must allocate 4*T face-entries
        GSimpleSpatialHash<uint32,tetsolid3_face_id_type> ssh( num_faces,
                                                               dimensions,
                                                               aabb.GetMin(), aabb.GetMax() );
        // GEO_LOG( "EditableTriSurfaceShape3::RebuildTopology() SH dim = (%d,%d,%d)", dimensions[0], dimensions[1], dimensions[2] );

        // First pass: classify tet-face barycenters
        ssh.BeginClassify();
        {
            for( unsigned int it_tet=0; it_tet<m_addT.size(); it_tet++ )
            {
                tetrahedron_type &tet( m_addT[it_tet] );
                for( unsigned int it_neighbour=0; it_neighbour<4; it_neighbour++ )
                    ssh.Classify( std::make_pair( it_tet, it_neighbour ),
                                  mal::Rcp<Real>(3.0f)*( m_addV[ tet.m_vecVID[(it_neighbour+1)%4] ].m_Pos
                                                         + m_addV[ tet.m_vecVID[(it_neighbour+2)%4] ].m_Pos
                                                         + m_addV[ tet.m_vecVID[(it_neighbour+3)%4] ].m_Pos ) );
                //\todo SHIIT... it's not so easy... see how to cycle over a tet, use it also in DAPD! ... consider map 4x3 so that given a vit = [0..3] we have the OPPOSITE vit in CCW face order!!
            }
        }
        ssh.EndClassify();

        // Second pass, test overlaps
        /* For each face, we test its barycenter and test all potential overlaps for the opposite face.
           - Opposite faces have EXACTLY the same barycenter, which ensures classification in the SAME cell
           - However, we DO NOT compare midpoints, we explicitly compare face VID to identify opposite faces.
        */
        uint32 num_tests_sh(0);
        ssh.BeginTestAdd();
        {
            std::vector<tetsolid3_face_id_type> vec_potential_overlaps;
            for( unsigned int it_tet=0; it_tet<m_addT.size(); it_tet++ )
            {
                editable_tetrahedron_type &tet( m_addT[it_tet] );
                // GEO_LOG_WARNING( "RebuildTopology() for tet[%d] with vid = (%d,%d,%d,%d) and ntit = (%d,%d,%d,%d)",
                //                  it_tet,
                //                  tet.m_vecVID[0], tet.m_vecVID[1], tet.m_vecVID[2], tet.m_vecVID[3],
                //                  tet.m_vecNTID[0], tet.m_vecNTID[1], tet.m_vecNTID[2], tet.m_vecNTID[3] );
                for( unsigned int it_neighbour=0; it_neighbour<4; it_neighbour++ )
                {
                    if( ssh.TestAdd( std::make_pair( it_tet, it_neighbour ),
                                     mal::Rcp<Real>(3.0f)*( m_addV[ tet.m_vecVID[(it_neighbour+1)%4] ].m_Pos
                                                            + m_addV[ tet.m_vecVID[(it_neighbour+2)%4] ].m_Pos
                                                            + m_addV[ tet.m_vecVID[(it_neighbour+3)%4] ].m_Pos ),
                                     GEO_TETSS_EPSILON_LENGTH,  //GEO_TETSS_EPSILON_LENGTH_SQ, //\todo CRITICAL: this SHOULD BE GEO_TETSS_EPSILON_LENGTH, not squared... MANY topology problems could derive from this shit!!
                                     vec_potential_overlaps ) )
                    {
                        // GEO_LOG_WARNING( "RebuildTopology() for tet[%d].ntit[%d] has %d potential matches", it_tet, it_neighbour, (int)vec_potential_overlaps.size() );
                        // Given face (it_tet,it_neighbour) check if any EXISTING potential overlap is the opposite face
                        //\note That (it_tet,it_neighbour) is *NEVER* in vec_potential_overlaps (as it's incrementally added in Test())
                        for( uint32 it_pov=0; it_pov < vec_potential_overlaps.size(); it_pov++ )
                        {
                            tetsolid3_face_id_type nfid( vec_potential_overlaps[it_pov] );
                            // count shared vid with potential overlapping face, 3 => neighbour
                            uint32 num_shared_vid(0);
                            for( unsigned int it_vif=0; it_vif<3; it_vif++ )
                            {
                                uint32 vid( tet.m_vecVID[(it_neighbour+it_vif+1)%4] ); //(+1,+2,+3) %4 are the vif of the opposite face to the vertex it_neighbour
                                bool b_shared(false);
                                for( unsigned int it_nvif=0; it_nvif<3 && !b_shared; it_nvif++ )
                                    b_shared = vid == m_addT[nfid.first].m_vecVID[(nfid.second+it_nvif+1)%4];
                                if( b_shared )
                                    num_shared_vid++;
                            }
                            if( num_shared_vid == 3 )
                            {
                                tet.m_vecNTID[it_neighbour] = nfid.first;
                                m_addT[nfid.first].m_vecNTID[nfid.second] = it_tet;
                                // GEO_LOG_WARNING( "RebuildTopology() ACCEPTING tet[%d].ntit[%d] = %d with vid = (%d,%d,%d,%d)",
                                //                  it_tet, it_neighbour, nfid.first,
                                //                  m_addT[nfid.first].m_vecVID[0], m_addT[nfid.first].m_vecVID[1], m_addT[nfid.first].m_vecVID[2], m_addT[nfid.first].m_vecVID[3] );
                            }
                            else
                            {
                                // GEO_LOG_WARNING( "RebuildTopology() DISCARDING tet[%d].ntit[%d] = %d with vid = (%d,%d,%d,%d)",
                                //                  it_tet, it_neighbour, nfid.first,
                                //                  m_addT[nfid.first].m_vecVID[0], m_addT[nfid.first].m_vecVID[1], m_addT[nfid.first].m_vecVID[2], m_addT[nfid.first].m_vecVID[3] );
                            }
                            num_tests_sh++;
                        }
                    }
                }
            }
        }
        ssh.EndTestAdd();
        // GEO_LOG( "EditableTetSolidShape3::RebuildTopology() SH tests = %d << %d", num_tests_sh, (num_faces * (num_faces-1)) / 2 );
    }

    // Rebuild BF from Tet topology
    m_addBF.clear();
    uint32 num_bf(0);
    for( unsigned int it_tet=0; it_tet<m_addT.size(); it_tet++ )
        for( unsigned int it_neighbour=0; it_neighbour<4; it_neighbour++ )
            if( cTetSolid3_InvalidFeatureIndex == m_addT[it_tet].m_vecNTID[it_neighbour] )
            {
                num_bf++;
                uint32 vec_vid[3] = { m_addT[it_tet].m_vecVID[(it_neighbour+1)%4],
                                      m_addT[it_tet].m_vecVID[(it_neighbour+2)%4],
                                      m_addT[it_tet].m_vecVID[(it_neighbour+3)%4] };
                if( mal::Dot( m_addV[ m_addT[it_tet].m_vecVID[it_neighbour] ].m_Pos - m_addV[vec_vid[0]].m_Pos,
                              mal::Cross(m_addV[vec_vid[1]].m_Pos-m_addV[vec_vid[0]].m_Pos,
                                         m_addV[vec_vid[2]].m_Pos-m_addV[vec_vid[0]].m_Pos) ) <= 0 )
                    m_addBF.push_back( editable_boundary_face_type( vec_vid[0], vec_vid[1], vec_vid[2] ) );
                else
                    m_addBF.push_back( editable_boundary_face_type( vec_vid[0], vec_vid[2], vec_vid[1] ) );
            }
    // GEO_LOG("EditableTetSolidShape3::RebuildTopology() has %d BF", num_bf );

    //\todo Compute topology stuff (Connected,SimpleConnected,Closed...)
}

/* Rebuild layers from m_addX data, assuming Topology is up-to-date
 */
void EditableTetSolidShape3::RebuildLayers()
{
    const unsigned int cMaxLayers = m_addT.size();
    //1) Compute per-vertex layer from boundary
    std::vector<unsigned int> vecVertexLayer( m_addV.size(), cMaxLayers );
    // Set boundary V layers to 0
    for( unsigned int it_bf=0; it_bf<m_addBF.size(); it_bf++ )
        for( unsigned int it_vibf=0; it_vibf<3; it_vibf++ )
            vecVertexLayer[ m_addBF[it_bf].m_vecVID[it_vibf] ] = 0;
    // Propagate layer from boundary (\todo O(MaxLayers*E) iteration could be optimized by O(V+E) flooding)
    for( unsigned int i=0; i<cMaxLayers; i++ )
        for( unsigned int it_tet=0; it_tet<m_addT.size(); it_tet++ )
            for( int it_vit1=0; it_vit1<4; it_vit1++ )
                for( int it_vit2=0; it_vit2<4; it_vit2++ )
                    vecVertexLayer[ m_addT[it_tet].m_vecVID[it_vit1] ] = mal::Min( vecVertexLayer[ m_addT[it_tet].m_vecVID[it_vit1] ],
                                                                                   vecVertexLayer[ m_addT[it_tet].m_vecVID[it_vit2] ]+1 );
    //2) Compute per-tetrahedron layer from min adjacent vertex layer
    std::vector<unsigned int> vecTetrahedronLayer( m_addT.size(), cMaxLayers );
    for( unsigned int it_tet=0; it_tet<m_addT.size(); it_tet++ )
        for( int it_vit=0; it_vit<4; it_vit++ )
            vecTetrahedronLayer[it_tet] = mal::Min( vecTetrahedronLayer[it_tet], vecVertexLayer[ m_addT[it_tet].m_vecVID[it_vit] ] );
    //3) Sort tetrahedrons per layer and generate array of layer descriptors [0..num_layers-1] \todo BOUNDARY P do not change!!
    unsigned int max_layer(0);
    for( unsigned int it_tet=0; it_tet<m_addT.size(); it_tet++ )
        if( max_layer < vecTetrahedronLayer[it_tet] )
            max_layer = vecTetrahedronLayer[it_tet];
    unsigned int num_layers(max_layer+1);
    // Add Tetrahedrons to layer-lists
    std::vector< std::vector< tetsolid3_feature_index_type > > vecLayers( num_layers );
    for( unsigned int it_tet=0; it_tet<m_addT.size(); it_tet++ )
        vecLayers[ vecTetrahedronLayer[it_tet] ].push_back( it_tet );
    // Build contiguous layers and compute TID mapping
    m_addL.resize( num_layers ); //sets all layers to 0 elements
    std::vector< tetsolid3_feature_index_type > vecOld2NewTID( m_addT.size() );
    std::vector< editable_tetrahedron_type > addT_sorted;
    unsigned int last_tid(0);
    for( unsigned int it_l=0; it_l<vecLayers.size(); it_l++ )
    {
        // GEO_LOG_WARNING( "Layer[%d] F%d,S%ld", it_l, last_tid, vecLayers[it_l].size() );
        m_addL[ it_l ].m_FirstTID = last_tid;
        m_addL[ it_l ].m_NumTetrahedrons = vecLayers[it_l].size();
        for( unsigned int it_til=0; it_til<vecLayers[it_l].size(); it_til++ )
        {
            vecOld2NewTID[ vecLayers[it_l][it_til] ] = last_tid++;
            addT_sorted.push_back( m_addT[ vecLayers[it_l][it_til] ] );
            // GEO_LOG_WARNING( "%d => %d", vecLayers[it_l][it_til], last_tid-1 );
        }
    }
    // Swap unsorted with sorted
    std::swap( m_addT, addT_sorted );
    //4) Remap NTID: No need to rebuild full topology, just reindex NTID
    for( unsigned int it_tet=0; it_tet<m_addT.size(); it_tet++ )
        for( int it_ntit=0; it_ntit<4; it_ntit++ )
            if( m_addT[it_tet].m_vecNTID[it_ntit] != cTetSolid3_InvalidFeatureIndex )
                m_addT[it_tet].m_vecNTID[it_ntit] = vecOld2NewTID[ m_addT[it_tet].m_vecNTID[it_ntit] ];
}

bool EditableTetSolidShape3::FixDegenerateVertices()
{
    //\todo Weld coincident vertices, see EditableTriSurfaceShape3::FixDegenerateVertices()
    //IMPORTANT: This WILL BREAK FixNonManifoldFeatures() if split V are not geometrically separated to be > epsilon_length

    // Compute V degree (topology is NOT available, only m_addT data is up to date)
    std::vector< int > vecVertexDegree( m_addV.size(), 0 );
    for( unsigned int it_tet=0; it_tet<m_addT.size(); it_tet++ )
        for( int it_vit=0; it_vit<4; it_vit++ )
            vecVertexDegree[ m_addT[it_tet].m_vecVID[it_vit] ]++;
    // Remove V with 0-degree and remap to their new index
    std::vector< tetsolid3_feature_index_type > vecOld2NewVID( m_addV.size() );
    unsigned int num_removed(0);
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
    {
        if( vecVertexDegree[it_v] > 0 )
        {
            vecOld2NewVID[it_v] = it_v - num_removed;
            m_addV[it_v-num_removed] = m_addV[it_v];
        }
        else
        {
            vecOld2NewVID[it_v] = cTetSolid3_InvalidFeatureIndex;
            num_removed++;
        }
    }
    m_addV.resize( m_addV.size() - num_removed ); //crop last num_removed V
    // Remap T vid
    for( unsigned int it_tet=0; it_tet<m_addT.size(); it_tet++ )
        for( int it_vit=0; it_vit<4; it_vit++ )
            m_addT[it_tet].m_vecVID[it_vit] = vecOld2NewVID[ m_addT[it_tet].m_vecVID[it_vit] ];
    // Remap BF vid => Nothing to do, BF are computed AFTERWARDS in RebuildTopology()

    if( num_removed > 0 )
        GEO_LOG( "EditableTetSolidShape3::FixDegenerateVertices() removed %u", num_removed );
    return num_removed > 0;
}

bool EditableTetSolidShape3::FixDegenerateTetrahedrons() { return false; }


//-----------------------------------------------------------------------------
// Free functions
//-----------------------------------------------------------------------------
Real ComputeVolume( const TetSolidShape3& tss3, const Transform3& transform, const Vec3* vec_sdof )
{
    const Vec3* actual_sdof( ( 0 != vec_sdof ) ? vec_sdof : tss3.GetVecDefaultSDOF() );
    Real volume_6(0);
#ifdef __DISABLED_CODE
    for( unsigned int it_tet=0; it_tet<tss3.GetNumT(); it_tet++ )
    {
        Vec3 vec_tet_pos[4] = { tss3.V_Pos( tss3.T_VID(it_tet,0), actual_sdof ),
                                tss3.V_Pos( tss3.T_VID(it_tet,1), actual_sdof ),
                                tss3.V_Pos( tss3.T_VID(it_tet,2), actual_sdof ),
                                tss3.V_Pos( tss3.T_VID(it_tet,3), actual_sdof ) };
        volume_6 += mal::Det(mal::GMat3x3_From_Columns( vec_tet_pos[1]-vec_tet_pos[0],
                                                        vec_tet_pos[2]-vec_tet_pos[0],
                                                        vec_tet_pos[3]-vec_tet_pos[0] ));
    }
#else // faster version using BF
    for( unsigned int it_bf=0; it_bf<tss3.GetNumBF(); it_bf++ )
        volume_6 += mal::Det(mal::GMat3x3_From_Columns( tss3.V_Pos( tss3.BF_VID(it_bf,0), actual_sdof ),
                                                        tss3.V_Pos( tss3.BF_VID(it_bf,1), actual_sdof ),
                                                        tss3.V_Pos( tss3.BF_VID(it_bf,2), actual_sdof ) ) );
#endif
    return volume_6/6;
}

//-----------------------------------------------------------------------------
//---- Make_TetSolidShape3_XXX Implementation
//-----------------------------------------------------------------------------

/* Minimal cube tetrahedralization => 5 tet
   \see "Simplicial Subdivisions and Sampling Artifacts"

   - CONVENTION: vid 0..3 and 4..7 represent the back/front z-parallel
     quads, each one with clockwise (as seen looking along -Z
     direction) enumerated vertices, starting from the -x,-y corner.
*/
void Make_TetSolidShape3_Box( EditableTetSolidShape3 &etss3, const Vec3 &half_sizes )
{
    etss3.Clear();
    etss3.BeginEdition();
    {
        tetsolid3_feature_index_type vid0 = etss3.AddVertex( Vec3( -half_sizes.x(), -half_sizes.y(), -half_sizes.z() ) );
        tetsolid3_feature_index_type vid1 = etss3.AddVertex( Vec3( -half_sizes.x(),  half_sizes.y(), -half_sizes.z() ) );
        tetsolid3_feature_index_type vid2 = etss3.AddVertex( Vec3(  half_sizes.x(),  half_sizes.y(), -half_sizes.z() ) );
        tetsolid3_feature_index_type vid3 = etss3.AddVertex( Vec3(  half_sizes.x(), -half_sizes.y(), -half_sizes.z() ) );
        tetsolid3_feature_index_type vid4 = etss3.AddVertex( Vec3( -half_sizes.x(), -half_sizes.y(),  half_sizes.z() ) );
        tetsolid3_feature_index_type vid5 = etss3.AddVertex( Vec3( -half_sizes.x(),  half_sizes.y(),  half_sizes.z() ) );
        tetsolid3_feature_index_type vid6 = etss3.AddVertex( Vec3(  half_sizes.x(),  half_sizes.y(),  half_sizes.z() ) );
        tetsolid3_feature_index_type vid7 = etss3.AddVertex( Vec3(  half_sizes.x(), -half_sizes.y(),  half_sizes.z() ) );
        etss3.AddTetrahedron( vid0, vid1, vid4, vid3 ); //Around v0
        etss3.AddTetrahedron( vid2, vid1, vid3, vid6 ); //Around v2
        etss3.AddTetrahedron( vid5, vid1, vid6, vid4 ); //Around v5
        etss3.AddTetrahedron( vid6, vid1, vid3, vid4 ); //Around v6
        etss3.AddTetrahedron( vid7, vid3, vid4, vid6 ); //Around v7
    }
    etss3.EndEdition();
}

/* Box tetrahedralization with given resolution
   - Each box element is split into 5 tets as in the previous
     function, but sub-box faces need to be swapped at odd indices for
     their global tetraedralization to match seamlessly.
     - If swaps are not performed, adjacent sub-boxes have incompatible
       tetraedralizations (actually, triangulations of their shared
       quad.
   \see "Simplicial Subdivisions and Sampling Artifacts", the "parity"
        rule mentioned there seems to guarantee continuity, as I've
        deduced.
*/
void Make_TetSolidShape3_Box( EditableTetSolidShape3 &etss3, const Vec3 &half_sizes, unsigned int num_x, unsigned int num_y, unsigned int num_z )
{
    GEO_ASSERT( !mal::IsNaN(half_sizes) );
    GEO_ASSERT( num_x > 1 && num_y > 1 && num_z > 1 );
    Vec3 sizes( 2*half_sizes );
    etss3.Clear();
    etss3.BeginEdition();
    {
        // Add vertices, linear vid(i,j,k) = [i*dimensions[1]*dimensions[2] + j*dimensions[2] + k]
        for( unsigned int i=0; i<num_x; i++ )
            for( unsigned int j=0; j<num_y; j++ )
                for( unsigned int k=0; k<num_z; k++ )
                    etss3.AddVertex( -half_sizes + geo::Vec3( sizes[0] * float(i)/(num_x-1),
                                                              sizes[1] * float(j)/(num_y-1),
                                                              sizes[2] * float(k)/(num_z-1) ) );
        // For each cube (dimensions-1), add 5 tets, swapping cube faces along dimensions with odd indices
        for( unsigned int i=0; i<num_x-1; i++ )
            for( unsigned int j=0; j<num_y-1; j++ )
                for( unsigned int k=0; k<num_z-1; k++ )
                {
                    uint32 vid0 = i*num_y*num_z + j*num_z + k; //---
                    uint32 vid1 = i*num_y*num_z + (j+1)*num_z + k; //-+-
                    uint32 vid2 = (i+1)*num_y*num_z + (j+1)*num_z + k; //++-
                    uint32 vid3 = (i+1)*num_y*num_z + j*num_z + k; //+--
                    uint32 vid4 = i*num_y*num_z + j*num_z + (k+1); //--+
                    uint32 vid5 = i*num_y*num_z + (j+1)*num_z + (k+1); //-++
                    uint32 vid6 = (i+1)*num_y*num_z + (j+1)*num_z + (k+1); //+++
                    uint32 vid7 = (i+1)*num_y*num_z + j*num_z + (k+1); //+-+
                    // Swap Quads
                    int num_mirrors(0);
                    if( i%2 == 1 )
                    {
                        num_mirrors++;
                        std::swap( vid0, vid3 );
                        std::swap( vid1, vid2 );
                        std::swap( vid5, vid6 );
                        std::swap( vid4, vid7 );
                    }
                    if( j%2 == 1 )
                    {
                        num_mirrors++;
                        std::swap( vid0, vid1 );
                        std::swap( vid3, vid2 );
                        std::swap( vid4, vid5 );
                        std::swap( vid7, vid6 );
                    }
                    if( k%2 == 1 )
                    {
                        num_mirrors++;
                        std::swap( vid0, vid4 );
                        std::swap( vid1, vid5 );
                        std::swap( vid2, vid6 );
                        std::swap( vid3, vid7 );
                    }
                    // If #mirrors is odd, the tets are inverted
                    // (det<0), so we swap the last 2vid \todo WE
                    // COULD FIX THIS inside AddTetrahedron, as tets
                    // do NOT HAVE, by now, a canonical "orientation"
                    // as triangles... but this may change if some
                    // reasonable order is discovered, anyway...
                    if( num_mirrors%2 == 0 )
                    {
                        etss3.AddTetrahedron( vid0, vid1, vid4, vid3 ); //Around v0
                        etss3.AddTetrahedron( vid2, vid1, vid3, vid6 ); //Around v2
                        etss3.AddTetrahedron( vid5, vid1, vid6, vid4 ); //Around v5
                        etss3.AddTetrahedron( vid6, vid1, vid3, vid4 ); //Around v6
                        etss3.AddTetrahedron( vid7, vid3, vid4, vid6 ); //Around v7
                    }
                    else
                    {
                        etss3.AddTetrahedron( vid0, vid1, vid3, vid4 ); //Around v0
                        etss3.AddTetrahedron( vid2, vid1, vid6, vid3 ); //Around v2
                        etss3.AddTetrahedron( vid5, vid1, vid4, vid6 ); //Around v5
                        etss3.AddTetrahedron( vid6, vid1, vid4, vid3 ); //Around v6
                        etss3.AddTetrahedron( vid7, vid3, vid6, vid4 ); //Around v7
                    }
                }
    }
    etss3.EndEdition();
}

//----------------------------------------------------------------
// Create DCR \todo move elsewhere
//----------------------------------------------------------------

// per-patch and per-element data
struct dcr3_patch_id_type
{
    uint32 m_EID; //Element id
    int m_PIE; //Patch-in-element
    dcr3_patch_id_type() : m_EID(0xFFFFFFFF), m_PIE(-1) {}
    dcr3_patch_id_type( uint32 eid, int pie ) : m_EID(eid), m_PIE(pie) {}
};
struct dcr3_editable_patch_type
{
    dcr3_editable_patch_type() {}
    std::vector<trisurface3_feature_index_type> m_vecTID;
    std::vector<trisurface3_edge_id_type> m_vecBEID;
    std::vector<trisurface3_feature_index_type> m_vecBVID;
};
struct dcr3_editable_element_type
{
    dcr3_editable_element_type() {}
    std::vector< dcr3_editable_patch_type > m_vecP;
};

/* Build the DCR for a TriSurfaceShape3 embedded in a TetSolidShape3
   \pre The surface CAN be open
   \pre The surface MUST be completely enclosed in the TetSS Layer[0]
   \pre The surface MAY NOT be clipped against the TetSS

   \note Non-clipped surface: As we classify triangle barycenters into
   a single DCR.E, and therefore this SHOULD WORK even for non-clipped
   surfaces. This may be useful if clipping results too expensive and
   we can assume small surface inaccuracies across DCR.E borders,
   which can be reasonable for highly detailed meshes). HOWEVER, the
   DCR.E will NOT be guaranteed to strictly contain all classified
   geometry. To do so, we should enlarge the BV(DCR.E) in a
   non-trivial way (it will depend on neighbour DCR.E configuration
   and the actual deformation method Barycentric/MLS/SubdK...

   Per-element partition Vs cross-patch/element topology can be
   handled in different ways:
   a) If we split the V in each patch boundary, then cross-patch
      topology may remain consistent regarding neighbour TID, BUT
      cross-patch neighbour triangles will NOT share global VID, as
      they've been split at the border to guarantee partition.
      - Splitting V will INCREASE once again the #V, which may have
        already been increased at Clip...
      - Splitting V in the non-clipped surface case is DANGEROUS, as
        the split-V will be EXTERNAL to one of the DCR.E that contain
        DCR.T that refer any of the split-V copies. EACH of these
        split-V copies will deform DIFFERENTLY, according to its
        containing DCR.E... which would cause
        cracks/self-intersections. The good news is that this split-V
        COULD be stored in a different sub-array per-DCR.E and
        transformed using linear-blend-skinning using neighbour DCR.E
        barycentric transforms... however, this ENFORCES a specific
        deformation method on the DCR, while the Clip-only version
        would not...
   b) Alternatively, we can avoid splitting V at cross-element
      boundaries, but store each of these cross-patch boundary curves
      explicitly. Boundary curves would then index global boundary VID,
      which WOULD NOT be part of any patch/element, but stored
      separatedly in a boundary-curves array. During collision
      detection, when a DCR.E is overlapped and needs to be "detailed"
      (transform V, sample points, etc...), ALL DCR.E boundary curves
      should also be considered.
      - However, cross-element VID cannot NOT be stored in a single
        subarray per-element either, only in patch-pair subarrays
        (whicn can be shared if reversed)
      - This solution requires adding the "boundary curve" concept to
        the DCR CD algorithms (eg: avoid re-transforming boundary
        curve vertices when both DCR.E sharing it are already
        transformed)

    IMPORTANT: By now, we choose option a) without any specific
    handling for the non-Clipped case (thus cracks will appear). We'll
    allocate split-V at the END of per-DCR.E vertex subarrays and, in
    the future, we'll consider handling them specifically. We will
    also IGNORE the potential errors due to the lack of strict
    BV(DCR.E) containment guarantee.

    BUG DCR.E > SIM.layer[0]
      En 3D, TOT I QUE m'asseguro d'eliminar tots els SIM.E que són
    completament exteriors a la Surface CLP, per problemes de
    precisió, crec, acaba passant que un cop clipejada la CLP contra
    el Solid:
    - Alguns SIM.E de la layer[0] NO CONTENEN (el baricentre de) cap CLP.T
      => Aquestos SIM.E HAURIEN d'haver estat etiquetats com a
         EXTERIORS i eliminats prèviament...
    - Com a CONSEQÜENCIA, alguns SIM.E que no són a la layer[0]
      CONTENEN (de forma exacta, perfect match) algun (bary) CLP.T i,
      per tant, acaben esdevenint DCR.E, de manera que #DCR.E >
      #layer[0]
    Solucions:
    - Fer passada pre-DCR, post-CLP per eliminar els SIM.Layer[0].E
      que no contenen cap CLP.T
      => aixo requereix conèixer CLP, per tant no es pot fer en la
         creació/fitting de SIM, que és on s'eliminen els SIM.E exteriors
    - Si arran d'això apareixen CLP.T sense perfect match, triar el
      closest SIM.E, pero PROJECTAR els vtx del CLP.T per tal que el
      DCR.F/V corresponent SI QUE estigui dins del closest SIM.E.

    BUG:
    - Found triangles with VID outside their containing element vecV[First,First+NumVertices]
    - NPID not yet available
*/
DCR_TetSolidShape3* Create_DCR_TetSolidShape3_From_TriSurfaceShape3( const TetSolidShape3& solid, const Transform3& solid_tr, const Vec3* vec_solid_sdof,
                                                                     const TriSurfaceShape3& surface, const Transform3& surface_tr, const Vec3* vec_surface_sdof )
{
    // Build a local BVH with solid's given transform and DOF
    // \todo Consider requiring valid and up-to-date solid.GetBVH() instead
    BVH_TetSolidShape3 bvh;
    bvh.Rebuild_TopDown( solid.GetNumT(), //TEMP solid.L_NumT(0), \todo this SHOULD WORK, but leaves plenty of surface.T without a perfect match solid.T...
                         boost::bind<void>( &GEBV_TetSolidShape3_E<BVH_TetSolidShape3::entry_index_type,BVH_TetSolidShape3::bv_type>,
                                            &solid, solid_tr, vec_solid_sdof,
                                            _1, _2) );

    Transform3 inv_solid_tr(mal::Inverse(solid_tr));
    // All solid.E are initially empty, all surface.T unflooded
    std::vector<dcr3_editable_element_type> vecElements( solid.GetNumT() );
    std::vector<bool> vecIsOpenTri( surface.GetNumT(), true );

    // Classify surface.T into solid.E using T.barycenter and solid.BVH, add to DCR.E lists
    uint32 max_tet_id(0);
    for( unsigned int it_tri=0; it_tri<surface.GetNumT(); it_tri++ )
    {
        // In the surface.T is open, it'll be the seed of a new patch DCR.P
        if( vecIsOpenTri[it_tri] )
        {
            // Find solid.E containing surface.T
            std::vector< BVH_TetSolidShape3::entry_index_type > vecOverlaps;
            Vec3 tri_barycenter( surface_tr*surface.T_Barycenter(it_tri,vec_surface_sdof) );
            BVH_TetSolidShape3::bv_type tri_bv(tri_barycenter);
            tri_bv.Extend(0.01f);
            uint32 overlap_tet_id( cTetSolid3_InvalidFeatureIndex );
            if( bvh.Test( tri_bv, vecOverlaps ) )
            {
                uint32 closest_tet_id( cTetSolid3_InvalidFeatureIndex );
                Real closest_barycentric_distance_sq( mal::Infinity<Real>() );
                uint32 num_perfect_matches(0);
                for( unsigned int it_o=0;
                     it_o<vecOverlaps.size(); //TEMP: Trying to fix #DCR.E > #layer[0]  && overlap_tet_id == cTetSolid3_InvalidFeatureIndex;
                     it_o++ )
                {
                    uint32 tet_id( vecOverlaps[it_o] );
                    /* Pure overlap requirement fails for large models (eg: dragon10K)
                    if( np::Overlap_Point3_Tetrahedron3( tri_barycenter,
                                                         solid_tr*solid.V_Pos( solid.T_VID(tet_id,0), vec_solid_sdof ),
                                                         solid_tr*solid.V_Pos( solid.T_VID(tet_id,1), vec_solid_sdof ),
                                                         solid_tr*solid.V_Pos( solid.T_VID(tet_id,2), vec_solid_sdof ),
                                                         solid_tr*solid.V_Pos( solid.T_VID(tet_id,3), vec_solid_sdof ) ) )
                        overlap_tet_id = tet_id;
                    */
                    // Compute element barycentric coordinates matrix in global coords
                    Vec3 vec_element_nodes[4] = { solid_tr*solid.V_Pos( solid.T_VID(tet_id,0), vec_solid_sdof ),
                                                  solid_tr*solid.V_Pos( solid.T_VID(tet_id,1), vec_solid_sdof ),
                                                  solid_tr*solid.V_Pos( solid.T_VID(tet_id,2), vec_solid_sdof ),
                                                  solid_tr*solid.V_Pos( solid.T_VID(tet_id,3), vec_solid_sdof ) };
                    Mat4x4 Bs( 1, 1, 1, 1,
                               vec_element_nodes[0].x(), vec_element_nodes[1].x(), vec_element_nodes[2].x(), vec_element_nodes[3].x(),
                               vec_element_nodes[0].y(), vec_element_nodes[1].y(), vec_element_nodes[2].y(), vec_element_nodes[3].y(),
                               vec_element_nodes[0].z(), vec_element_nodes[1].z(), vec_element_nodes[2].z(), vec_element_nodes[3].z() );
                    Mat4x4 invBs( mal::Inverse( Bs ) );
                    // Compute vertex barycentric coords (may be out of range)
                    Vec4 barycentric_coords( invBs * mal::Concat(1,tri_barycenter) );
                    Vec4 clamped_barycentric_coords = mal::Clamp( barycentric_coords, Vec4(0,0,0,0), Vec4(1,1,1,1) );
                    Real barycentric_distance_sq = mal::NormSq( barycentric_coords - clamped_barycentric_coords );
                    // Accept perfect match or save best candidate
                    if( barycentric_distance_sq == 0 )
                    {
                        num_perfect_matches++;
                        overlap_tet_id = mal::Min(overlap_tet_id,tet_id); //TEMP: Trying to fix #DCR.E > #layer[0] overlap_tet_id = tet_id;
                    }
                    else if( closest_tet_id == cTetSolid3_InvalidFeatureIndex
                             || barycentric_distance_sq < closest_barycentric_distance_sq )
                    {
                        closest_tet_id = tet_id;
                        closest_barycentric_distance_sq = barycentric_distance_sq;
                    }
                }
                if( overlap_tet_id == cTetSolid3_InvalidFeatureIndex )
                {
                    overlap_tet_id = closest_tet_id;
                    GEO_LOG_WARNING( "Create_DCR_TetSolidShape3_From_TriSurfaceShape3(): surface.T[%u] has no strictly overlapping solid.T, using closest solid.T[%u]",
                                     it_tri, closest_tet_id );
                }
                if( num_perfect_matches > 1 )
                    GEO_LOG_WARNING( "Create_DCR_TetSolidShape3_From_TriSurfaceShape3(): surface.T[%u] has >1 perfect matches, we selected the smallest %u",
                                     it_tri, closest_tet_id );
            }
            GEO_ASSERT(overlap_tet_id != cTetSolid3_InvalidFeatureIndex);
            max_tet_id = mal::Max( max_tet_id, overlap_tet_id );
            // Gather overlapping tet vertices
            Vec3 vec_tet_pos[4] = { solid_tr*solid.V_Pos( solid.T_VID(overlap_tet_id,0), vec_solid_sdof ),
                                    solid_tr*solid.V_Pos( solid.T_VID(overlap_tet_id,1), vec_solid_sdof ),
                                    solid_tr*solid.V_Pos( solid.T_VID(overlap_tet_id,2), vec_solid_sdof ),
                                    solid_tr*solid.V_Pos( solid.T_VID(overlap_tet_id,3), vec_solid_sdof ) };
            // Create new DCR.P
            vecElements[overlap_tet_id].m_vecP.push_back( dcr3_editable_patch_type() );
            dcr3_editable_patch_type& ep( vecElements[overlap_tet_id].m_vecP.back() );
            // Flood from surface.T while inside solid.E
            std::vector<trisurface3_feature_index_type> stackTID;
            stackTID.push_back( it_tri );
            do
            {
                // pop triangle
                trisurface3_feature_index_type tid( stackTID.back() );
                stackTID.pop_back();
                // process if still open, as it may have been closed through a different route while on stack
                if( vecIsOpenTri[tid] )
                {
                    // add tri to ep and close it
                    ep.m_vecTID.push_back( tid );
                    vecIsOpenTri[tid] = false;
                    // flood to existing neighbours if inside same tet
                    for( unsigned int it_ntit=0; it_ntit<3; it_ntit++ )
                    {
                        uint32 ntid( surface.T_NTID(tid,it_ntit) );
                        if( ntid != cTriSurface3_InvalidFeatureIndex )
                        {
                            if( vecIsOpenTri[ntid] )
                            {
                                // Flood if inside DCR.E
                                if( np::Overlap_Point3_Tetrahedron3( surface_tr*surface.T_Barycenter(ntid,vec_surface_sdof),
                                                                     vec_tet_pos[0], vec_tet_pos[1], vec_tet_pos[2], vec_tet_pos[3] ) )
                                    stackTID.push_back( ntid ); //\note it may be already closed, but we'll NEED to re-check when processed anyway
                                else // element boundary otherwise
                                    ep.m_vecBEID.push_back( trisurface3_edge_id_type(tid,it_ntit) );
                            }
                            else if( ep.m_vecTID.end() == std::find( ep.m_vecTID.begin(), ep.m_vecTID.end(), ntid ) ) //closed and NOT in the same DCR.P ==> boundary
                                ep.m_vecBEID.push_back( trisurface3_edge_id_type(tid,it_ntit) );

                            /* This version adds BE when trying to re-open a closed ntid REGARDLESS of its containing DCR.P.
                               It SEEMED to work, but it created TOO MANY PATCHES and VERTICES due to artificial boundaries
                            //if open neighbour inside, flood
                            if( vecIsOpenTri[ntid]
                                && np::Overlap_Point3_Tetrahedron3( surface_tr*surface.T_Barycenter(ntid,vec_surface_sdof),
                                                                    vec_tet_pos[0], vec_tet_pos[1], vec_tet_pos[2], vec_tet_pos[3] ) )
                                stackTID.push_back( ntid ); //\note it may be already closed, but we'll NEED to re-check when processed anyway
                            else // element boundary
                                ep.m_vecBEID.push_back( trisurface3_edge_id_type(tid,it_ntit) );
                            */

                            /* BUG!! Due to numeric fuckups,
                               Overlap_Point3_Tetrahedron3() may be
                               inconsistent and DCR.P with a SINGLE
                               triangle have been observed to add 1 or
                               2 BE instead of 3!!

                               => THIS is probably the cause of
                               out-of-element-range VID in triangles...

                               WE COULD check if neighbour is closed
                               and, if so, check if it's in THIS
                               patch, and if NOT, then it must be in a
                               different previous one, thus the shared
                               edge MUST be a boundary, regardless of
                               the results of
                               Overlap_Point3_Tetrahedron3()
                            // if neighbour inside, flood
                            if( np::Overlap_Point3_Tetrahedron3( surface_tr*surface.T_Barycenter(ntid,vec_surface_sdof),
                                                                 vec_tet_pos[0], vec_tet_pos[1], vec_tet_pos[2], vec_tet_pos[3] ) )
                                stackTID.push_back( ntid ); //\note it may be already closed, but we'll NEED to re-check when processed anyway
                            else // element boundary
                                ep.m_vecBEID.push_back( trisurface3_edge_id_type(tid,it_ntit) );
                            */
                        }
                        else // global boundary
                            ep.m_vecBEID.push_back( trisurface3_edge_id_type(tid,it_ntit) );
                    }
                }
            }
            while( !stackTID.empty() );
        }
    }
    // Split shared boundary edges to guarantee self-contained DCR.E->P->T->V
    unsigned int num_patches(0);
    unsigned int num_split_e(0);
    for( unsigned int it_e=0; it_e<vecElements.size(); it_e++ )
    {
        //\todo multi-patch solid.E SHOULD be reported or, better, split (if solid was non-const...)
        dcr3_editable_element_type& element( vecElements[it_e] );
        num_patches += element.m_vecP.size();
        if( element.m_vecP.size() > 1 )
            GEO_LOG_WARNING("Element[%u] has %d > 1 patches", it_e, (int)element.m_vecP.size() );
        // Split boundary-edge origin VID for all inter-patch boundary edges (NOT for global boundary edges)
        for( auto patch : element.m_vecP )
            for( auto beid : patch.m_vecBEID )
                if( surface.T_NTID( beid.first, beid.second ) != cTriSurface3_InvalidFeatureIndex )
                {
                    patch.m_vecBVID.push_back( surface.T_VID(beid.first,beid.second) );
                    num_split_e++;
                }
        //TEMP: Checking single-face DCR.P have correct number of BE
        for( unsigned int it_p=0; it_p<element.m_vecP.size(); it_p++ )
            GEO_LOG_ERROR_IF( element.m_vecP[it_p].m_vecTID.size() == 1 && element.m_vecP[it_p].m_vecBEID.size() != 3,
                              "DCR.E[%u].P[%u] has 1 Triangle but %u != 3 Boundary Edges!",
                              it_e, it_p, (uint32)element.m_vecP[it_p].m_vecBEID.size() );
    }

    // LOG
    // {
    //     for( unsigned int it_e=0; it_e<max_tet_id+1; it_e++ )
    //     {
    //         const dcr3_editable_element_type& ee( vecElements[it_e] );
    //         GEO_LOG("dcr3_editable_element_type[%u]: #P = %u",
    //                 it_e, (uint32)ee.m_vecP.size() );
    //         for( unsigned int it_p=0; it_p<ee.m_vecP.size(); it_p++ )
    //         {
    //             const dcr3_editable_patch_type& ep( ee.m_vecP[it_p] );
    //             GEO_LOG("dcr3_editable_patch_type[%u]: #T = %u, #BE = %u",
    //                     it_p, (uint32)ep.m_vecTID.size(), (uint32)ep.m_vecBEID.size() );
    //             for( unsigned int it_t=0; it_t<ep.m_vecTID.size(); it_t++ )
    //             {
    //                 uint32 tid( ep.m_vecTID[it_t] );
    //                 GEO_LOG("T[%u] = %u", it_t, tid );
    //             }
    //         }
    //     }
    // }

    // \todo HERE we could reorder BoundaryEdges CCW so that their
    // splitVID are correlative and allow fast access to the geometric
    // DCR.P boundary. This would be useful for VIZ and have no
    // runtime impact

    // Bake the DCR
    DCR_TetSolidShape3* pDCR = new DCR_TetSolidShape3();
    // \todo Ideally only layer[0] elements contain patches, but we'll
    // tolerate otherwise by now
    pDCR->m_NumElements = max_tet_id+1;
    if( pDCR->m_NumElements != solid.L_NumT(0) )
        GEO_LOG_ERROR( "Create_DCR_TetSolidShape3_From_TriSurfaceShape3() DCR uses <= %u elements and TetSS layer[0] has %u",
                       pDCR->m_NumElements, solid.L_NumT(0) );
    pDCR->m_NumPatches = num_patches;
    pDCR->m_NumTriangles = surface.GetNumT();
    // Alloc DCR arrays
    pDCR->m_vecE = new DCR_TetSolidShape3::Element[ pDCR->m_NumElements ];
    pDCR->m_vecP = new DCR_TetSolidShape3::Patch[ pDCR->m_NumPatches ];
    pDCR->m_vecT = new DCR_TetSolidShape3::Triangle[ pDCR->m_NumTriangles ];

    //IMPORTANT: Vertices are NOT allocated yet, as boundary V will be relocated and potentially split
    std::vector<Vec3> addV;
    addV.reserve( surface.GetNumV() + num_split_e ); //upper bound

    /*Baking:
      for each element in TetSS that contains patches
        for each patch
           add cloned boundary vertices to global array
           add triangles to global array
             add unique non-boundary vertices to global array
           add patch
    */
    unsigned int last_patch_index(0);
    unsigned int last_triangle_index(0);
    unsigned int last_vertex_index(0);
    std::vector<trisurface3_feature_index_type> mapVID( surface.GetNumV(), cTriSurface3_InvalidFeatureIndex ); //TEMP: pinta malament... DCR.V > surface.V!!!
    std::vector<trisurface3_feature_index_type> mapTID( surface.GetNumT(), cTriSurface3_InvalidFeatureIndex );
    for( unsigned int it_e=0; it_e<pDCR->m_NumElements; it_e++ )
    {
        // GEO_LOG("it_e %u",it_e);
        const dcr3_editable_element_type& ee( vecElements[it_e] );
        DCR_TetSolidShape3::Element& ed( pDCR->m_vecE[it_e] );
        ed.m_FirstPID = last_patch_index;
        ed.m_NumPatches = 0;
        ed.m_FirstVID = last_vertex_index;
        ed.m_NumVertices = 0;
        /* Assign sequential indices to P,T and V. In case of V, as
           they are shared among T, assign them an index for each T,
           but ONLY if not yet assigned
           (==cTriSurface3_InvalidFeatureIndex), therefore we'll asign
           a unique sequential index to all V in a single pass through
           all T.
        */
        for( unsigned int it_patch=0; it_patch<ee.m_vecP.size(); it_patch++ )
        {
            // GEO_LOG("it_patch %u/%u",it_patch,last_patch_index);
            const dcr3_editable_patch_type& ep( ee.m_vecP[it_patch] );
            DCR_TetSolidShape3::Patch& pd( pDCR->m_vecP[last_patch_index] );
            pd.m_EID = it_e;
            pd.m_FirstTID = last_triangle_index;
            pd.m_NumTriangles = 0;
#define __SPLIT_BOUNDARY
#ifdef __SPLIT_BOUNDARY
            // Add ALL boundary vertices (both global and inter-patch) at the beggining of ed.m_vecV and at the end of addV/mapVID
            unsigned int last_global_vid( mapVID.size() );
            for( unsigned int it_be=0; it_be<ep.m_vecBEID.size(); it_be++ )
            {
                trisurface3_edge_id_type eid( ep.m_vecBEID[it_be] );
                uint32 vid( surface.T_VID(eid.first,eid.second) );
                addV.push_back( surface.V_Pos_0(vid) );
                mapVID.push_back(last_vertex_index);
                last_vertex_index++;
                ed.m_NumVertices++;
                // GEO_LOG("it_be %u, VID %u/%u",it_be,vid,last_vertex_index);
            }
#endif
            // Add patch triangles and vertices
            for( unsigned int it_triangle=0; it_triangle<ep.m_vecTID.size(); it_triangle++ )
            {
                uint32 tid( ep.m_vecTID[it_triangle] );
                // GEO_LOG("it_triangle %u/%u, TID %u",it_triangle,last_triangle_index,tid);
                GEO_ASSERT( tid < pDCR->m_NumTriangles );
                DCR_TetSolidShape3::Triangle& tri( pDCR->m_vecT[last_triangle_index] );
                // Copy ntid
                for( unsigned int it_ntit=0; it_ntit<3; it_ntit++ )
                    if( surface.T_NTID(tid,it_ntit) != cTriSurface3_InvalidFeatureIndex )
                        tri.SetNTID(it_ntit,surface.T_NTID(tid,it_ntit));
                    else //DCR has its SPECIFIC cInvalidTID definition (with less bits than trimesh TID!)
                        tri.SetNTID(it_ntit,DCR_TetSolidShape3::Triangle::cInvalidTID);
                // Assign linear TID and save into global array
                for( unsigned int it_vit=0; it_vit<3; it_vit++ )
                {
                    uint32 vid( surface.T_VID(tid,it_vit) );
                    // Find boundary vid
                    uint32 boundary_vid( cTriSurface3_InvalidFeatureIndex );
#ifdef __SPLIT_BOUNDARY
                    for( unsigned int it_be=0;
                         it_be<ep.m_vecBEID.size() && boundary_vid == cTriSurface3_InvalidFeatureIndex;
                         it_be++ )
                        if( surface.T_VID(ep.m_vecBEID[it_be].first,ep.m_vecBEID[it_be].second) == vid )
                            boundary_vid = last_global_vid + it_be;
#endif
                    // If not existing boundary V, add V with sequential VID
                    if( boundary_vid == cTriSurface3_InvalidFeatureIndex )
                    {
                        // GEO_LOG("it_vit %u, vid %u",it_vit,vid);
                        tri.SetVID(it_vit,vid);
                        GEO_ASSERT( vid < surface.GetNumV() );
                        //Assign linear VID and save into global array if NOT already existing
                        //GEO_ASSERT( mapVID[vid] == cTriSurface3_InvalidFeatureIndex );//MUST be a unique V
                        if( mapVID[vid] == cTriSurface3_InvalidFeatureIndex )
                        {
                            // GEO_LOG("vid/%u",last_vertex_index);
                            addV.push_back( surface.V_Pos_0(vid) );
                            mapVID[vid] = last_vertex_index;
                            last_vertex_index++;
                            ed.m_NumVertices++;
                        }
                    }
                    else
                    {
                        // GEO_LOG("it_vit %u, bvid %u",it_vit,boundary_vid);
                        tri.SetVID(it_vit,boundary_vid);
                        //\todo raises, but result seems ok... GEO_ASSERT( boundary_vid < addV.size() );
                    }
                }
                mapTID[tid] = last_triangle_index;
                last_triangle_index++;
                pd.m_NumTriangles++;
            }
            //\note Patch topology added in a later pass
            last_patch_index++;
            ed.m_NumPatches++;
        }
    }
    GEO_ASSERT( last_patch_index == pDCR->m_NumPatches );
    GEO_ASSERT( last_triangle_index == pDCR->m_NumTriangles );
    GEO_ASSERT( last_vertex_index == addV.size() );

    // Save addV into m_vecV
    pDCR->m_NumVertices = addV.size();
    pDCR->m_vecV = new Vec3[ pDCR->m_NumVertices ];
    memcpy( pDCR->m_vecV, &addV[0], pDCR->m_NumVertices*sizeof(Vec3) );

    // Reindex DCR.T.VID and DCR.T.NTID using mapVID and mapTID
    for( unsigned int it_tri=0; it_tri<pDCR->m_NumTriangles; it_tri++ )
    {
        DCR_TetSolidShape3::Triangle& tri( pDCR->m_vecT[it_tri] );
        for( unsigned int it_vit=0; it_vit<3; it_vit++ )
            tri.SetVID(it_vit,mapVID[tri.GetVID(it_vit)]);
        for( unsigned int it_ntit=0; it_ntit<3; it_ntit++ )
            if( tri.GetNTID(it_ntit) != DCR_TetSolidShape3::Triangle::cInvalidTID )
                tri.SetNTID( it_ntit, mapTID[tri.GetNTID(it_ntit)] );
    }

#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_FIXED_LENGTH
    /* Bake patch topology
       for each DCR.P.T
         for each NTID
           if NTID is in a different patch NPID
             if NPID is not already in vecNPID
                Add NPID to vecNPID
    */
    /* \todo CONSIDER Using global m_vecNPID array and a
       DCR.P.FirstIndexNPID/NumNeighbours sub-array per patch instead
       of reserving Patch::cMaxNeighbours large enough to avoid
       overflow
    */
    uint32 max_neighbours(0);
    uint32 max_triangles(0);
    for( unsigned int it_p=0; it_p<pDCR->m_NumPatches; it_p++ )
    {
        DCR_TetSolidShape3::Patch& pd( pDCR->m_vecP[it_p] );
        pd.m_NumNeighbours = 0;
        uint32 num_neighbours(0);
        for( unsigned int it_tip=0; it_tip<pd.m_NumTriangles; it_tip++ )
        {
            const DCR_TetSolidShape3::Triangle& tri( pDCR->m_vecT[pd.m_FirstTID+it_tip] );
            for( unsigned int it_ntit=0; it_ntit<3; it_ntit++ )
            {
                uint32 ntid( tri.GetNTID(it_ntit) );
                if( ntid < pd.m_FirstTID || ntid >= pd.m_FirstTID+pd.m_NumTriangles )
                {
                    uint32 npid( DCR_TetSolidShape3_Find_PID_From_TID(pDCR,ntid) );
                    GEO_ASSERT(npid != 0xFFFFFFFF);
                    bool bFound(false);
                    for( unsigned int it_npip=0; !bFound && it_npip<pd.m_NumNeighbours; it_npip++ )
                        bFound = pd.m_vecNPID[it_npip] == npid;
                    if( !bFound )
                    {
                        if( pd.m_NumNeighbours < DCR_TetSolidShape3::Patch::cMaxNeighbours )
                            pd.m_vecNPID[pd.m_NumNeighbours++] = npid;
                        num_neighbours++; //count unadded neighbours for reference
                    }
                }
            }
        }
        max_neighbours = mal::Max(num_neighbours,max_neighbours);
        max_triangles = mal::Max(pd.m_NumTriangles,max_neighbours);
    }
    GEO_LOG( "Create_DCR_TetSolidShape3_From_TriSurfaceShape3() Max DCR.P.#Neighbours = %u, DCR.P.#Triangles = %u", max_neighbours, max_triangles ); //TEMP: verbose
    GEO_LOG_ERROR_IF( max_neighbours > DCR_TetSolidShape3::Patch::cMaxNeighbours,
                      "Create_DCR_TetSolidShape3_From_TriSurfaceShape3() DCR has patches with %u > %d cMaxNeighbours",
                      max_neighbours, DCR_TetSolidShape3::Patch::cMaxNeighbours );
#endif
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
    /* Bake patch topology
       for each DCR.P.T
         for each NTID
           if NTID is in a different patch NPID
             if NPID is not already in vecNPID
                Add NPID to vecNPID
    */
    uint32 max_neighbours(0);
    uint32 max_triangles(0);
    std::vector<uint32> vecNPIP;
    for( unsigned int it_p=0; it_p<pDCR->m_NumPatches; it_p++ )
    {
        DCR_TetSolidShape3::Patch& pd( pDCR->m_vecP[it_p] );
        pd.m_FirstIndexNPID = vecNPIP.size();
        pd.m_NumNeighbours = 0;
        uint32 num_neighbours(0);
        for( unsigned int it_tip=0; it_tip<pd.m_NumTriangles; it_tip++ )
        {
            const DCR_TetSolidShape3::Triangle& tri( pDCR->m_vecT[pd.m_FirstTID+it_tip] );
            for( unsigned int it_ntit=0; it_ntit<3; it_ntit++ )
            {
                uint32 ntid( tri.GetNTID(it_ntit) );
                if( ntid < pd.m_FirstTID || ntid >= pd.m_FirstTID+pd.m_NumTriangles )
                {
                    uint32 npid( pDCR->Find_PID_From_TID(ntid) );
                    GEO_ASSERT(npid != 0xFFFFFFFF);
                    bool bFound(false);
                    for( unsigned int it_npip=0; !bFound && it_npip<pd.m_NumNeighbours; it_npip++ )
                        bFound = vecNPIP[pd.m_FirstIndexNPID+it_npip] == npid;
                    if( !bFound )
                    {
                        vecNPIP.push_back(npid);
                        pd.m_NumNeighbours++;
                    }
                }
            }
        }
        max_neighbours = mal::Max(pd.m_NumNeighbours,max_neighbours);
        max_triangles = mal::Max(pd.m_NumTriangles,max_neighbours);
    }
    pDCR->m_NumNPID = vecNPIP.size();
    pDCR->m_vecNPID = new uint32[pDCR->m_NumNPID];
    for( unsigned int it_npid=0; it_npid<pDCR->m_NumNPID; it_npid++ )
        pDCR->m_vecNPID[it_npid] = vecNPIP[it_npid];
    GEO_LOG( "Create_DCR_TetSolidShape3_From_TriSurfaceShape3() #DCR.NPID = %u. NPID/P = %f, Max DCR.P.#Neighbours = %u, DCR.P.#Triangles = %u",
             pDCR->m_NumNPID, pDCR->m_NumPatches > 0 ? float(pDCR->m_NumNPID)/pDCR->m_NumPatches : 0, max_neighbours, max_triangles );
#endif

    // Per-DCR.E pass to compute BDOP and BSlab using already flattened subarray m_vecV[ed.m_FirstVID..ed.m_NumVertices]
    for( unsigned int it_e=0; it_e<pDCR->m_NumElements; it_e++ )
    {
        DCR_TetSolidShape3::Element& ed( pDCR->m_vecE[it_e] );
        for( int i=0; i<4; i++ ) ed.m_BDOP[i] = Interval::Empty();
        ed.m_BDOP_BestSlabIdx = 0;
        // Gather element nodes
        uint32 vec_tet_vid[4] = { solid.T_VID(it_e,0), solid.T_VID(it_e,1), solid.T_VID(it_e,2), solid.T_VID(it_e,3) };
        Vec3 vec_tet_pos[4] = { solid.V_Pos_0( vec_tet_vid[0] ), solid.V_Pos_0( vec_tet_vid[1] ), solid.V_Pos_0( vec_tet_vid[2] ), solid.V_Pos_0( vec_tet_vid[3] ) };
        // Compute element barycentric coordinates matrix and its inverse
        Mat4x4 Bs( 1, 1, 1, 1,
                   vec_tet_pos[0].x(), vec_tet_pos[1].x(), vec_tet_pos[2].x(), vec_tet_pos[3].x(),
                   vec_tet_pos[0].y(), vec_tet_pos[1].y(), vec_tet_pos[2].y(), vec_tet_pos[3].y(),
                   vec_tet_pos[0].z(), vec_tet_pos[1].z(), vec_tet_pos[2].z(), vec_tet_pos[3].z() );
        Mat4x4 invBs( mal::Inverse( Bs ) );
        for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
        {
            // Compute barycentric coords
            uint32 it_v( ed.m_FirstVID + it_vie );
            Vec4 b( invBs * mal::Concat(1,pDCR->m_vecV[it_v]) );
            /*\note As we add b[i] to the BDOP, geometrically
              the BDOP.Min() plane will be the *farthest*
              from v[i], and BDOP.Max() will be the *closest*
            */
            if( it_vie == 0 ) //Init
                for( int it_axis=0; it_axis<4; it_axis++ )
                    ed.m_BDOP[it_axis].Set( b[it_axis] );
            else //Grow
                for( int it_axis=0; it_axis<4; it_axis++ )
                    ed.m_BDOP[it_axis].Merge( b[it_axis] );
        }

        // \todo WHY DOES THIS HAPPEN??? It seems that some V/T are assigned to non-layer[0] elements, which are then NOT ADDED to valid layer[0] elements that remain empty...
        if( ed.m_NumVertices == 0 )
            GEO_LOG_WARNING("Create_DCR_TetSolidShape3_From_TriSurfaceShape3() DCR.E[%u] has NO VERTICES", it_e );

        // Select best BDOP => B-Slab with smallest volume
        if( !ed.m_BDOP[0].IsEmpty()
            && !ed.m_BDOP[1].IsEmpty()
            && !ed.m_BDOP[2].IsEmpty()
            && !ed.m_BDOP[3].IsEmpty() )
        {
            geo::Vec3 bdop_pos[6];
            geo::Real volume_bdop[6] = { mal::Infinity<geo::Real>() };
            ed.m_BDOP_BestSlabIdx = 0;
            for( int i=0; i<4; i++ )
            {
                bdop_pos[0] = ed.m_BDOP[i].Min() * vec_tet_pos[i] + (1-ed.m_BDOP[i].Min()) * vec_tet_pos[(i+1)%4];
                bdop_pos[1] = ed.m_BDOP[i].Min() * vec_tet_pos[i] + (1-ed.m_BDOP[i].Min()) * vec_tet_pos[(i+2)%4];
                bdop_pos[2] = ed.m_BDOP[i].Min() * vec_tet_pos[i] + (1-ed.m_BDOP[i].Min()) * vec_tet_pos[(i+3)%4];
                bdop_pos[3] = ed.m_BDOP[i].Max() * vec_tet_pos[i] + (1-ed.m_BDOP[i].Max()) * vec_tet_pos[(i+1)%4];
                bdop_pos[4] = ed.m_BDOP[i].Max() * vec_tet_pos[i] + (1-ed.m_BDOP[i].Max()) * vec_tet_pos[(i+2)%4];
                bdop_pos[5] = ed.m_BDOP[i].Max() * vec_tet_pos[i] + (1-ed.m_BDOP[i].Max()) * vec_tet_pos[(i+3)%4];
                // Trapezoid volume = Big + 0.5*(Big-Small)
                geo::Real h = mal::Abs( mal::Dot( bdop_pos[3]-bdop_pos[0],
                                                  mal::Normalized( mal::Cross(bdop_pos[1]-bdop_pos[0],bdop_pos[2]-bdop_pos[0]) ) ) );
                Real area_big( 0.5f*mal::Norm(mal::Cross(bdop_pos[1]-bdop_pos[0],bdop_pos[2]-bdop_pos[0])) );
                Real area_small( 0.5f*mal::Norm(mal::Cross(bdop_pos[4]-bdop_pos[3],bdop_pos[5]-bdop_pos[3])) );
                if( area_big < area_small ) std::swap(area_big,area_small);
                volume_bdop[i] = h*area_big - 0.5f*h*(area_big-area_small);
                if( volume_bdop[i] < volume_bdop[ed.m_BDOP_BestSlabIdx] )
                    ed.m_BDOP_BestSlabIdx = i;
            }
        }
        else //\todo This should NEVER happen if all DCR.E contain >0 DCR.V, but for SOME REASON does...
        {
            GEO_ASSERT( ed.m_NumVertices == 0 ); //ENSURE that this happens ONLY when the DCR.E is empty...
            GEO_LOG_WARNING("Create_DCR_TetSolidShape3_From_TriSurfaceShape3() DCR.E[%u] has empty BDOP, setting to whole DCR.E", it_e );
            ed.m_BDOP[0] = Interval(0,1);
            ed.m_BDOP[1] = Interval(0,1);
            ed.m_BDOP[2] = Interval(0,1);
            ed.m_BDOP[3] = Interval(0,1);
            ed.m_BDOP_BestSlabIdx = 0;
        }
        // Compute per-patch aggregate stuff
        for( unsigned int it_p=ed.m_FirstPID; it_p<ed.m_FirstPID+ed.m_NumPatches; it_p++ )
        {
            DCR_TetSolidShape3::Patch& pd( pDCR->m_vecP[it_p] );
            // Vector area = \sum N*dA
            pd.m_VectorArea = Vec3::Zero();
            for( unsigned int it_t=pd.m_FirstTID; it_t<pd.m_FirstTID+pd.m_NumTriangles; it_t++ )
            {
                const DCR_TetSolidShape3::Triangle& tri( pDCR->m_vecT[it_t] );
                pd.m_VectorArea += mal::Cross( pDCR->m_vecV[tri.GetVID(1)] - pDCR->m_vecV[tri.GetVID(0)],
                                               pDCR->m_vecV[tri.GetVID(2)] - pDCR->m_vecV[tri.GetVID(0)] );
            }
            pd.m_VectorArea *= Real(0.5);
        }
    }

    // LOG baked elements
    // {
    //     for( unsigned int it_e=0; it_e<pDCR->m_NumElements; it_e++ )
    //     {
    //         const DCR_TetSolidShape3::Element& ed( pDCR->m_vecE[it_e] );
    //         GEO_LOG("Element[%u]: #P = %u [%u..%u), #V = %u [%u..%u)",//, #T = %u (%u..%u), ",
    //                 it_e, ed.m_NumPatches, ed.m_FirstPID, ed.m_FirstPID+ed.m_NumPatches,
    //                 ed.m_NumVertices, ed.m_FirstVID, ed.m_FirstVID + ed.m_NumVertices );
    //         for( unsigned int it_v=ed.m_FirstVID; it_v<ed.m_FirstVID+ed.m_NumVertices; it_v++ )
    //         {
    //             const Vec3& v( pDCR->m_vecV[it_v] );
    //             GEO_LOG("V[%u] = (%f,%f,%f)", it_v, v.x(), v.y(), v.z() );
    //         }
    //         for( unsigned int it_p=ed.m_FirstPID; it_p<ed.m_FirstPID+ed.m_NumPatches; it_p++ )
    //         {
    //             const DCR_TetSolidShape3::Patch& pd( pDCR->m_vecP[it_p] );
    //             GEO_LOG("Patch[%u]: #T = %u [%u..%u), ",
    //                     it_p, pd.m_NumTriangles, pd.m_FirstTID, pd.m_FirstTID+pd.m_NumTriangles );
    //             for( unsigned int it_t=pd.m_FirstTID; it_t<pd.m_FirstTID+pd.m_NumTriangles; it_t++ )
    //             {
    //                 const DCR_TetSolidShape3::Triangle& tri( pDCR->m_vecT[it_t] );
    //                 GEO_LOG("T[%u] = (%u,%u,%u)", it_t, tri.m_vecVID[0], tri.m_vecVID[1], tri.m_vecVID[2] );
    //             }
    //         }
    //     }
    // }

    // TEMPORAL: check to ensure valid ids...
    GEO_LOG_ERROR_IF( !Check_DCR_TetSolidShape3(*pDCR), "Check_DCR_TetSolidShape3() failed!" );

#ifdef __ENABLE_TETSS_DCR_PATCH_BDT
    // Bake per-patch BDT
    {
        std::vector<uint32> mapTID( pDCR->m_NumTriangles, 0xFFFFFFFF );
        uint32 processed_tid(0);
        for( unsigned int it_e=0; it_e<pDCR->m_NumElements; it_e++ )
        {
            const DCR_TetSolidShape3::Element& ed( pDCR->m_vecE[it_e] );
            // Get inverse barycentric transform

            // Gather element nodes
            uint32 vec_tet_vid[4] = { solid.T_VID(it_e,0), solid.T_VID(it_e,1), solid.T_VID(it_e,2), solid.T_VID(it_e,3) };
            Vec3 vec_tet_pos[4] = { solid.V_Pos_0( vec_tet_vid[0] ), solid.V_Pos_0( vec_tet_vid[1] ), solid.V_Pos_0( vec_tet_vid[2] ), solid.V_Pos_0( vec_tet_vid[3] ) };
            // Compute element barycentric coordinates matrix and its inverse
            Mat4x4 Bm( 1, 1, 1, 1,
                       vec_tet_pos[0].x(), vec_tet_pos[1].x(), vec_tet_pos[2].x(), vec_tet_pos[3].x(),
                       vec_tet_pos[0].y(), vec_tet_pos[1].y(), vec_tet_pos[2].y(), vec_tet_pos[3].y(),
                       vec_tet_pos[0].z(), vec_tet_pos[1].z(), vec_tet_pos[2].z(), vec_tet_pos[3].z() );
            Mat4x4 invBm( mal::Inverse( Bm ) ); //\todo invBm could be precomputed in DCR.ED!!
            // Compute per-patch BDT
            for( unsigned int it_pie=0; it_pie<ed.m_NumPatches; it_pie++ )
            {
                uint32 pid( ed.m_FirstPID + it_pie );
                // GEO_LOG("DCR.P[%u].BDT",pid);
                const DCR_TetSolidShape3::Patch& pd( pDCR->m_vecP[pid] );

                typedef std::pair<DCR_TetSolidShape3::Triangle,uint32> patch_segment_type;
                std::vector< patch_segment_type > vecPS;

                // Gather Patch-Triangle sub-array with their current TID
                for( unsigned int it_tip=0; it_tip<pd.m_NumTriangles; it_tip++ )
                {
                    uint32 tid( pd.m_FirstTID + it_tip );
                    vecPS.push_back( std::make_pair(pDCR->m_vecT[tid],tid) );
                    processed_tid++;
                }
                GEO_ASSERT(!vecPS.empty());

                // Generate BDT
                // GEO_LOG( "**** Creating DCR.P[%u].BDT with %u DCR.S", pid, pd.m_NumTriangles );
                std::vector<DCR_TetSolidShape3::Patch::bdt_node_tip_range> stackBDTN;
                stackBDTN.push_back( DCR_TetSolidShape3::Patch::bdt_node_tip_range(0,vecPS.size()) ); //root-range
                unsigned int num_bdtn(0);
                while( !stackBDTN.empty() )
                {
                    // Pop PS subarray
                    DCR_TetSolidShape3::Patch::bdt_node_tip_range bdtn_sr( stackBDTN.back() );
                    stackBDTN.pop_back();

                    // Compute BDOP for current PS subarray
                    bv::BDOP4 bdop;
                    for( unsigned int it_tin=bdtn_sr.first; it_tin<bdtn_sr.second; it_tin++ )
                    {
                        //Get the original TID for a PS in the sub-array and use it to retrieve its geometry
                        uint32 tid( vecPS[it_tin].second );
                        const DCR_TetSolidShape3::Triangle& td( pDCR->m_vecT[tid] );
                        Vec4 vec_b[3] = { invBm * mal::Concat(1,pDCR->m_vecV[td.GetVID(0)]),
                                          invBm * mal::Concat(1,pDCR->m_vecV[td.GetVID(1)]),
                                          invBm * mal::Concat(1,pDCR->m_vecV[td.GetVID(2)]) };
                        for( int i=0; i<3; i++ ) //vertices
                            for( int j=0; j<4; j++ ) //bdop-axis
                                bdop[j].Merge( mal::Clamp01(vec_b[i][j]) );
                    }
                    // Select longest axis \todo Consider better heuristics
                    unsigned int best_axis(0);
                    for( unsigned int i=1; i<4; i++ )
                        if( bdop[i].Length() > bdop[best_axis].Length() )
                            best_axis = i;

                    // Sort PS subarray along chosen axis
                    std::sort( vecPS.begin() + bdtn_sr.first, vecPS.begin() + bdtn_sr.second,
                               [pDCR, invBm, best_axis](const patch_segment_type& ps1, const patch_segment_type& ps2)
                               {
                                   Vec4 ps1_b( invBm * mal::Concat(1,
                                                                   Real(0.333333) * (pDCR->m_vecV[pDCR->m_vecT[ps1.second].GetVID(0)]
                                                                                     + pDCR->m_vecV[pDCR->m_vecT[ps1.second].GetVID(1)]
                                                                                     + pDCR->m_vecV[pDCR->m_vecT[ps1.second].GetVID(2)]) ) );
                                   Vec4 ps2_b( invBm * mal::Concat(1,
                                                                   Real(0.333333) * (pDCR->m_vecV[pDCR->m_vecT[ps2.second].GetVID(0)]
                                                                                     + pDCR->m_vecV[pDCR->m_vecT[ps2.second].GetVID(1)]
                                                                                     + pDCR->m_vecV[pDCR->m_vecT[ps2.second].GetVID(2)]) ) );
                                   return ps1_b[best_axis] < ps2_b[best_axis];
                                   // return ps1.second > ps2.second;
                               }
                        );

                    // Create BDTN for current PS subarray at the first PS
                    DCR_TetSolidShape3::Triangle& bdtn( vecPS[bdtn_sr.first].first );
                    bdtn.BDTN_Init( bdop );
                    // TEMP: DEBUG
                    // {
                    //     //TEMP: Log node
                    //     bv::BDOP4 bdop_r( bdtn.BDTN_BDOPq() );
                    //     DCR_TetSolidShape3::Triangle::BDOP4q bdop_q( bdtn.BDTN_BDOPq() );
                    //     GEO_LOG("BDTN[%u] best_axis %u, range [%u,%u]", num_bdtn, best_axis, bdtn_sr.first, bdtn_sr.second );
                    //     GEO_LOG("- BDOP  [%f,%f]x[%f,%f]x[%f,%f]x[%f,%f]", bdop[0].Min(), bdop[0].Max(), bdop[1].Min(), bdop[1].Max(), bdop[2].Min(), bdop[2].Max(), bdop[3].Min(), bdop[3].Max() );
                    //     GEO_LOG("- BDOPr [%f,%f]x[%f,%f]x[%f,%f]x[%f,%f]", bdop_r[0].Min(), bdop_r[0].Max(), bdop_r[1].Min(), bdop_r[1].Max(), bdop_r[2].Min(), bdop_r[2].Max(), bdop_r[3].Min(), bdop_r[3].Max() );
                    //     GEO_LOG("- BDOPq [%d,%d]x[%d,%d]x[%d,%d]x[%d,%d]", bdop_q[0].Min(), bdop_q[0].Max(), bdop_q[1].Min(), bdop_q[1].Max(), bdop_q[2].Min(), bdop_q[2].Max(), bdop_q[3].Min(), bdop_q[3].Max() );
                    //     //TEMP: Check conservative quantization
                    //     for( int i=0; i<4; i++ )
                    //     {
                    //         GEO_ASSERT( bdop_r[i].Min() <= bdop[i].Min() && bdop_r[i].Max() >= bdop[i].Max() );
                    //         GEO_LOG_WARNING_IF( !mal::ApproxEq( bdop_r[i].Min(), bdop[i].Min(), mal::Rcp<Real>(1<<DCR_TetSolidShape3::Triangle::eBDT_IntervalQuantizationBits) )
                    //                             || !mal::ApproxEq( bdop_r[i].Max(), bdop[i].Max(), mal::Rcp<Real>(1<<DCR_TetSolidShape3::Triangle::eBDT_IntervalQuantizationBits) ),
                    //                             "Low BDOPq precision" );
                    //     }
                    // }
                    num_bdtn++;

                    // Skip the BDTN PS and split the remaining PS subarray at the median
                    int length_sr( bdtn_sr.second - bdtn_sr.first );
                    int remaining_sr( length_sr - 1 );
                    if( remaining_sr > 1 )
                    {
                        // GEO_LOG("Stacking 2 L/R");
                        stackBDTN.push_back( DCR_TetSolidShape3::Patch::bdt_node_tip_range( bdtn_sr.first+1, bdtn_sr.first+1+remaining_sr/2 ) );
                        stackBDTN.push_back( DCR_TetSolidShape3::Patch::bdt_node_tip_range( bdtn_sr.first+1+remaining_sr/2, bdtn_sr.second ) );
                    }
                    else if( remaining_sr == 1 )
                    {
                        // GEO_LOG("Stacking 1 L");
                        stackBDTN.push_back( DCR_TetSolidShape3::Patch::bdt_node_tip_range( bdtn_sr.first+1, bdtn_sr.first+2 ) );
                    }
                }

                // Bake Element Triangle sub-array into global segment sub-array and fill mapTID
                for( unsigned int it_tip=0; it_tip<vecPS.size(); it_tip++ )
                {
                    uint32 old_tid( vecPS[it_tip].second );
                    uint32 new_tid( pd.m_FirstTID + it_tip );
                    pDCR->m_vecT[new_tid] = vecPS[it_tip].first;
                    mapTID[old_tid] = new_tid;
                    processed_tid++;
                }
            }
        }

        /*TEMP: Check correct remapping
          GEO_ASSERT( processed_tid == 2*pDCR->m_NumTriangles );
          uint32 num_unmapped_tid(0);
          for( unsigned int it_s=0; it_s<pDCR->m_NumTriangles; it_s++ )
          if( mapTID[it_s] == 0xFFFFFFFF )
          num_unmapped_tid++;
          GEO_LOG_ERROR_IF( num_unmapped_tid > 0, "NumUnmappedTID = %u", num_unmapped_tid );
          GEO_ASSERT( num_unmapped_tid == 0 );
        */

        // Remap TID in DCR.T.m_vecNTID
        for( unsigned int it_tri=0; it_tri<pDCR->m_NumTriangles; it_tri++ )
        {
            DCR_TetSolidShape3::Triangle& tri( pDCR->m_vecT[it_tri] );
            for( unsigned int it_ntit=0; it_ntit<3; it_ntit++ )
                if( tri.GetNTID(it_ntit) != DCR_TetSolidShape3::Triangle::cInvalidTID )
                    tri.SetNTID( it_ntit, mapTID[tri.GetNTID(it_ntit)] );
        }
        //TEMP: UNNECESSARY...
        // for( unsigned int it_p=0; it_p<pDCR->m_NumPatches; it_p++ )
        // {
        //     DCR_TetSolidShape3::Patch& patch( pDCR->m_vecP[it_p] );
        //     for( unsigned int it_npip=0; it_npip<patch.m_NumNeighbours; it_npip++ )
        //         p.m_vecNPID[it_npip] = mapPID[ patch.m_vecNPID[it_npip] ];
        // }

        //TEMP: DEBUG LOG
        // for( unsigned int it_p=0; it_p < pDCR->m_NumPatches; it_p++ )
        // {
        //     uint32 num_visited_nodes(0);
        //     const geo::DCR_TetSolidShape3::Patch& pd( pDCR->m_vecP[it_p] );
        //     typedef std::pair< DCR_TetSolidShape3::Patch::bdt_node_tip_range,uint32> stack_entry_type;
        //     std::vector< stack_entry_type > stackBDTN;
        //     stackBDTN.push_back( stack_entry_type(DCR_TetSolidShape3::Patch::bdt_node_tip_range(0,pd.m_NumTriangles),0) );
        //     GEO_LOG("DCR.P[%u].BDT",it_p);
        //     while( !stackBDTN.empty() )
        //     {
        //         num_visited_nodes++;
        //         // Pop PS subarray
        //         stack_entry_type se( stackBDTN.back() );
        //         stackBDTN.pop_back();
        //         DCR_TetSolidShape3::Patch::bdt_node_tip_range node_sr( se.first );
        //         // Get node
        //         const DCR_TetSolidShape3::Triangle& bdtn( pDCR->m_vecT[pd.m_FirstTID+node_sr.first] );
        //         bv::BDOP4 bdop( bdtn.BDTN_BDOPq() );
        //         GEO_LOG("BDTN level %u, range [%u,%u], volume %f", se.second, node_sr.first, node_sr.second, bdop[0].Length()*bdop[1].Length()*bdop[2].Length()*bdop[3].Length() );
        //         GEO_LOG("- BDOP [%f,%f]x[%f,%f]x[%f,%f]x[%f,%f]", bdop[0].Min(), bdop[0].Max(), bdop[1].Min(), bdop[1].Max(), bdop[2].Min(), bdop[2].Max(), bdop[3].Min(), bdop[3].Max() );
        //         // Recurse
        //         int length_sr( node_sr.second - node_sr.first );
        //         int remaining_sr( length_sr - 1 );
        //         if( remaining_sr > 1 )
        //         {
        //             stackBDTN.push_back( stack_entry_type( DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_sr.first+1, node_sr.first+1+remaining_sr/2 ), se.second+1 ) );
        //             stackBDTN.push_back( stack_entry_type( DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_sr.first+1+remaining_sr/2, node_sr.second ), se.second+1) );
        //         }
        //         else if( remaining_sr == 1 )
        //             stackBDTN.push_back( stack_entry_type( DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_sr.first+1, node_sr.first+2 ), se.second+1) );
        //     }
        //     GEO_LOG_ERROR_IF( num_visited_nodes != pd.m_NumTriangles, "Recursive visit skips some nodes %u < %u!!", num_visited_nodes, pd.m_NumTriangles );
        // }
        //
    }
#endif //__ENABLE_MSS_DCR_PATCH_BDT

    // TEMPORAL: Re-check to ensure that BDOP-Tree has not broken anything
    GEO_LOG_ERROR_IF( !Check_DCR_TetSolidShape3(*pDCR), "Check_DCR_TetSolidShape3() failed!" );

    return pDCR;
}

bool Check_DCR_TetSolidShape3( const DCR_TetSolidShape3& dcr )
{
    bool bOk(true);
    unsigned int num_open_edges(0);
    for( unsigned int it_e=0; it_e<dcr.m_NumElements; it_e++ )
    {
        const DCR_TetSolidShape3::Element& ed( dcr.m_vecE[it_e] );

        // GEO_LOG("Element[%u]: #P = %u [%u..%u), #V = %u [%u..%u)",//, #T = %u (%u..%u), ",
        //         it_e, ed.m_NumPatches, ed.m_FirstPID, ed.m_FirstPID+ed.m_NumPatches,
        //         ed.m_NumVertices, ed.m_FirstVID, ed.m_FirstVID + ed.m_NumVertices );

        // for( unsigned int it_v=ed.m_FirstVID; it_v<ed.m_FirstVID+ed.m_NumVertices; it_v++ )
        // {
        //     const Vec3& v( dcr.m_vecV[it_v] );
        //     // GEO_LOG("V[%u] = (%f,%f,%f)", it_v, v.x(), v.y(), v.z() );
        // }
        for( unsigned int it_p=ed.m_FirstPID; it_p<ed.m_FirstPID+ed.m_NumPatches; it_p++ )
        {
            const DCR_TetSolidShape3::Patch& pd( dcr.m_vecP[it_p] );
            /*\todo ENABLE when NPID become available
            if( pd.GetNPID(0) < 0 || pd.GetNPID(0) >= dcr.m_NumPatches
                || pd.GetNPID(1) < 0 || pd.GetNPID(1) >= dcr.m_NumPatches
                || pd.GetNPID(2) < 0 || pd.GetNPID(2) >= dcr.m_NumPatches )
                GEO_LOG_ERROR( "DCR.E[%u].P[%u] = DCR.P[%u] has NPID (%u,%u,%u,%u) outside [%u...%u)",
                               it_e, it_p-ed.m_FirstPID, it_p,
                               pd.GetNPID(0), pd.GetNPID(1), pd.GetNPID(2), pd.GetNPID(3),
                               0, dcr.m_NumPatches );
            */
            for( unsigned int it_t=pd.m_FirstTID; it_t<pd.m_FirstTID+pd.m_NumTriangles; it_t++ )
            {
                const DCR_TetSolidShape3::Triangle& tri( dcr.m_vecT[it_t] );
                if( tri.GetVID(0) < ed.m_FirstVID || tri.GetVID(0) >= ed.m_FirstVID+ed.m_NumVertices
                    || tri.GetVID(1) < ed.m_FirstVID || tri.GetVID(1) >= ed.m_FirstVID+ed.m_NumVertices
                    || tri.GetVID(2) < ed.m_FirstVID || tri.GetVID(2) >= ed.m_FirstVID+ed.m_NumVertices )
                {
                    GEO_LOG_ERROR( "DCR.E[%u].P[%u].T[%u] = DCR.T[%u] has VID (%u,%u,%u) outside [%u...%u) #DCR.P.T = %u",
                                   it_e, it_p-ed.m_FirstPID, it_t-pd.m_FirstTID, it_t,
                                   tri.GetVID(0), tri.GetVID(1), tri.GetVID(2),
                                   ed.m_FirstVID, ed.m_FirstVID+ed.m_NumVertices, pd.m_NumTriangles );
                    bOk = false;
                }
                if( tri.GetNTID(0) < 0 || (tri.GetNTID(0) != DCR_TetSolidShape3::Triangle::cInvalidTID && tri.GetNTID(0) >= dcr.m_NumTriangles)
                    || tri.GetNTID(1) < 0 || (tri.GetNTID(1) != DCR_TetSolidShape3::Triangle::cInvalidTID && tri.GetNTID(1) >= dcr.m_NumTriangles)
                    || tri.GetNTID(2) < 0 || (tri.GetNTID(2) != DCR_TetSolidShape3::Triangle::cInvalidTID && tri.GetNTID(2) >= dcr.m_NumTriangles) )
                {
                    GEO_LOG_ERROR( "DCR.E[%u].P[%u].T[%u] = DCR.T[%u] has NTID (%u,%u,%u) outside [%u...%u)",
                                   it_e, it_p-ed.m_FirstPID, it_t-pd.m_FirstTID, it_t,
                                   tri.GetNTID(0), tri.GetNTID(1), tri.GetNTID(2),
                                   0, dcr.m_NumTriangles );
                    bOk = false;
                }
                for( unsigned int it_ntit=0; it_ntit<3; it_ntit++ )
                    if( tri.GetNTID(it_ntit) == DCR_TetSolidShape3::Triangle::cInvalidTID )
                        num_open_edges++;
            }
        }
    }
    GEO_LOG_ERROR_IF( num_open_edges > 0, "DCR has %u open edges", num_open_edges );
    return bOk && num_open_edges == 0;;
}


#ifdef __GEO_TETSS_ENABLE_BVH

/* Specialization of BV(BDOP) with tight AABB(BDOP) computation based
   on the papers:
   - [2007] Adaptive Deformations with Fast Tight Bounds
   - [2008] Tight and efficient surface bounds in meshless animation
*/
template <>
void GEBV_TetSolidShape3_BDOP( const TetSolidShape3* p_solid, const DCR_TetSolidShape3* p_dcr, const Transform3& tr, const TetSolidShape3::sdof_type* vec_node_pos,
                               tetsolid3_feature_index_type eid, bv::AABB3& bv )
{
    GEO_ASSERT( 0 != p_dcr && eid >= 0 && eid < p_dcr->m_NumElements );
    const DCR_TetSolidShape3::Element& ed( p_dcr->m_vecE[eid] );
    // compute global node pos
    Vec3 node_pos[4] = { tr*vec_node_pos[ p_solid->T_VID(eid,0) ],
                         tr*vec_node_pos[ p_solid->T_VID(eid,1) ],
                         tr*vec_node_pos[ p_solid->T_VID(eid,2) ],
                         tr*vec_node_pos[ p_solid->T_VID(eid,3) ] };
    // sort node indices along X,Y
    int snie_x[4] = { 0, 1, 2, 3 };
    int snie_y[4] = { 0, 1, 2, 3 };
    int snie_z[4] = { 0, 1, 2, 3 };
    //We use Sort4_Bubbble, 2x faster than generic std::sort() and simpler/shorter code than Sort4_HC while slightly (5%) slower
    for( int i=0; i<3; i++ )
        for( int j=0; j<3; j++ )
        {
            if( node_pos[snie_x[j]].x() > node_pos[snie_x[j+1]].x() ) std::swap(snie_x[j],snie_x[j+1]);
            if( node_pos[snie_y[j]].y() > node_pos[snie_y[j+1]].y() ) std::swap(snie_y[j],snie_y[j+1]);
            if( node_pos[snie_z[j]].z() > node_pos[snie_z[j+1]].z() ) std::swap(snie_z[j],snie_z[j+1]);
        }
    /* STL generic sort was 2x expensive than Sort4_Bubbble
       std::sort( snie_x, snie_x+4, [node_pos](int a, int b){ return node_pos[a].x() < node_pos[b].x(); } );
       std::sort( snie_y, snie_y+4, [node_pos](int a, int b){ return node_pos[a].y() < node_pos[b].y(); } );
       std::sort( snie_z, snie_z+4, [node_pos](int a, int b){ return node_pos[a].z() < node_pos[b].z(); } );
    */
    // compute min/max \sa Eq(8) in [2007] Adaptive Deformations with Fast Tight Bounds \todo Consider manual optimization
    bv.SetMinMax( Vec3( ed.m_BDOP[snie_x[0]].Max()*node_pos[snie_x[0]].x() + ed.m_BDOP[snie_x[3]].Min()*node_pos[snie_x[3]].x() + ed.m_BDOP[snie_x[2]].Min()*node_pos[snie_x[2]].x()
                        + (1-ed.m_BDOP[snie_x[0]].Max()-ed.m_BDOP[snie_x[3]].Min()-ed.m_BDOP[snie_x[2]].Min()) * node_pos[snie_x[1]].x(),
                        ed.m_BDOP[snie_y[0]].Max()*node_pos[snie_y[0]].y() + ed.m_BDOP[snie_y[3]].Min()*node_pos[snie_y[3]].y() + ed.m_BDOP[snie_y[2]].Min()*node_pos[snie_y[2]].y()
                        + (1-ed.m_BDOP[snie_y[0]].Max()-ed.m_BDOP[snie_y[3]].Min()-ed.m_BDOP[snie_y[2]].Min()) * node_pos[snie_y[1]].y(),
                        ed.m_BDOP[snie_z[0]].Max()*node_pos[snie_z[0]].z() + ed.m_BDOP[snie_z[3]].Min()*node_pos[snie_z[3]].z() + ed.m_BDOP[snie_z[2]].Min()*node_pos[snie_z[2]].z()
                        + (1-ed.m_BDOP[snie_z[0]].Max()-ed.m_BDOP[snie_z[3]].Min()-ed.m_BDOP[snie_z[2]].Min()) * node_pos[snie_z[1]].z() ),
                  Vec3( ed.m_BDOP[snie_x[3]].Max()*node_pos[snie_x[3]].x() + ed.m_BDOP[snie_x[0]].Min()*node_pos[snie_x[0]].x() + ed.m_BDOP[snie_x[1]].Min()*node_pos[snie_x[1]].x()
                        + (1-ed.m_BDOP[snie_x[3]].Max()-ed.m_BDOP[snie_x[0]].Min()-ed.m_BDOP[snie_x[1]].Min()) * node_pos[snie_x[2]].x(),
                        ed.m_BDOP[snie_y[3]].Max()*node_pos[snie_y[3]].y() + ed.m_BDOP[snie_y[0]].Min()*node_pos[snie_y[0]].y() + ed.m_BDOP[snie_y[1]].Min()*node_pos[snie_y[1]].y()
                        + (1-ed.m_BDOP[snie_y[3]].Max()-ed.m_BDOP[snie_y[0]].Min()-ed.m_BDOP[snie_y[1]].Min()) * node_pos[snie_y[2]].y(),
                        ed.m_BDOP[snie_z[3]].Max()*node_pos[snie_z[3]].z() + ed.m_BDOP[snie_z[0]].Min()*node_pos[snie_z[0]].z() + ed.m_BDOP[snie_z[1]].Min()*node_pos[snie_z[1]].z()
                        + (1-ed.m_BDOP[snie_z[3]].Max()-ed.m_BDOP[snie_z[0]].Min()-ed.m_BDOP[snie_z[1]].Min()) * node_pos[snie_z[2]].z() ) );
}

/* Specialization of BV(BDOP) with tight DOP3K6(BDOP) computation, same as AABB(BDOP) */
template <>
void GEBV_TetSolidShape3_BDOP( const TetSolidShape3* p_solid, const DCR_TetSolidShape3* p_dcr, const Transform3& tr, const TetSolidShape3::sdof_type* vec_node_pos,
                               tetsolid3_feature_index_type eid, bv::DOP3_K6& bv )
{
    GEO_ASSERT( 0 != p_dcr && eid >= 0 && eid < p_dcr->m_NumElements );
    const DCR_TetSolidShape3::Element& ed( p_dcr->m_vecE[eid] );
    // compute global node pos
    Vec3 node_pos[4] = { tr*vec_node_pos[ p_solid->T_VID(eid,0) ],
                         tr*vec_node_pos[ p_solid->T_VID(eid,1) ],
                         tr*vec_node_pos[ p_solid->T_VID(eid,2) ],
                         tr*vec_node_pos[ p_solid->T_VID(eid,3) ] };
    // sort node indices along X,Y
    int snie_x[4] = { 0, 1, 2, 3 };
    int snie_y[4] = { 0, 1, 2, 3 };
    int snie_z[4] = { 0, 1, 2, 3 };
    for( int i=0; i<3; i++ )
        for( int j=0; j<3; j++ )
        {
            if( node_pos[snie_x[j]].x() > node_pos[snie_x[j+1]].x() ) std::swap(snie_x[j],snie_x[j+1]);
            if( node_pos[snie_y[j]].y() > node_pos[snie_y[j+1]].y() ) std::swap(snie_y[j],snie_y[j+1]);
            if( node_pos[snie_z[j]].z() > node_pos[snie_z[j+1]].z() ) std::swap(snie_z[j],snie_z[j+1]);

        }
    // compute min/max \sa Eq(8) in [2007] Adaptive Deformations with Fast Tight Bounds \todo Consider manual optimization
    bv.SetInterval<0>(
        Interval( ed.m_BDOP[snie_x[0]].Max()*node_pos[snie_x[0]].x() + ed.m_BDOP[snie_x[3]].Min()*node_pos[snie_x[3]].x() + ed.m_BDOP[snie_x[2]].Min()*node_pos[snie_x[2]].x()
                  + (1-ed.m_BDOP[snie_x[0]].Max()-ed.m_BDOP[snie_x[3]].Min()-ed.m_BDOP[snie_x[2]].Min()) * node_pos[snie_x[1]].x(),
                  ed.m_BDOP[snie_x[3]].Max()*node_pos[snie_x[3]].x() + ed.m_BDOP[snie_x[0]].Min()*node_pos[snie_x[0]].x() + ed.m_BDOP[snie_x[1]].Min()*node_pos[snie_x[1]].x()
                  + (1-ed.m_BDOP[snie_x[3]].Max()-ed.m_BDOP[snie_x[0]].Min()-ed.m_BDOP[snie_x[1]].Min()) * node_pos[snie_x[2]].x() ) );
    bv.SetInterval<1>(
        Interval( ed.m_BDOP[snie_y[0]].Max()*node_pos[snie_y[0]].y() + ed.m_BDOP[snie_y[3]].Min()*node_pos[snie_y[3]].y() + ed.m_BDOP[snie_y[2]].Min()*node_pos[snie_y[2]].y()
                  + (1-ed.m_BDOP[snie_y[0]].Max()-ed.m_BDOP[snie_y[3]].Min()-ed.m_BDOP[snie_y[2]].Min()) * node_pos[snie_y[1]].y(),
                  ed.m_BDOP[snie_y[3]].Max()*node_pos[snie_y[3]].y() + ed.m_BDOP[snie_y[0]].Min()*node_pos[snie_y[0]].y() + ed.m_BDOP[snie_y[1]].Min()*node_pos[snie_y[1]].y()
                  + (1-ed.m_BDOP[snie_y[3]].Max()-ed.m_BDOP[snie_y[0]].Min()-ed.m_BDOP[snie_y[1]].Min()) * node_pos[snie_y[2]].y() ) );
    bv.SetInterval<2>(
        Interval( ed.m_BDOP[snie_z[0]].Max()*node_pos[snie_z[0]].z() + ed.m_BDOP[snie_z[3]].Min()*node_pos[snie_z[3]].z() + ed.m_BDOP[snie_z[2]].Min()*node_pos[snie_z[2]].z()
                  + (1-ed.m_BDOP[snie_z[0]].Max()-ed.m_BDOP[snie_z[3]].Min()-ed.m_BDOP[snie_z[2]].Min()) * node_pos[snie_z[1]].z(),
                  ed.m_BDOP[snie_z[3]].Max()*node_pos[snie_z[3]].z() + ed.m_BDOP[snie_z[0]].Min()*node_pos[snie_z[0]].z() + ed.m_BDOP[snie_z[1]].Min()*node_pos[snie_z[1]].z()
                  + (1-ed.m_BDOP[snie_z[3]].Max()-ed.m_BDOP[snie_z[0]].Min()-ed.m_BDOP[snie_z[1]].Min()) * node_pos[snie_z[2]].z() ) );
}

/* Specialization of BV(BDOP) with tight AABB(BDOP) computation based
   on the papers:
   - [2007] Adaptive Deformations with Fast Tight Bounds
   - [2008] Tight and efficient surface bounds in meshless animation
*/
template <>
void GEBV_TetSolidShape3_BDOP( const TetSolidShape3* p_solid, const DCR_TetSolidShape3* p_dcr, const Transform3& tr, const TetSolidShape3::sdof_type* vec_node_pos,
                               tetsolid3_feature_index_type eid, bv::DOP3_K14& bv )
{
    GEO_ASSERT( 0 != p_dcr && eid >= 0 && eid < p_dcr->m_NumElements );
    const DCR_TetSolidShape3::Element& ed( p_dcr->m_vecE[eid] );
    // compute global node pos
    Vec3 node_pos[4] = { tr*vec_node_pos[ p_solid->T_VID(eid,0) ],
                         tr*vec_node_pos[ p_solid->T_VID(eid,1) ],
                         tr*vec_node_pos[ p_solid->T_VID(eid,2) ],
                         tr*vec_node_pos[ p_solid->T_VID(eid,3) ] };
    //\todo FIRST 3 axis are shared with DOP3_K6, try to reuse code!!
    // sort node indices along X,Y
    int snie_x[4] = { 0, 1, 2, 3 };
    int snie_y[4] = { 0, 1, 2, 3 };
    int snie_z[4] = { 0, 1, 2, 3 };
    for( int i=0; i<3; i++ )
        for( int j=0; j<3; j++ )
        {
            if( node_pos[snie_x[j]].x() > node_pos[snie_x[j+1]].x() ) std::swap(snie_x[j],snie_x[j+1]);
            if( node_pos[snie_y[j]].y() > node_pos[snie_y[j+1]].y() ) std::swap(snie_y[j],snie_y[j+1]);
            if( node_pos[snie_z[j]].z() > node_pos[snie_z[j+1]].z() ) std::swap(snie_z[j],snie_z[j+1]);
        }
    // compute min/max \sa Eq(8) in [2007] Adaptive Deformations with Fast Tight Bounds \todo Consider manual optimization
    bv.SetInterval<0>(
        Interval( ed.m_BDOP[snie_x[0]].Max()*node_pos[snie_x[0]].x() + ed.m_BDOP[snie_x[3]].Min()*node_pos[snie_x[3]].x() + ed.m_BDOP[snie_x[2]].Min()*node_pos[snie_x[2]].x()
                  + (1-ed.m_BDOP[snie_x[0]].Max()-ed.m_BDOP[snie_x[3]].Min()-ed.m_BDOP[snie_x[2]].Min()) * node_pos[snie_x[1]].x(),
                  ed.m_BDOP[snie_x[3]].Max()*node_pos[snie_x[3]].x() + ed.m_BDOP[snie_x[0]].Min()*node_pos[snie_x[0]].x() + ed.m_BDOP[snie_x[1]].Min()*node_pos[snie_x[1]].x()
                  + (1-ed.m_BDOP[snie_x[3]].Max()-ed.m_BDOP[snie_x[0]].Min()-ed.m_BDOP[snie_x[1]].Min()) * node_pos[snie_x[2]].x() ) );
    bv.SetInterval<1>(
        Interval( ed.m_BDOP[snie_y[0]].Max()*node_pos[snie_y[0]].y() + ed.m_BDOP[snie_y[3]].Min()*node_pos[snie_y[3]].y() + ed.m_BDOP[snie_y[2]].Min()*node_pos[snie_y[2]].y()
                  + (1-ed.m_BDOP[snie_y[0]].Max()-ed.m_BDOP[snie_y[3]].Min()-ed.m_BDOP[snie_y[2]].Min()) * node_pos[snie_y[1]].y(),
                  ed.m_BDOP[snie_y[3]].Max()*node_pos[snie_y[3]].y() + ed.m_BDOP[snie_y[0]].Min()*node_pos[snie_y[0]].y() + ed.m_BDOP[snie_y[1]].Min()*node_pos[snie_y[1]].y()
                  + (1-ed.m_BDOP[snie_y[3]].Max()-ed.m_BDOP[snie_y[0]].Min()-ed.m_BDOP[snie_y[1]].Min()) * node_pos[snie_y[2]].y() ) );
    bv.SetInterval<2>(
        Interval( ed.m_BDOP[snie_z[0]].Max()*node_pos[snie_z[0]].z() + ed.m_BDOP[snie_z[3]].Min()*node_pos[snie_z[3]].z() + ed.m_BDOP[snie_z[2]].Min()*node_pos[snie_z[2]].z()
                  + (1-ed.m_BDOP[snie_z[0]].Max()-ed.m_BDOP[snie_z[3]].Min()-ed.m_BDOP[snie_z[2]].Min()) * node_pos[snie_z[1]].z(),
                  ed.m_BDOP[snie_z[3]].Max()*node_pos[snie_z[3]].z() + ed.m_BDOP[snie_z[0]].Min()*node_pos[snie_z[0]].z() + ed.m_BDOP[snie_z[1]].Min()*node_pos[snie_z[1]].z()
                  + (1-ed.m_BDOP[snie_z[3]].Max()-ed.m_BDOP[snie_z[0]].Min()-ed.m_BDOP[snie_z[1]].Min()) * node_pos[snie_z[2]].z() ) );
    // Extra axis 3..6
    for( int a=3; a<7; a++ )
    {
        int snie_a[4] = { 0, 1, 2, 3 };
        Real vec_proj_a[4] = { bv::GKDOP_Dot<3,14>(node_pos[0],a), bv::GKDOP_Dot<3,14>(node_pos[1],a), bv::GKDOP_Dot<3,14>(node_pos[2],a), bv::GKDOP_Dot<3,14>(node_pos[3],a) };
        for( int i=0; i<3; i++ ) for( int j=0; j<3; j++ ) if( vec_proj_a[snie_a[j]] > vec_proj_a[snie_a[j+1]] ) std::swap(snie_a[j],snie_a[j+1]);
        bv.SetInterval( a,
                        Interval( ed.m_BDOP[snie_a[0]].Max()*vec_proj_a[snie_a[0]] + ed.m_BDOP[snie_a[3]].Min()*vec_proj_a[snie_a[3]] + ed.m_BDOP[snie_a[2]].Min()*vec_proj_a[snie_a[2]]
                                  + (1-ed.m_BDOP[snie_a[0]].Max()-ed.m_BDOP[snie_a[3]].Min()-ed.m_BDOP[snie_a[2]].Min()) * vec_proj_a[snie_a[1]],
                                  ed.m_BDOP[snie_a[3]].Max()*vec_proj_a[snie_a[3]] + ed.m_BDOP[snie_a[0]].Min()*vec_proj_a[snie_a[0]] + ed.m_BDOP[snie_a[1]].Min()*vec_proj_a[snie_a[1]]
                                  + (1-ed.m_BDOP[snie_a[3]].Max()-ed.m_BDOP[snie_a[0]].Min()-ed.m_BDOP[snie_a[1]].Min()) * vec_proj_a[snie_a[2]] ) );
    }
}

/* Specialization of BV(BDOP) for Sphere(BDOP)
   \todo Fallback to BV(E)
*/
template <>
void GEBV_TetSolidShape3_BDOP( const TetSolidShape3* p_solid, const DCR_TetSolidShape3* p_dcr, const Transform3& tr, const TetSolidShape3::sdof_type* vec_node_pos,
                               tetsolid3_feature_index_type eid, bv::Sphere3& bv )
{
    GEBV_TetSolidShape3_E<tetsolid3_feature_index_type,bv::Sphere3>( p_solid, tr, vec_node_pos, eid, bv );
}

#endif //__GEO_TETSS_ENABLE_BVH

} //namespace geo
