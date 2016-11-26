#ifndef GEO_SHAPE_TRISURFACESHAPE3_H
#define GEO_SHAPE_TRISURFACESHAPE3_H

#include <Geo/shape/IShape.h>
#include <Geo/bv/bv.h> // BV types for ComputeBV()
#include <vector>

#define __GEO_TRISS_ENABLE_BVH
#ifdef __GEO_TRISS_ENABLE_BVH
#  include <Geo/bv/GBoundingVolumeHierarchy.h> //TEMPORAL: quick and dirty BVH
#endif

namespace geo
{

// We use specific TSS V/T counts (uint32) instead of
// geo::feature_index_type (uint16) to allow for LARGE triangle meshes
// (for viz)
typedef uint32 trisurface3_feature_index_type;
typedef uint32 trisurface3_feature_count_type;
static const uint32 cTriSurface3_InvalidFeatureIndex = 0xFFFFFFFF;

// Req by iterator_around_vertex_ccw \todo This maybe NOT BE REQUIRED AT ALL by now... do it when used
struct triangle_vertex_topology
{
    trisurface3_feature_index_type m_TID; //Smallest incident TID
};
struct triangle_topology
{
    trisurface3_feature_index_type m_vecVID[3];
    trisurface3_feature_index_type m_vecNTID[3];
};

typedef std::pair<trisurface3_feature_index_type,int> trisurface3_edge_id_type;

#ifdef __GEO_TRISS_ENABLE_BVH
//TEMPORAL: Consider generic IBVH so that each MSS can have a user-specified or best fitting one
typedef bv::GBoundingVolumeHierarchy_ST_DG<bv::AABB3,trisurface3_feature_index_type> BVH_TriSurfaceShape3;
#endif

} //namespace geo

namespace geo
{

/* Triangle Solid 3D
   Stores:
   - Vertices
   - Triangles
*/
class TriSurfaceShape3: public IShape3
{
public:
    enum EConstants
    {
        cVersion = 2 //Added m_Flags
        //cVersion = 1 //Initial version
    };
    enum EFlags
    {
        eFlag_None      = 0,
        eFlag_Closed    = 1<<0,
        eFlag_Connected = 1<<1,
        //eFlag_Convex    = 1<<2 \todo
    };

public:
    //! \name Mandatory generic interface for IShape
    //@{
    static const unsigned int cDimension = 3;
    typedef mal::GTransform<Real,3> transform_type;
    typedef mal::GVec<Real,3> sdof_type;
    typedef mal::GVec<Real,3> vec_type;
    //@}

public:
    struct vertex_type: public triangle_vertex_topology {};
    struct triangle_type: public triangle_topology {};

public:
    TriSurfaceShape3();
    ~TriSurfaceShape3();

    //!\name Dynamic IShape implementation
    //@{
    EShapeType GetType() const { return eShape_TriSurface3; }
    unsigned int GetNumSDOF() const { return m_NumV; }
    const sdof_type* GetVecDefaultSDOF() const { return m_vecPoints; }
    const Real* GetVecDefaultDOF() const { return reinterpret_cast<const Real*>(m_vecPoints); }
    void ComputeBVD( bv::BoundingVolume3& bv, const transform_type& transform, const sdof_type* vec_sdof ) const;

    IDomainSampler3* CreateDomainSampler() const { return 0; } //\todo
    //@}

    //! Init from external arrays
    void SetBakedData( bool b_shared,
                       Flags32 flags,
                       uint32 num_vertices, uint32 num_triangles,
                       const Vec3* vec_points, const vertex_type* vec_vertices, const triangle_type* vec_triangles );
#ifdef __GEO_TRISS_ENABLE_BVH
    inline void SetBakedBVH_StrictlyNonshared_UglyHack( BVH_TriSurfaceShape3* p_nonshared_bvh ) { m_pBVH = p_nonshared_bvh; } //never shared, breaks SetBakedData() strict const enforcement
#endif

    inline Flags32 GetFlags() const { return m_Flags; }
    inline uint32 GetNumV() const { return m_NumV; }
    inline uint32 GetNumT() const { return m_NumT; }

    inline const Vec3* GetVecPoints() const { return m_vecPoints; }
    inline const vertex_type* GetVecV() const { return m_vecV; }
    inline const triangle_type* GetVecT() const { return m_vecT; }

    // \todo
    inline bool IsClosed() const { return m_Flags.Test(eFlag_Closed); }
    //inline bool IsSimplyConnected() const { return m_Flags.Test(eFlag_SimplyConnected); } //!< Check! false if there are holes
    inline bool IsConnected() const { return m_Flags.Test(eFlag_Connected); }
    inline bool IsConvex() const { return false; } //\todo

    // Vertex access
    inline Vec3 V_Pos_0( uint32 vid ) const { return m_vecPoints[vid]; }
    inline Vec3 V_Pos( uint32 vid, const sdof_type* vec_sdof ) const { return vec_sdof[vid]; }

    // Triangle access
    inline uint32 T_VID( uint32 tid, uint32 vit ) const { return m_vecT[tid].m_vecVID[vit]; }
    inline uint32 T_NTID( uint32 tid, uint32 ntit ) const { return m_vecT[tid].m_vecNTID[ntit]; }
    Vec3 T_Barycenter_0( uint32 tid ) const;
    Vec3 T_Barycenter( uint32 tid, const sdof_type* vec_sdof ) const;
    void T_Edge_0( uint32 tid, uint32 eit, Vec3& p0, Vec3& p1 ) const;

    //\todo iterator_around_vertex_ccw

    /*\todo Consider GenerateNormals, GenerateTangentSpace,
      GenerateUV... do NOT store them as members. Better separate
      these methods in a TriSurfaceOps namespace, not as member
      funcs */

#ifdef __GEO_TRISS_ENABLE_BVH
    inline BVH_TriSurfaceShape3* GetBVH() const { return m_pBVH; } //\note non-const until topology/geometry are split, because we need to refit it!
#endif

protected:
    void ClearBakedData();

protected:

    Flags32 m_Flags;

    uint32 m_NumV;
    uint32 m_NumT;

    const Vec3* m_vecPoints;
    const vertex_type* m_vecV;
    const triangle_type* m_vecT;

    uint32* m_pBuffer; //!< If exists, the shape is NOT SHARED

#ifdef __GEO_TRISS_ENABLE_BVH
    BVH_TriSurfaceShape3* m_pBVH; //strictly nonshared, by now...
#endif

    /*\todo
private:
    class DomainSampler;
    */
};

/* Editable Triangle Solid Shape 3D
*/
class EditableTriSurfaceShape3: public TriSurfaceShape3
{
public:
    EditableTriSurfaceShape3();
    ~EditableTriSurfaceShape3();

    void Clear();
    void Set( const TriSurfaceShape3& tss3 );

    void BeginEdition(); //\note Incremental, preserves BakedData if any, call Clear() to reset completely
    bool EndEdition();
    inline bool IsBeingEdited() const { return m_IsBeingEdited; }

    trisurface3_feature_index_type AddVertex( const Vec3& point );
    trisurface3_feature_index_type AddTriangle( trisurface3_feature_index_type vid0, trisurface3_feature_index_type vid1, trisurface3_feature_index_type vid2 );

    void Transform( const Transform3& tr );
    void Scale( Real scale );

    bool TryToCloseHoles( Real epsilon_length = mal::Infinity<Real>() ); //\todo EXPERIMENTAL, mostly for escalunya/Clip_TriSurfaceShape3_TetSolidShape3()
    unsigned int Simplify_EdgeCollapse( Real epsilon_length ); //\todo EXPERIMENTAL, mostly for escalunya/Clip_TriSurfaceShape3_TetSolidShape3()

//    inline bool IsClosed_EDITABLE() const { return m_IsClosed; }

#ifdef __GEO_TRISS_ENABLE_BVH
    bool AddBVH(); //bvh type!?
#endif

private:

    void ClearEditData();
    void RebuildBakedData();
    void RebuildAABB();
    bool FixDegeneracies(); //\note Return true if degeneracied found
    bool RebuildTopology(); //\note Return true if degeneracies found and new iter required
    void RebuildClosed();

    bool FixDegenerateVertices();
    bool FixDegenerateTriangles();
    bool RemoveIsolatedTriangles();

private:
    struct editable_vertex_type: public vertex_type
    {
        inline editable_vertex_type() : m_Pos(0,0,0) {} //\note Required by m_addV.resize()
        inline editable_vertex_type( const Vec3& pos )
        : m_Pos(pos) {}
        Vec3 m_Pos;
    };
    struct editable_triangle_type: public triangle_type
    {
        inline editable_triangle_type() {} //\note Required by m_addT.resize()
        inline editable_triangle_type( trisurface3_feature_index_type vid0, trisurface3_feature_index_type vid1, trisurface3_feature_index_type vid2 )
            {
                m_vecVID[0] = vid0; m_vecVID[1] = vid1; m_vecVID[2] = vid2;
                ResetNeigbours();
            }
        inline void ResetNeigbours()
            {
                m_vecNTID[0] = cTriSurface3_InvalidFeatureIndex;
                m_vecNTID[1] = cTriSurface3_InvalidFeatureIndex;
                m_vecNTID[2] = cTriSurface3_InvalidFeatureIndex;
            }
    };

    bool m_IsBeingEdited;
    std::vector<editable_vertex_type> m_addV;
    std::vector<editable_triangle_type> m_addT;

    bv::AABB3 m_AABB;
    bool m_IsClosed;
    bool m_IsConnected;
};


//\name Free functions
//@{
Real ComputeVolume( const TriSurfaceShape3& tss3, const Transform3& transform, const Vec3* vec_sdof );
//@}

//!\name Make methods
//@{
void Make_TriSurfaceShape3_Box( EditableTriSurfaceShape3& emss, const Vec3& half_sizes );
//\todo void Make_TriSurfaceShape3_Box( EditableTriSurfaceShape3& emss, const Vec3& half_sizes, int num_x, int num_y, int num_z );
//@}

#ifdef __GEO_TRISS_ENABLE_BVH
/*! BV(Element)
  \note generic implementation for any supported BV
  \todo Add pre-transformed vec_node_pos version
  IMPORTANT: Both std::bind and boost::bind REQUIRE a copy-constructor
             for const TriSurfaceShape3& surface param, which causes
             *LOTS* of trouble, therefore we pass a const
             TriSurfaceShape3* p_shape to avoid them.
*/
template <typename EIT, typename BVT>
inline void GEBV_TriSurfaceShape3_E( const TriSurfaceShape3* p_surface, const TriSurfaceShape3::transform_type& tr, const TriSurfaceShape3::sdof_type* vec_node_pos,
                                     EIT eid, BVT& bv )
{
    GEO_ASSERT( 0 != vec_node_pos );
    bv = BVT( tr * p_surface->V_Pos( p_surface->T_VID(eid,0), vec_node_pos ) );
    bv.Merge( tr * p_surface->V_Pos( p_surface->T_VID(eid,1), vec_node_pos ) );
    bv.Merge( tr * p_surface->V_Pos( p_surface->T_VID(eid,2), vec_node_pos ) );
};

#  ifdef __GEO_TRISS_ENABLE_BVH_TODO //\todo PORT MSS2 stuff to TSS3
/*! BV(BDOP)
  \note Has AABB2 and DOP2_K8 specializations, and falls back to BV(E) for Sphere2
  \todo Add pre-transformed vec_node_pos version
  \pre 0 <= eid << p_dcr->m_NumElements
*/
template <typename EIT, typename BVT>
void GEBV_MeshSolidShape2_BDOP( const MeshSolidShape2* p_mesh, const DCR_MeshSolidShape2* p_dcr, const Transform2& tr, const MeshSolidShape2::sdof_type* vec_node_pos,
                                EIT eid, BVT& bv );

/*! BV(BSlab)
  \note Has AABB2 and DOP2_K8 specializations, and falls back to BV(E) for Sphere2
  \todo Add pre-transformed vec_node_pos version
  \pre 0 <= eid << p_dcr->m_NumElements
*/
template <typename EIT, typename BVT>
void GEBV_MeshSolidShape2_BSlab( const MeshSolidShape2* p_mesh, const DCR_MeshSolidShape2* p_dcr, const Transform2& tr, const MeshSolidShape2::sdof_type* vec_node_pos,
                                 EIT eid, BVT& bv );

/*! Safe version of GEBV_MeshSolidShape2_BDOP */
template <typename EIT, typename BVT>
void GEBV_MeshSolidShape2_BDOP_Safe( const MeshSolidShape2* p_mesh, const Transform2& tr, const MeshSolidShape2::sdof_type* vec_node_pos,
                                     EIT eid, BVT& bv )
{
    // Check if safe, fall-back to BV(E) if not
    const DCR_MeshSolidShape2* pDCR( p_mesh->GetDCR() );
    if( !pDCR
        || eid < 0
        || eid >= pDCR->m_NumElements ) { GEBV_MeshSolidShape2_E<EIT,BVT>( p_mesh, tr, vec_node_pos, eid, bv ); return; }
    const DCR_MeshSolidShape2::ElementData& ed( pDCR->m_vecED[eid] );
    if( ed.m_BDOP[0].IsEmpty()
        || ed.m_BDOP[1].IsEmpty()
        || ed.m_BDOP[2].IsEmpty() ) { GEBV_MeshSolidShape2_E<EIT,BVT>( p_mesh, tr, vec_node_pos, eid, bv ); return; }
    // Call fast version if safe
    GEBV_MeshSolidShape2_BDOP<EIT,BVT>( p_mesh, pDCR, tr, vec_node_pos, eid, bv );
}

/*! Safe version of GEBV_MeshSolidShape2_BSlab */
template <typename EIT, typename BVT>
void GEBV_MeshSolidShape2_BSlab_Safe( const MeshSolidShape2* p_mesh, const Transform2& tr, const MeshSolidShape2::sdof_type* vec_node_pos,
                                      EIT eid, BVT& bv )
{
    // Check if safe, fall-back to BV(E) if not
    const DCR_MeshSolidShape2* pDCR( p_mesh->GetDCR() );
    if( !pDCR
        || eid < 0
        || eid >= pDCR->m_NumElements ) { GEBV_MeshSolidShape2_E<EIT,BVT>( p_mesh, tr, vec_node_pos, eid, bv ); return; }
    const DCR_MeshSolidShape2::ElementData& ed( pDCR->m_vecED[eid] );
    if( ed.m_BDOP[0].IsEmpty()
        || ed.m_BDOP[1].IsEmpty()
        || ed.m_BDOP[2].IsEmpty() ) { GEBV_MeshSolidShape2_E<EIT,BVT>( p_mesh, tr, vec_node_pos, eid, bv ); return; }
    // Call fast version if safe
    GEBV_MeshSolidShape2_BSlab<EIT,BVT>( p_mesh, pDCR, tr, vec_node_pos, eid, bv );
}
#  endif //__GEO_TRISS_ENABLE_BVH_TODO
#endif //__GEO_TRISS_ENABLE_BVH

} //namespace geo

#endif // GEO_SHAPE_TRISURFACESHAPE3_H
