#ifndef GEO_SHAPE_MESHSOLIDSHAPE2_H
#define GEO_SHAPE_MESHSOLIDSHAPE2_H

#include <Geo/shape/IShape.h>
#include <Geo/bv/bv.h> // BV types for ComputeBV()
#include <Geo/bv/GBDOP.h>
#include <vector>

//\todo CONNECT HE_Sym to BP edges?? This is UNTESTED, may change iterator behaviour, don't enable it without checking it thoroughly
//#define __ENABLE_LINK_BOUNDARY_SYMMETRIC_HE

#define __GEO_MSS_ENABLE_BVH
#ifdef __GEO_MSS_ENABLE_BVH
#  include <Geo/bv/GBoundingVolumeHierarchy.h> //TEMPORAL: quick and dirty BVH
#endif

#define __ENABLE_MSS_DCR_PATCH_TOPOLOGY //should be enabled by default
#define __ENABLE_MSS_DCR_PATCH_BDT //should be enabled by default

// dimension-independent mesh topology
namespace geo
{

#ifdef __GEO_MSS_ENABLE_BVH
//TEMPORAL: Consider generic IBVH so that each MSS can have a user-specified or best fitting one
//typedef bv::GBoundingVolumeHierarchy_ST_DG<bv::AABB2,feature_index_type> BVH_MeshSolidShape2;
typedef bv::GBoundingVolumeHierarchy_ST_DG<bv::DOP2_K8,feature_index_type> BVH_MeshSolidShape2;
//typedef bv::GBoundingVolumeHierarchy_ST_DG<bv::Sphere2,feature_index_type> BVH_MeshSolidShape2;
#endif

struct vertex_topology
{
    /*
    vertex_topology() : m_OutHEID(cInvalidFeatureIndex), m_NumEdges(0) {}
    vertex_topology( feature_index_type out_heid, feature_count_type num_edges ) : m_OutHEID(out_heid), m_NumEdges(num_edges) {}
    */
    feature_index_type m_OutHEID;
    feature_count_type m_NumEdges; //! Redundant but useful
};
struct polygon_topology
{
    /*
    polygon_topology() : m_FirstHEID(cInvalidFeatureIndex), m_NumEdges(0) {}
    polygon_topology( feature_index_type first_heid, feature_count_type num_edges ) : m_FirstHEID(first_heid), m_NumEdges(num_edges) {}
    */
    feature_index_type m_FirstHEID;
    feature_count_type m_NumEdges; //!< Redundant but useful
};
struct half_edge_topology
{
    /*
    half_edge_topology()
    : m_OriginVID(cInvalidFeatureIndex), m_NextHEID(cInvalidFeatureIndex), m_SymHEID(cInvalidFeatureIndex), m_LeftPId(cInvalidFeatureIndex) {}
    half_edge_topology( feature_index_type origin_vid, feature_index_type next_heid, feature_index_type sym_heid, feature_index_type right_pid )
    : m_OriginVID(origin_vid), m_NextHEID(next_heid), m_SymHEID(sym_heid), m_LeftPId(right_pid) {}
    */
    feature_index_type m_OriginVID;
    feature_index_type m_NextHEID;
    feature_index_type m_SymHEID;
    feature_index_type m_LeftPID;
};

} //namespace geo

namespace geo
{

//----------------------------------------------------------------
// DCR_MeshSolidShape2 \todo Move elsewhere!
//----------------------------------------------------------------

class IDCR2
{
public:
    IDCR2() {}
    virtual ~IDCR2() {}
};

/*! Detailed Collision Representation for MSS2

  For each MSS2 polygon/element E_i, DCR stores a compact array of
  Vertices and Edges contained in E_i, completely decoupled from other
  E_j.
  - All V inside a E_i can be transformed at once, independently of other E_j
    => V at the boundary between E_i are DUPLICATED to ensure this property.
  - All E inside a E_i can be iterated independently of other E_j
  - Topology
    - Is LOCAL topology IS REQUIRED?
      => YES, for CD we need to navigate the boundary surface (eg:
         surface point generation, intersection mapping...), otherwise
         CD algorithms would need to filter redundant info (repeated
         edges, vertices...) at runtime.
      => This DISALLOWS topology-less flattened feature storage (ex:
         3vtx per triangle, 2vtx per edge), which is better for
         HW acceleration in some contexts.
    - Is GLOBAL cross-patch topology also REQUIRED?
      => YES. Even if we use MSS/TSS elements as BV for DCR geometry,
         CD methods should be able to cross per-element patch
         boundaries into neighbour elements/patches, specially if only
         MSS/TSS Crust elements are considered in BV overlapping tests
         (otherwise, it would be possible for a CD algorithm to be
         unaware of a whole Crust element being completely inside
         another TSS/MSS, without overlapping its Crust. Intersection
         mapping, for example, must be free to navigate the TSS/MSS
         boundary for this completely enclosed elements/patches)
         => This DISALLOWS storing per-patch topology using local
            per-patch indices (relative to patch base VID/EID) as
            smaller integral type (ex: uint8 allows 256 V...), which
            would otherwise be a nice optimization to save memory.
    => MOREOVER, CD navigation NEEDS to be aware of element crossing
       if per-element vertices are transformed lazily. This is EASY if
       global feature indexing is used, as neighbour indices can be
       compared with "current Element" global vertex index range
       [first, first+count) and, if outside, transform all vertices in
       the newly entered element.
       => Neighbour elements should be checked first for matching, but
          a flood may be necessary if CD navigation can cross several
          elements at once.

  - \todo We WILL NOT assume that each E_i only contains a single
    connected boundary patch, but we'll store the features
    sequentially in sub-arrays per-patch, and add some per-patch
    iteration functionality. Ideally 1 Element <=> 1 connected Patch,
    but this needs to be enforced at TSS/MSS generation, not here..

  - \todo If we guarantee that m_NumPatches = 1, the PatchData attribs
    could be embedded into ElementData and remove m_vecPD.

  - DCR geometry is stored relative to MSS DOF, and affected at
    runtime by MSS object transform and deformable DOF.

  - Deformed DCR geometry is NOT PERSISTENT, and could be stored in
    np::Context::m_ScratchPad during the specific Contact/Overlap
    function.

  - \todo Consider storing per-patch bounding-volume (BV) in element
    (baricentric) coordinates. Only bv::GSimplex<D> BV preserves
    complexity when deformed baricentrically (AABB become general
    quads/hexadrons, KDOP even worse, spheres become distorted, etc..).
    - Computing the "smallest enclosing simplex" (SES) for a cloud of
      points is not trivial. Even worse, the SES can extend outside
      the element, which invalidates its baricentric deformation.
    - Testing simplex pairs is more expensive than testing
      KDOP/Spheres
    => Convexity is preserved by baricentric deformation, therefore we
       could use the exact or *simplified* convex-hull (CH) of the patch,
       instead of the SES, which is guaranteed to be inside the element.
       - Storing arbitrary per-patch CH is troublesome due to dynamic
         sizes, but achievable with a first_ch_vid/num_ch_vid
         structure and a global ch_vertices array
       - Deformed CH vs CH overlap test is expensive for nontrivial
         CH... unless we use GJK with some persistence, which further
         complicates every-fucking-thing....
    => BKDOPS (Barycentric-KDOPS)
       - KDOP in barycentric coords... why not?!
       - Baricentric coord axis should be enough (3 in 2D, 4 in 3D)
         => BDOP (with fixed K=2*(Dimension+1))
       - Deformed BDOPS can be seen as KDOPS with different/deformed
         Directions, but affine/baricentric transforms preserve linear
         geometry (points, segments, triangles, etc...)
       - If we compute the baricentric vertices of BDOP_1 and BDOP_2,
         we CAN use SAT to test them.
         - In 2D, just per-Direction tests (3+3) are required
           (corresponding to BDOP_i face normals). At most 6 points
           need to be stored for each BDOP_i.
         - In 3D, cross-edge-pair directions (A LOT) would be required
           for an EXACT test, but for early-out BV purposes these can
           be avoided, thus, per-Direction (4+4) tests are enough. At
           most 14 points need to be stored for each BDOP_i.
       - IMPORTANT: As the axis/directions are the SAME as the
         triangle/element edge/face normals, the BDOP pair test can
         include the exact element pair test as an early out test
         (using the 3 or 4 element vertices in the SAT) or directly
         substitute it with tighter bounds.
       - Using affine transforms from e1 to e2 (that include
         translations, rotations and baricentric deformations) we can
         test e1 linear geometry features f1 (points, segments,
         triangles) against e2 BDOP by computing their BDOP locally in
         e2, which should discard most feature1 vs feature2
         pairs... However, this would require transforming f1 vertices
         to ALL potentially overlapping e2 local coordinates.
         => For stochastic this would not be that expensive, because
            points need to be generated and transformed on the fly
            anyway and no global caching of transformed vertices
            should be applied.
       ==> All in all, BDOPs are a generalization of AABB, NOT KDOP,
           in barycentric coords, as they do not use arbitrary
           directions but only coordinate axis. All range tests become
           single barycentric coordinate tests with no need for
           dot-product projection ops along non-coordinate axis.
       ==> The paper ""2007 Adaptive Deformations with Fast Tight
           Bounds" did something similar, but uses the BDOP to compute
           a conservative but tight AABB in global coords that bounds
           the embedded geometry, not the element.
*/
class DCR_MeshSolidShape2: public IDCR2
{
public:
    enum EConstants
    {
        cVersion = 4 // 4) Added patch topology and BDOP-Tree officially
        //cVersion = 3 // 3) Renamed attributes
        //cVersion = 2 // 2) Added m_BDOP_BestSlabIdx
        //cVersion = 1 // 1) Added m_BDOP
        //cVersion = 0
    };
public:
    struct Element
    {
        // Sub-array of global patches in this element, unconnected
        feature_index_type m_FirstPID;
        feature_count_type m_NumPatches;
        // Sub-array of global V in this element, to transform all at once
        uint32 m_FirstVID;
        uint32 m_NumVertices; //\todo Consider 12:20 format for a single uint32
        // Barycentric-DOP
        bv::BDOP3 m_BDOP;
        uint32 m_BDOP_BestSlabIdx; //2b would be enough!
        /* Barycentric-DOP
        Vec2 m_BDOP_Vertices[6];
        */
    };
    struct Patch
    {
        uint32 m_EID; //Useful to identify the DCR.E when navigating to neighbour patches, faster than searching new EID, even if topologically
        // Sub-array of global S in this patch
        uint32 m_FirstSID;
        uint32 m_NumSegments; //\todo Consider 12:20 format for a single uint32
#ifdef __ENABLE_MSS_DCR_PATCH_TOPOLOGY
        // Topology IMPORTANT m_vecSID[0,1] != m_FirstSID,m_FirstSID+m_NumSegments-1 (sub-array order is NOT NECESSARILY topological CCW order if BDT is used)
        inline Patch() { m_EID = 0xFFFFFFFF; m_FirstSID = 0xFFFFFFFF; m_NumSegments = 0xFFFFFFFF; m_vecSID[0] = 0xFFFFFFFF; m_vecSID[1] = 0xFFFFFFFF; m_vecNPID[0] = 0xFFFFFFFF; m_vecNPID[1] = 0xFFFFFFFF; } //TEMP: Fill to test!
        uint32 m_vecSID[2]; //begin/end CCW SID, this IS NOT REQUIRED, but allows skipping the whole Patch, possible due to explicit CCW order (this is NOT possible in DCR_TetSolidShape3)
        uint32 m_vecNPID[2]; //prev/next CCW PID: There can be either 0 or 2 NPID (0=closed curve, 2=4-way tube with on exit through each tetrahedron face, mark with -1 invalid entries
#else
        uint32 m_DUMMY[4];
#endif

#ifdef __ENABLE_MSS_DCR_PATCH_BDT
        typedef std::pair<uint32,uint32> bdt_node_sip_range;
#endif
        //\todo Store per-patch aggregate data (ex: avgpos, avgnormal|integrated N*dA, total_length...)
    };
    struct Segment
    {
        /*\todo VID can be relative to DCR.E.m_FirstVID, thus, 16b should be enough
          Layout:
          - uint32 m_vecVID[2];
          - uint32 m_vecSID[2];

          Discussion:
          - ID could use only 20b (1M), freeing 4x12b = 48b
          - a) We can use 6x8b = 48b to store 6 quantized BDOP min/max
            - The split is implicit at N/2, but DCR.F sorting and split axis are ARBITRARY, no need to store it at all as we save whole BDOP
            - This forces clipping rays against BDOP[3] at each node, and L/R sorting would not be possible (well, it WOULD be possible for a given BDOP axis)
          - b) We can use 2x8b to store 2 q-BSlab min/max +2b to store BSlab axis +20b to store R-offset (L=P+1) => 16+2+20=38b
            - Easier to factor 2x8b q-BSlab + 2b axis + 14b offset => Extract upper 4x8b
            - c) SIMPLER: We can use only the upper 4x8b: 2x15b to store high-res q-BSlab min/max and 2b to store the axis. VID/NSID will be 24b then, more than enough!
          - We do NOT NEED to flag leaf-nodes specifically, the tree-traversal can decide to test "all features in the current subtree" at once at ANY POINT, for example, when #F is below a certain threshold.
          - \note Min/Max intervals CANNOT BE EMPTY, but cannot use GInterval either (it requires signed<T>)

          - IMPORTANT: If we only save BDOP[best_axis] in a BSTN, and
            each L/R child can select its own best_axis, THERE IS NO
            L/R culling from the single parent BSlab!! Either we
            retrieve L/R children and compute their "split range", or
            we're DOOMED...  => TriTree stores separate rangeNode and
            rangeSplit
            => Thus, either we save the L/R split_range in each BSTN
               parent node OR we store the whole BDOP so that L/R child
               will have updated intervals in ALL axis, not just their
               best one
        */
        /*
          Overall, saving the whole BDOP in each BDTNode seems the simplest and most efficient solution
          - BSTN are self-contained
          - Uniform data representation
          - RayVsBDOP only tested once per node, clipped and propagated to children
        */
        inline Segment() //TEMP: Fill it with invalid to check BST, remove when sufficiently tested!
            {
                m_vecVID[0] = 0xFFFFFFFF;
                m_vecVID[1] = 0xFFFFFFFF;
                m_vecNSID[0] = 0xFFFFFFFF;
                m_vecNSID[1] = 0xFFFFFFFF;
            }

#ifdef __ENABLE_MSS_DCR_PATCH_BDT
        enum EConstants { eBDT_IntervalQuantizationBits = 8,
                          eBDT_MaskVID = 0x000FFFFF,
                          eBDT_MaskNSID = 0x000FFFFF };
        inline uint32 GetVID( int i ) const { return m_vecVID[i] & eBDT_MaskVID; }
        inline uint32 GetNSID( int i ) const { return m_vecNSID[i] & eBDT_MaskNSID; }
        inline void SetVID( int i, uint32 vid ) { GEO_ASSERT(0 == (vid & ~eBDT_MaskVID)); m_vecVID[i] = (m_vecVID[i] & ~eBDT_MaskVID) | (vid & eBDT_MaskVID); }
        inline void SetNSID( int i, uint32 nsid ) { GEO_ASSERT(0 == (nsid & ~eBDT_MaskNSID)); m_vecNSID[i] = (m_vecNSID[i] & ~eBDT_MaskNSID) | (nsid & eBDT_MaskNSID); }
        //\name BarycentricDOPTreeNode
        //@{
        /* Quantized BDOP3 with 8b quantization but int32 interval
           values, as GInterval REQUIRES signed T
        */
        class BDOP3q
        {
        public:
            typedef mal::GInterval<int32> interval_q_type;
            inline BDOP3q() {}
            // Init from BDOP
            inline BDOP3q( const bv::BDOP3& bdop3 )
                {
                    for( int i=0; i<3; i++ )
                        m_vecInterval[i].SetMinMax( mal::Clamp<int32>( mal::Floor(mal::Clamp01(bdop3[i].Min())*0x000000FF), 0, 0x000000FF ),
                                                    mal::Clamp<int32>( mal::Ceil(mal::Clamp01(bdop3[i].Max())*0x000000FF), 0, 0x000000FF ) );
                }
            // Cast to BDOP
            inline operator bv::BDOP3() const
                {
                    bv::BDOP3 bdop3;
                    const Real q2f( mal::Rcp<Real>( (1<<eBDT_IntervalQuantizationBits)-1 ) );
                    for( int i=0; i<3; i++ )
                        bdop3[i].SetMinMax( q2f*m_vecInterval[i].Min(), q2f*m_vecInterval[i].Max() );
                    return bdop3;
                }
            const interval_q_type& operator[](int i) const { return m_vecInterval[i]; }
            interval_q_type& operator[](int i) { return m_vecInterval[i]; }
        private:
            interval_q_type m_vecInterval[3];
        };

        inline BDOP3q BDTN_BDOPq() const
            {
                BDOP3q bdop3q;
                bdop3q[0].SetMinMax( (m_vecVID[0] & 0xFF000000) >> 24, (m_vecVID[1] & 0xFF000000) >> 24 );
                bdop3q[1].SetMinMax( (m_vecNSID[0] & 0xFF000000) >> 24, (m_vecNSID[1] & 0xFF000000) >> 24 );
                bdop3q[2].SetMinMax( ((m_vecVID[0] & 0x00F00000) >> 16) | ((m_vecVID[1] & 0x00F00000) >> 20),
                                     ((m_vecNSID[0] & 0x00F00000) >> 16) | ((m_vecNSID[1] & 0x00F00000) >> 20) );
                return bdop3q;
            }

        inline void BDTN_Init( const bv::BDOP3& bdop3 )
            {
                BDOP3q bdop3q( bdop3 );
                // Save bdop3q into the right bits...
                // BDOP3q axis 0,1 min/max go into upper 8b
                m_vecVID[0] = (m_vecVID[0] & eBDT_MaskVID) | ((bdop3q[0].Min() & 0x000000FF) << 24);
                m_vecVID[1] = (m_vecVID[1] & eBDT_MaskVID) | ((bdop3q[0].Max() & 0x000000FF) << 24);
                m_vecNSID[0] = (m_vecNSID[0] & eBDT_MaskNSID) | ((bdop3q[1].Min() & 0x000000FF) << 24);
                m_vecNSID[1] = (m_vecNSID[1] & eBDT_MaskNSID) | ((bdop3q[1].Max() & 0x000000FF) << 24);
                // BDOP3q axis 2 min/max are in 4b packets at the [20..24) range spread accross vecVID[2] and vecNSID[2]
                m_vecVID[0] = (m_vecVID[0] & 0xFF0FFFFF) | ((bdop3q[2].Min() & 0x000000F0) << 16);
                m_vecVID[1] = (m_vecVID[1] & 0xFF0FFFFF) | ((bdop3q[2].Min() & 0x0000000F) << 20);
                m_vecNSID[0] = (m_vecNSID[0] & 0xFF0FFFFF) | ((bdop3q[2].Max() & 0x000000F0) << 16);
                m_vecNSID[1] = (m_vecNSID[1] & 0xFF0FFFFF) | ((bdop3q[2].Max() & 0x0000000F) << 20);
            }
        //@}
#else
        inline uint32 GetVID( int i ) const { return m_vecVID[i]; }
        inline uint32 GetNSID( int i ) const { return m_vecNSID[i]; }
        inline void SetVID( int i, uint32 vid ) { m_vecVID[i] = vid; }
        inline void SetNSID( int i, uint32 nsid ) { m_vecNSID[i] = nsid; }
#endif

    private:
        uint32 m_vecVID[2];
        uint32 m_vecNSID[2]; //prev,next \todo NECESSARY, as extreme features in an element/patch MUST be connected to neighbour elements/patches
    };

    // const IObject2* m_pObject; //\todo I'm not sure... this would be like embedding, as the MSS_DCR actual DOF depend on the actual MSS DOF, but DCR can never be standalone...

    uint32 m_NumElements;
    uint32 m_NumPatches;
    uint32 m_NumSegments;
    uint32 m_NumVertices;

    Element* m_vecE; //Indices parallel to MSS elements \todo Should only contain entries layer[0] if the MSS is properly cropped
    Patch* m_vecP; //stores all patches in per-element contiguous subarrays, NOT in global topological order
    Segment* m_vecS;
    Vec2* m_vecV; // Relative to MSS tr/dof

    //\todo Store as uint8* m_pBuffer if shared is supported

    DCR_MeshSolidShape2()
    : m_NumElements(0), m_NumPatches(0), m_NumSegments(0), m_NumVertices(0)
    , m_vecE(0), m_vecP(0), m_vecS(0), m_vecV(0)
        {}
    ~DCR_MeshSolidShape2()
        {
            if(m_vecE) delete[] m_vecE; m_vecE = 0;
            if(m_vecP) delete[] m_vecP; m_vecP = 0;
            if(m_vecS) delete[] m_vecS; m_vecS = 0;
            if(m_vecV) delete[] m_vecV; m_vecV = 0;
        }

    void SetBakedData( bool b_shared,
                       uint32 num_e, uint32 num_p, uint32 num_s, uint32 num_v,
                       const Element* vec_ed, const Patch* vec_pd, const Segment* vec_s, const Vec2* vec_v )
        {
            if( b_shared ) GEO_LOG_WARNING("DCR_MeshSolidShape2::SetBakedData ignores b_shared by now...");
            m_NumElements = num_e;
            m_NumPatches = num_p;
            m_NumSegments = num_s;
            m_NumVertices = num_v;
            m_vecE = new Element[num_e];
            m_vecP = new Patch[num_p];
            m_vecS = new Segment[num_s];
            m_vecV = new Vec2[num_v];
            memcpy( m_vecE, vec_ed, sizeof(Element)*num_e );
            memcpy( m_vecP, vec_pd, sizeof(Patch)*num_p );
            memcpy( m_vecS, vec_s, sizeof(Segment)*num_s );
            memcpy( m_vecV, vec_v, sizeof(Vec2)*num_v );
        }
};

//----------------------------------------------------------------
// MeshSolidShape2
//----------------------------------------------------------------

/* 2D Polygon mesh representation.
   \todo Derive OR build a TriSolidShape2 class restricted to 2D
         triangulated polygonals, similar to TetSolidShape3
         => From a MeshSolidShape2 a compact TriSolidShape2 could be
            built, storing only Vertices and Faces

   \todo It should be possible to have the topology-only mesh
         representation dimension-independent and the geometry
         representation (V,normals, etc...) in dimension-specific
         subclasses or instances. CGAL has some stuff to do this, but
         too much abstraction will hurt speed, so FORGET IT.
*/
class MeshSolidShape2: public IShape2
{
public:
    enum EConstants
    {
        cVersion = 1 //Initial version, added together with Layers functionality
        //cVersion = 0 //FAKE version to keep compatibility with existing data, must remove and reexport when done with layers
    };

public:
    //! \name Mandatory generic interface for IShape
    //@{
    static const unsigned int cDimension = 2;
    typedef mal::GTransform<Real,2> transform_type;
    typedef mal::GVec<Real,2> sdof_type;
    typedef mal::GVec<Real,2> vec_type;
    //@}

public:
    struct vertex_type: public vertex_topology {}; //!< \note Could add 2d geometry data here
    struct polygon_type: public polygon_topology {}; //!< \note Could add 2d geometry data here
    struct half_edge_type: public half_edge_topology {}; //!< \note Could add 2d geometry data here
    struct polygon_layer_type
    {
        feature_index_type m_FirstPID;
        feature_count_type m_NumPolygons;
    };

public:
    MeshSolidShape2();
    ~MeshSolidShape2();

    //!\name Dynamic IShape implementation
    //@{
    EShapeType GetType() const { return eShape_MeshSolid2; }
    unsigned int GetNumSDOF() const { return m_NumV; }
    const sdof_type *GetVecDefaultSDOF() const { return m_vecPoints; }
    const Real *GetVecDefaultDOF() const { return reinterpret_cast<const Real*>(m_vecPoints); }
    void ComputeBVD( bv::BoundingVolume2& bv, const transform_type& transform, const sdof_type *vec_sdof ) const;

    IDomainSampler2* CreateDomainSampler() const;
    //@}

    //! Init from external arrays
    void SetBakedData( bool b_shared,
                       uint32 num_v, uint32 num_p, uint32 num_he, uint32 num_boundary_p, uint32 num_boundary_he, uint32 num_layers,
                       const Vec2* vec_points, const vertex_type* vec_v, const polygon_type* vec_p, const half_edge_type* vec_he, const polygon_layer_type* vec_l );
    //TEMPORAL: hacks to add annotations to non-editable MSS
    inline void SetBakedDCR_StrictlyNonshared_UglyHack( DCR_MeshSolidShape2* p_nonshared_dcr ) { m_pDCR = p_nonshared_dcr; } //never shared, breaks SetBakedData() strict const enforcement
#ifdef __GEO_MSS_ENABLE_BVH
    inline void SetBakedBVH_StrictlyNonshared_UglyHack( BVH_MeshSolidShape2* p_nonshared_bvh ) { m_pBVH = p_nonshared_bvh; } //never shared, breaks SetBakedData() strict const enforcement
#endif

    inline uint32 GetNumV() const { return m_NumV; }
    inline uint32 GetNumP() const { return m_NumP; }
    inline uint32 GetNumHE() const { return m_NumHE; }

    inline uint32 GetNumBoundaryP() const { return m_NumBoundaryP; }
    inline uint32 GetNumBoundaryHE() const { return m_NumBoundaryHE; }
    inline uint32 GetNumL() const { return m_NumL; }

    inline uint32 GetNumAllocV() const { return m_NumAllocV; }
    inline uint32 GetNumAllocP() const { return m_NumAllocP; }
    inline uint32 GetNumAllocHE() const { return m_NumAllocHE; }

    inline const Vec2 *GetVecPoints() const { return m_vecPoints; }
    inline const vertex_type *GetVecV() const { return m_vecV; }
    inline const polygon_type *GetVecP() const { return m_vecP; }
    inline const half_edge_type *GetVecHE() const { return m_vecHE; }
    inline const polygon_layer_type *GetVecL() const { return m_vecL; }

    // \todo
    inline bool IsClosed() const { return true; } //!< Must be closed by construction
    inline bool IsSimplyConnected() const { return true; } //!< Check! false if there are holes
    inline bool IsConnected() const { return true; } //!< Check!
    inline bool IsConvex() const { return false; } //!< Not so useful for solid meshes...

    // Vertex access
    inline Vec2 V_Pos_0( uint32 vid ) const { return m_vecPoints[vid]; }
    inline Vec2 V_Pos( uint32 vid, const sdof_type *vec_sdof ) const { return vec_sdof[vid]; }
    inline uint32 V_OutHEID( uint32 vid ) const { return m_vecV[vid].m_OutHEID; }
    inline uint32 V_NumEdges( uint32 vid ) const { return m_vecV[vid].m_NumEdges; }

    // Polygon access
    inline uint32 P_FirstHEID( uint32 pid ) const { return m_vecP[pid].m_FirstHEID; }
    inline uint32 P_NumEdges( uint32 pid ) const { return m_vecP[pid].m_NumEdges; }
    Vec2 P_Barycenter_0( uint32 pid ) const;
    Vec2 P_Barycenter( uint32 pid, const sdof_type *vec_sdof ) const;
    uint32 P_VecVID( uint32 pid, uint32* vec_vid, uint32 max_vid ) const;

    // Boundary polygon access
    inline uint32 BP_FirstHEID( uint32 bpid ) const { return P_FirstHEID( m_NumP + bpid ); }
    inline uint32 BP_NumEdges( uint32 bpid ) const { return P_NumEdges( m_NumP + bpid ); }
    Vec2 BP_Barycenter_0( uint32 bpid ) const { return P_Barycenter_0( m_NumP + bpid ); }
    Vec2 BP_Barycenter( uint32 bpid, const sdof_type *vec_sdof ) const{ return P_Barycenter( m_NumP + bpid, vec_sdof ); }

    // Layers access
    inline uint32 L_FirstPID( unsigned int lid ) const { return m_vecL[lid].m_FirstPID; }
    inline uint32 L_NumP( unsigned int lid ) const { return m_vecL[lid].m_NumPolygons; }

    // HalfEdge access (\note Works both for P and BP)
    inline uint32 HE_Next( uint32 heid ) const { return m_vecHE[heid].m_NextHEID; }
    inline uint32 HE_Prev( uint32 heid ) const {
        /*\note Assuming consecutive HE along a P perimeter in CCW
          order. Otherwise it would require search for Prev along face
          HE, very expensive O(N_e)
        */
        uint32 left_pid( HE_LeftPID(heid) );
        uint32 num_edges( P_NumEdges(left_pid) );
        uint32 heid0( P_FirstHEID(left_pid) );
        return heid0 + ((heid - heid0 - 1 + num_edges) % num_edges );
    }
    inline uint32 HE_Sym( uint32 heid ) const { return m_vecHE[heid].m_SymHEID; }
    inline uint32 HE_OriginVID( uint32 heid ) const { return m_vecHE[heid].m_OriginVID; }
    inline uint32 HE_FinalVID( uint32 heid ) const { return HE_OriginVID( HE_Next(heid) ); }
    inline uint32 HE_LeftPID( uint32 heid ) const { return m_vecHE[heid].m_LeftPID; }
    inline uint32 HE_RightPID( uint32 heid ) const { return ( cInvalidFeatureIndex != HE_Sym(heid) )
                                                            ? HE_LeftPID( HE_Sym(heid) )
                                                            : cInvalidFeatureIndex; }
    inline Vec2 HE_Direction_0( uint32 heid ) const { return mal::Normalized( V_Pos_0(HE_FinalVID(heid)) - V_Pos_0(HE_OriginVID(heid)) ); }
    inline Vec2 HE_Direction( uint32 heid, const sdof_type *vec_sdof ) const
        { return mal::Normalized( V_Pos(HE_FinalVID(heid),vec_sdof) - V_Pos(HE_OriginVID(heid),vec_sdof) ); }


/* \todo replicate this functionality but prettier :D
    inline Vec2 P_Normal( uint32 pid ) const { return Vec2(m_vecP[pid].m_Normal); }
    inline Real P_CoeffD( uint32 pid ) const { return Real(m_vecP[pid].m_CoeffD); }
    inline Vec2 P_Point0( uint32 pid ) const { return V_Pos( HE_OriginVID( P_FirstHEID( pid ) ) ); }
    inline Real P_Depth( uint32 pid ) const { return m_vecP[pid].m_Depth; }
*/

    // Get annotations or 0 if none
    /*\todo Annotations are built for specific shape params and
     default DOF, which they can refer implicitly assuming constant
     indices/values, if annotation-specific data is required to deal
     with runtime DOF, this annotation-DOF should be allocated by
     IObject. Ex: Transformed/Deformable vertex positions in a DCR
    */
    inline const DCR_MeshSolidShape2* GetDCR() const { return m_pDCR; }
#ifdef __GEO_MSS_ENABLE_BVH
    inline BVH_MeshSolidShape2* GetBVH() const { return m_pBVH; } //\note non-const until topology/geometry are split, because we need to refit it!
#endif

public:
    /*! Iterates over all vertex-incident Polygons in CCW (++,forward)
        or CW (--,backwards) order:

        Perimeter vertices can only be iterated CCW

        Iteration finishes whenever the a loop is closed (back to
        initial PID) or an open edge is found.

        The current HEID() is Outgoing from VID(), thus,
        HE_FinalVID(HEID()) is the Final VID

        \notes
        it_ccw++ => m_CurrentHEID = Sym( Prev( m_CurrentHEID ) )
        it_ccw-- => m_CurrentHEID = Next( Sym( m_CurrentHEID ) )
        it_ccw.HEID() = m_CurrentHEID
        it_ccw.PID() = HE_LeftPID( m_CurrentHEID ) )
    */
    class iterator_polygons_around_vertex_ccw
    {
    public:
        finline iterator_polygons_around_vertex_ccw( const MeshSolidShape2 *p_mss, uint32 vid )
        : m_pMSS(p_mss), m_VID(vid), m_IsValid(true), m_Type(eUnknown)
        {
            m_HEID = m_pMSS->V_OutHEID(m_VID);
        }

        finline uint32 VID() const { return m_VID; } //\note This is the central vid, not iterated vid!!
        finline uint32 HEID() const { return m_HEID; }
        finline uint32 PID() const { return m_pMSS->HE_LeftPID(m_HEID); }

        finline bool IsValid() const { return m_IsValid; }
        finline bool IsOpen() const { return eOpen == m_Type; }

        inline void operator++()
        {
            GEO_ASSERT( IsValid() );
            uint32 next_heid_ccw = m_pMSS->HE_Sym( m_pMSS->HE_Prev(m_HEID) );
            if( next_heid_ccw == cInvalidFeatureIndex ) //open?
            {
                m_Type = eOpen;
                m_IsValid = false;
            }
            else if( next_heid_ccw == m_pMSS->V_OutHEID(m_VID) ) //last?
            {
                m_Type = eClosed;
                m_IsValid = false;
            }
            else //advance!
                m_HEID = next_heid_ccw;
        }
        inline void operator--()
        {
            GEO_ASSERT( IsValid() );
            uint32 sym_heid = m_pMSS->HE_Sym( m_HEID );
            if( sym_heid == cInvalidFeatureIndex )
            {
                m_Type = eOpen;
                m_IsValid = false;
            }
            else
            {
                uint32 next_heid_cw = m_pMSS->HE_Next(sym_heid);
                if( next_heid_cw == m_pMSS->V_OutHEID(m_VID) ) //last?
                {
                    m_Type = eClosed;
                    m_IsValid = false;
                }
                else //advance!
                    m_HEID = next_heid_cw;
            }
        }
    private:
        const MeshSolidShape2 *m_pMSS;
        uint32 m_VID;
        uint32 m_HEID;
        bool m_IsValid;
        enum EType { eUnknown = 0, eOpen, eClosed } m_Type; //Open edge found while iterating
    };
    iterator_polygons_around_vertex_ccw GetIterator_PolygonsAroundVertexCCW( uint32 vid ) const { return iterator_polygons_around_vertex_ccw(this,vid); }
    //@}

    //\todo iterator_boundary: Iterate around a given boundary polygon, retrieve VID, HEID and PID

protected:
    void ClearBakedData();

protected:

    uint32 m_NumV;
    uint32 m_NumP;
    uint32 m_NumHE;

    uint32 m_NumBoundaryP;
    uint32 m_NumBoundaryHE;

    uint32 m_NumL;

    uint32 m_NumAllocV;  //= m_NumV
    uint32 m_NumAllocP;  //= m_NumP + m_NumBoundaryP
    uint32 m_NumAllocHE; //= m_NumHE + m_NumBoundaryHE

    const Vec2* m_vecPoints;
    const vertex_type* m_vecV;
    const polygon_type* m_vecP;
    const half_edge_type* m_vecHE;
    const polygon_layer_type* m_vecL;

    uint32* m_pBuffer; //!< If exists, the shape is NOT SHARED

    // Annotations
    DCR_MeshSolidShape2* m_pDCR; //strictly nonshared, by now...
#ifdef __GEO_MSS_ENABLE_BVH
    BVH_MeshSolidShape2* m_pBVH; //strictly nonshared, by now...
#endif

private:
    class DomainSampler;
};

/* Editable Mesh Solid Shape 2D
   \todo All ops should work with half-edge representation directly.
*/
class EditableMeshSolidShape2: public MeshSolidShape2
{
public:
    void operator=( const EditableMeshSolidShape2& emss ) = delete;
    EditableMeshSolidShape2( const EditableMeshSolidShape2& emss ) = delete;
public:
    EditableMeshSolidShape2();
    ~EditableMeshSolidShape2();

    void Clear();
    void Set( const MeshSolidShape2& mss2 );

    void BeginEdition(); //\note Incremental, preserves BakedData if any, call Clear() to reset completely
    bool EndEdition();
    inline bool IsBeingEdited() const { return m_IsBeingEdited; }

    feature_index_type AddVertex( const Vec2& point );
    feature_index_type AddPolygon3( feature_index_type vid0, feature_index_type vid1, feature_index_type vid2 );
    feature_index_type AddPolygonN( uint32 num_vid, feature_index_type *vec_vid );

    void SetVertex( uint32 vid, const Vec2& point );
    bool Subdivide( uint32 first_smooth_vid = 0 ); //1 Edge => 1 Tri, \todo Flags: None, LoopV, LoopE, ....
    //void Distort( Real magnitude, machine_uint_type seed = 666 ); //!< random-dir per-vertex distortion
    bool RemovePolygons( uint32 num_pid, const feature_index_type* vec_pid );

    // Add annotations
    bool AddDCR( const IShape2* p_shape, const Transform2& tr_s2mss );//, EDCRtype dcrt ); //Polygonal2, eDCR_Partition );
#ifdef __GEO_MSS_ENABLE_BVH
    bool AddBVH(); //bvh type!?
#endif
    /* \todo Add annotations: optional attributes, with automatic load/save if they exist
    bool AddDCR( DCR_MeshSolidShape2* p_dcr ); //This would decouple DCR creation from MSS API, allowing for ANY ext creation method... may be better than polymorphic AddDCR(IShape2)
    IBVH* AddBVH( BVHT, params... ); //Add BVH annotation
    */

private:

    void ClearEditData();
    void RebuildBakedData();
    bool FixDegeneracies();
    void RebuildTopology();
    void RebuildLayers();

    bool FixDegenerateVertices();
    bool FixDegenerateEdges(); //!< Collapse too-short and aligned edges
    bool FixDegeneratePolygons(); //!< Collapse too-small or pointy polygons
    feature_index_type FindV( const Vec2& pos, Real epsilon_sq ) const;
    feature_index_type FindHE( feature_index_type vid0, feature_index_type vid1 ) const;

private:
    struct editable_vertex_type: public vertex_type
    {
        editable_vertex_type() : m_Pos(0,0) { m_OutHEID = cInvalidFeatureIndex; m_NumEdges = 0; } //\note Required by m_addV.resize()
        editable_vertex_type( const Vec2& pos, feature_index_type out_heid, feature_count_type num_edges )
        : m_Pos(pos) { m_OutHEID = out_heid; m_NumEdges = num_edges; }
        Vec2 m_Pos;
    };
    struct editable_polygon_type: public polygon_type
    {
        editable_polygon_type( feature_index_type first_heid, feature_count_type num_edges )
            { m_FirstHEID = first_heid; m_NumEdges = num_edges; }
    };
    struct editable_half_edge_type: public half_edge_type
    {
        editable_half_edge_type( feature_index_type origin_vid, feature_index_type next_heid, feature_index_type sym_heid, feature_index_type left_pid )
            { m_OriginVID = origin_vid; m_NextHEID = next_heid; m_SymHEID = sym_heid; m_LeftPID = left_pid; }
    };

    bool m_IsBeingEdited;
    std::vector<editable_vertex_type> m_addV;
    std::vector<editable_polygon_type> m_addP;
    std::vector<editable_half_edge_type> m_addHE;

    std::vector<editable_polygon_type> m_addBP; //Boundary Polygons
    std::vector<editable_half_edge_type> m_addBHE; //Boundary HalfEdges

    struct editable_polygon_layer_type: public polygon_layer_type
    {
        editable_polygon_layer_type()
            { m_FirstPID = 0; m_NumPolygons = 0; } //\note Required by m_addL.resize()
        editable_polygon_layer_type( feature_index_type first_pid, feature_count_type num_polygons )
            { m_FirstPID = first_pid; m_NumPolygons = num_polygons; }
    };
    std::vector<editable_polygon_layer_type> m_addL;
};

//!\name Make methods
//@{
void Make_MeshSolidShape2_Box( EditableMeshSolidShape2& emss, const Vec2& half_sizes );
void Make_MeshSolidShape2_BoxWithHole( EditableMeshSolidShape2& emss, const Vec2& half_sizes, const Vec2& hole_half_sizes );
void Make_MeshSolidShape2_Box( EditableMeshSolidShape2& emss, const Vec2& half_sizes, unsigned int num_x, unsigned int num_y );
//@}

DCR_MeshSolidShape2* Create_DCR_MeshSolidShape2_From_PolygonalShape2( const MeshSolidShape2& mesh, const Transform2& mesh_tr, const Vec2* vec_mesh_sdof,
                                                                      const PolygonalShape2& polygonal, const Transform2& polygonal_tr, const Vec2* vec_polygonal_sdof );

#ifdef __GEO_MSS_ENABLE_BVH
/*! BV(Element)
  \note generic implementation for any supported BV
  \todo Add pre-transformed vec_node_pos version
  IMPORTANT: Both std::bind and boost::bind REQUIRE a copy-constructor
             for const MeshSolidShape2& mesh param, which causes
             *LOTS* of trouble, therefore we pass a const
             MeshSolidShape2* p_mesh to avoid them.
*/
template <typename EIT, typename BVT>
inline void GEBV_MeshSolidShape2_E( const MeshSolidShape2* p_mesh, const Transform2& tr, const MeshSolidShape2::sdof_type* vec_node_pos,
                                    EIT eid, BVT& bv )
{
    GEO_ASSERT( 0 != vec_node_pos );
    unsigned int it_he( p_mesh->P_FirstHEID(eid) );
    bv = BVT( tr * p_mesh->V_Pos( p_mesh->HE_OriginVID(it_he), vec_node_pos ) );
    do
    {
        it_he = p_mesh->HE_Next(it_he);
        bv.Merge( tr * p_mesh->V_Pos( p_mesh->HE_OriginVID(it_he), vec_node_pos ) );
    } while( it_he != p_mesh->P_FirstHEID(eid) );
    //bv.Extend( 2 );
};

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
    const DCR_MeshSolidShape2::Element& ed( pDCR->m_vecE[eid] );
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
    const DCR_MeshSolidShape2::Element& ed( pDCR->m_vecE[eid] );
    if( ed.m_BDOP[0].IsEmpty()
        || ed.m_BDOP[1].IsEmpty()
        || ed.m_BDOP[2].IsEmpty() ) { GEBV_MeshSolidShape2_E<EIT,BVT>( p_mesh, tr, vec_node_pos, eid, bv ); return; }
    // Call fast version if safe
    GEBV_MeshSolidShape2_BSlab<EIT,BVT>( p_mesh, pDCR, tr, vec_node_pos, eid, bv );
}

#endif //__GEO_MSS_ENABLE_BVH

} //namespace geo

#endif // GEO_SHAPE_MESHSOLIDSHAPE2_H
