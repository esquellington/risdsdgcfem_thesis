#ifndef GEO_SHAPE_TETSOLIDSHAPE3_H
#define GEO_SHAPE_TETSOLIDSHAPE3_H

#include <Geo/shape/IShape.h>
#include <Geo/bv/bv.h> // BV types for ComputeBV()
#include <Geo/bv/GBDOP.h>
#include <vector>

#define __GEO_TETSS_ENABLE_BVH
#ifdef __GEO_TETSS_ENABLE_BVH
#  include <Geo/bv/GBoundingVolumeHierarchy.h> //TEMPORAL: quick and dirty BVH
#endif

// #define __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_FIXED_LENGTH
#define __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH //should be enabled by default
#define __ENABLE_TETSS_DCR_PATCH_BDT //should be enabled by default
//#define __ENABLE_TETSS_DCR_PATCH_BDT_UNPACKED //TEMP: Adds 64b per Triangle, used to test correctness and speed of PACKED version, so far PACKED seems correct and EQUALLY FAST.

namespace geo
{

// We use specific TSS V/T counts (uint32) instead of
// geo::feature_index_type (uint16) to allow for LARGE triangle meshes
// (for viz)
typedef uint16 tetsolid3_feature_index_type;
typedef uint16 tetsolid3_feature_count_type;
static const uint16 cTetSolid3_InvalidFeatureIndex = 0xFFFF;

//\todo THIS IS NOT USED AT ALL but stored, RECONSIDER it and REMOVE
//if useless, will affect binary representation, but as we'll only
//REMOVE a field, it should be backwards compatible by just ignoring
//such field...
struct tetrahedron_vertex_topology
{
    uint8 m_Dummy; //unused, just testing cVersion = 2
};
struct tetrahedron_topology
{
    tetsolid3_feature_index_type m_vecVID[4];
    tetsolid3_feature_index_type m_vecNTID[4];
};
struct tetrahedron_boundary_face_topology
{
    tetsolid3_feature_index_type m_vecVID[3];
    //\todo tetsolid3_feature_index_type m_vecNBFID[3]; //topology may be necessary for CD, making this basically equivalent to TriSurfaceShape3
};

typedef std::pair<tetsolid3_feature_index_type,int> tetsolid3_face_id_type;

#ifdef __GEO_TETSS_ENABLE_BVH
//TEMPORAL: Consider generic IBVH so that each MSS can have a user-specified or best fitting one
//typedef bv::GBoundingVolumeHierarchy_ST_DG<bv::AABB3,tetsolid3_feature_index_type> BVH_TetSolidShape3;
// typedef bv::GBoundingVolumeHierarchy_ST_DG<bv::DOP3_K6,tetsolid3_feature_index_type> BVH_TetSolidShape3;
typedef bv::GBoundingVolumeHierarchy_ST_DG<bv::DOP3_K14,tetsolid3_feature_index_type> BVH_TetSolidShape3;
//typedef bv::GBoundingVolumeHierarchy_ST_DG<bv::Sphere3,tetsolid3_feature_index_type> BVH_TetSolidShape3;
#endif

} //namespace geo


namespace geo
{

//----------------------------------------------------------------
// DCR_TetSolidShape3 \todo Move elsewhere!
//----------------------------------------------------------------

class IDCR3
{
public:
    IDCR3() {}
    virtual ~IDCR3() {}
};

/*! Detailed Collision Representation for TetSolidShape3
  \see DCR_MeshSolidShape2 for an exhaustive comment
*/
class DCR_TetSolidShape3: public IDCR3
{
public:
    enum EConstants
    {
        cVersion = 4, //3) Added valid P.NPID with var-length sub-arrays __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
        //cVersion = 3, //3) Added valid P.NPID with fixed-length array __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_FIXED_LENGTH
        //cVersion = 2, //2) Added vector area, removed UNPACKED to reduce mem usage/cache impact
        //cVersion = 1, //1) Added m_EID and BDOP-Tree
        //cVersion = 0
    };
public:
    struct Element
    {
        // Sub-array of global patches in this element, unconnected
        uint32 m_FirstPID;
        uint32 m_NumPatches; //\todo Consider uint16 or even fitting first/num in 32b (eg: 24:8 or even 20:12)
        // Sub-array of global V in this element, to transform all at once
        uint32 m_FirstVID;
        uint32 m_NumVertices; //\todo Consider uint16 or even fitting first/num in 32b (eg: 24:8 or even 20:12)
        // Barycentric-DOP
        bv::BDOP4 m_BDOP; //\todo COULD be quantized [0,1] to uint8 (conservative min/max ==> floor/ceil)
        uint32 m_BDOP_BestSlabIdx; //2b would be enough!
        /* Barycentric-DOP
        Vec3 m_BDOP_Vertices[??];
        */
    };
    struct Patch
    {
        uint32 m_EID; //Useful to identify the DCR.E when navigating to neighbour patches, faster than searching new EID, even if topologically
        // Sub-array of global T in this patch
        uint32 m_FirstTID;
        uint32 m_NumTriangles; //\todo Consider uint16 or even fitting first/num in 32b (eg: 24:8 or even 20:12)
        //\todo CONSIDER m_FirstVID/m_NumVertices for patch sub-array
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_FIXED_LENGTH
        enum EConstants { cMaxNeighbours = 10 };
        uint32 m_NumNeighbours;
        uint32 m_vecNPID[cMaxNeighbours];
        inline Patch() //TEMP: Fill to test!
            {
                m_EID = 0xFFFFFF; m_FirstTID = 0xFFFFFFFF; m_NumTriangles = 0xFFFFFFFF;
                m_NumNeighbours = 0xFFFFFFFF;
                for( int i=0; i<cMaxNeighbours; i++ ) m_vecNPID[i] = 0xFFFFFFFF;
            }
#endif
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
        uint32 m_NumNeighbours;
        uint32 m_FirstIndexNPID;
        inline Patch() //TEMP: Fill to test!
            {
                m_EID = 0xFFFFFF; m_FirstTID = 0xFFFFFFFF; m_NumTriangles = 0xFFFFFFFF;
                m_NumNeighbours = 0xFFFFFFFF;
                m_FirstIndexNPID = 0xFFFFFFFF;
            }
#endif
#ifdef __ENABLE_TETSS_DCR_PATCH_BDT
        typedef std::pair<uint32,uint32> bdt_node_tip_range;
#endif
        //\todo Store per-patch aggregate data (ex: avgpos, avgnormal|integrated N*dA, total_area...)
        Vec3 m_VectorArea; //\sum N*dA
        //\todo Vec3 m_Centroid; float m_Area;
    };
    struct Triangle
    {
        inline Triangle() //TEMP: Fill it with invalid to check BST, remove when sufficiently tested!
            {
                for( int i=0; i<3; i++ ) m_vecVID[i] = 0xFFFFFFFF;
                for( int i=0; i<3; i++ ) m_vecNTID[i] = 0xFFFFFFFF;
#ifdef __ENABLE_TETSS_DCR_PATCH_BDT_UNPACKED
                for( int i=0; i<4; i++ ) for( int j=0; j<2; j++ ) m_vecBDOP[i][j] = 0xFF;
#endif
            }
#ifdef __ENABLE_TETSS_DCR_PATCH_BDT
        /*\todo
          - For large meshes (16M) Glboal VID/NTID could be stored in 24b, leaving 6x8b free.
          - We could also limit #tri/element to 2^16, thus saving DCR.E local VID in just 16b
          - Zero Memory BDOP-Tree:
            - 3x8b freed in m_vecNTID[3] can be used to save
              a) BSlab min/max [0..1] quantized to [0..255] in 2x8b
                (could be 2x7b with [0..127] quantization)
                IMPORTANT: Quantized min/max must use floor/ceil respectively
              b) #axis (2b) (could be implicit using alternate ordering...)
              c) Either 6b or 8b remain free in the upper 8b of the third m_vecNTID
            - If m_vecVID[3] uses 16b local VID, theres a 16b slot due
              to padding for the whole structure. This could be used
              to store 2x8b BSlab min/max.
              => THUS, if we use the 8b remaining in c) we could store
                 the min/max for all 4 BDOP axis on each BDOP-Tree
                 node
              - OR, we can use this 16b to store the R child offset
                and split the tree using SAH instead of a worse median
                cut.
              - NOTICE that we do not need to store #tri/node, this
                can be computed on the fly during a recursive
                traversal
              - The current BDOP can also be constructed recursively,
                starting from the top BDOP (explicit), updating a
                given axis range when entering a new node.
             - LEAF nodes are implicit if #tri/node = 1, but larger
               leaves are possible if we mark leaf nodes with 1b. This
               would allow SIMD 1 ray vs N tri
             - Using 22/21/20b per NTID, 4M/2M/1M tri are possible, and
               extra 6/9/12b are freed.
           - Diff with ZMT and NMH
             - They do NOT store bounds explicitly but implicitly in
               representative triangles, that need to cover the
               EXTREMES of the node BV and need to be checked
               unconditionally. therefore, are NOT CULLED as
               efficiently as in a standard BVH. In ZMT with numtri >>
               2 this seems particularly inefficient.
             - The "Box Cache" idea of ZMT sounds really really bad...
           - The BDOP-Tree should be also used for intersection queries.

           - IMPORTANT: The ray should be quantized ONCE and its
             q-BDOP intervals tested against successive node q-BDOP
             - \todo Should Ray q-BDOP be clipped when entering a node??
        */

        /*
          Overall, saving the whole BDOP in each BDTNode seems the simplest and most efficient solution
          - BSTN are self-contained
          - Uniform data representation
          - RayVsBDOP only tested once per node, clipped and propagated to children
          Layout:
          - VID and NTID have 21b (2M values)
          - 3x8b upper byte in each VID store BDOPq[0].Min/Max, BDOP[1].Min
            - Extract with mask 0xFF000000 and >> 24
          - 3x3b after upper byte in each VID store BDOP[1].Max (3x3b=9b, we discard 1b)
            - Extract with mask 0x00E00000 and >> 21 for 3b and 0x00C00000 and >> 22 for 2b
          - 3x8b upper byte in each NTID store BDOPq[2].Min/Max, BDOP[3].Min
          - 3x3b after upper byte in each VID store BDOP[3].Max (3x3b=9b, we discard 1b)
        */
        enum EConstants { eBDT_IntervalQuantizationBits = 8,
                          eBDT_MaskVID = 0x001FFFFF,
                          eBDT_MaskNTID = 0x001FFFFF,
                          cInvalidTID = eBDT_MaskNTID };
        inline uint32 GetVID( int i ) const { return m_vecVID[i] & eBDT_MaskVID; }
        inline uint32 GetNTID( int i ) const { return m_vecNTID[i] & eBDT_MaskNTID; }
        inline void SetVID( int i, uint32 vid ) { GEO_ASSERT(0 == (vid & ~eBDT_MaskVID)); m_vecVID[i] = (m_vecVID[i] & ~eBDT_MaskVID) | (vid & eBDT_MaskVID); }
        inline void SetNTID( int i, uint32 ntid ) { GEO_ASSERT(0 == (ntid & ~eBDT_MaskNTID)); m_vecNTID[i] = (m_vecNTID[i] & ~eBDT_MaskNTID) | (ntid & eBDT_MaskNTID); }
        //\name BarycentricDOPTreeNode
        //@{
        /* Quantized BDOP4 with 8b quantization but int32 interval
           values, as GInterval REQUIRES signed T
        */
        class BDOP4q
        {
        public:
            typedef mal::GInterval<int32> interval_q_type;
            inline BDOP4q() {}
            // Init from BDOP
            inline BDOP4q( const bv::BDOP4& bdop4 )
                {
                    for( int i=0; i<4; i++ )
                        m_vecInterval[i].SetMinMax( mal::Clamp<int32>( mal::Floor(mal::Clamp01(bdop4[i].Min())*0x000000FF), 0, 0x000000FF ),
                                                    mal::Clamp<int32>( mal::Ceil(mal::Clamp01(bdop4[i].Max())*0x000000FF), 0, 0x000000FF ) );
                }
            // Cast to BDOP
            inline operator bv::BDOP4() const
                {
                    bv::BDOP4 bdop4;
                    const Real q2f( mal::Rcp<Real>( (1<<eBDT_IntervalQuantizationBits)-1 ) );
                    for( int i=0; i<4; i++ )
                        bdop4[i].SetMinMax( q2f*m_vecInterval[i].Min(), q2f*m_vecInterval[i].Max() );
                    return bdop4;
                }
            const interval_q_type& operator[](int i) const { return m_vecInterval[i]; }
            interval_q_type& operator[](int i) { return m_vecInterval[i]; }
        private:
            interval_q_type m_vecInterval[4];
        };

        inline BDOP4q BDTN_BDOPq() const
            {
                BDOP4q bdop4q;
                bdop4q[0].SetMinMax( (m_vecVID[0] & 0xFF000000) >> 24, (m_vecVID[1] & 0xFF000000) >> 24 ); //upper 8b
                bdop4q[1].SetMinMax( (m_vecVID[2] & 0xFF000000) >> 24, //upper 8b
                                     ((m_vecVID[0] & 0x00E00000) >> (21-5)) //3b [21..24]
                                     | ((m_vecVID[1] & 0x00E00000) >> (21-2)) //3b [21..24]
                                     | ((m_vecVID[2] & 0x00C00000) >> 22) ); //2b [21..23]
                bdop4q[2].SetMinMax( (m_vecNTID[0] & 0xFF000000) >> 24, (m_vecNTID[1] & 0xFF000000) >> 24 );
                bdop4q[3].SetMinMax( (m_vecNTID[2] & 0xFF000000) >> 24, //upper 8b
                                     ((m_vecNTID[0] & 0x00E00000) >> (21-5)) //3b [21..24]
                                     | ((m_vecNTID[1] & 0x00E00000) >> (21-2)) //3b [21..24]
                                     | ((m_vecNTID[2] & 0x00C00000) >> 22) ); //2b [21..23]
                return bdop4q;
            }

        inline void BDTN_Init( const bv::BDOP4& bdop4 )
            {
                BDOP4q bdop4q( bdop4 );
                // Save bdop4q into the right bits...
                // BDOP4q axis 0,2 min/max go into upper 8b and axis 1,3 min too
                m_vecVID[0] = (m_vecVID[0] & eBDT_MaskVID) | ((bdop4q[0].Min() & 0x000000FF) << 24);
                m_vecVID[1] = (m_vecVID[1] & eBDT_MaskVID) | ((bdop4q[0].Max() & 0x000000FF) << 24);
                m_vecVID[2] = (m_vecVID[2] & eBDT_MaskVID) | ((bdop4q[1].Min() & 0x000000FF) << 24);
                m_vecNTID[0] = (m_vecNTID[0] & eBDT_MaskNTID) | ((bdop4q[2].Min() & 0x000000FF) << 24);
                m_vecNTID[1] = (m_vecNTID[1] & eBDT_MaskNTID) | ((bdop4q[2].Max() & 0x000000FF) << 24);
                m_vecNTID[2] = (m_vecNTID[2] & eBDT_MaskNTID) | ((bdop4q[3].Min() & 0x000000FF) << 24);
                // BDOP4q axis 1,3 max are in 3b+3b+2b packets at the [21..24) range spread accross vecVID[3] and vecNTID[3]
                m_vecVID[0] = (m_vecVID[0] & 0xFF1FFFFF) | ((bdop4q[1].Max() & 0x000000E0) << 16); //3b upper in 8b
                m_vecVID[1] = (m_vecVID[1] & 0xFF1FFFFF) | ((bdop4q[1].Max() & 0x0000001C) << 19); //3b mid in 8b
                m_vecVID[2] = (m_vecVID[2] & 0xFF3FFFFF) | ((bdop4q[1].Max() & 0x00000003) << 22); //2b lower in 8b
                m_vecNTID[0] = (m_vecNTID[0] & 0xFF1FFFFF) | ((bdop4q[3].Max() & 0x000000E0) << 16);
                m_vecNTID[1] = (m_vecNTID[1] & 0xFF1FFFFF) | ((bdop4q[3].Max() & 0x0000001C) << 19);
                m_vecNTID[2] = (m_vecNTID[2] & 0xFF3FFFFF) | ((bdop4q[3].Max() & 0x00000003) << 22);
#ifdef __ENABLE_TETSS_DCR_PATCH_BDT_UNPACKED
                for( int i=0; i<4; i++ )
                {
                    m_vecBDOP[i][0] = bdop4q[i].Min();
                    m_vecBDOP[i][1] = bdop4q[i].Max();
                }
#endif
            }
#ifdef __ENABLE_TETSS_DCR_PATCH_BDT_UNPACKED
        inline BDOP4q BDTN_BDOPq_UNPACKED() const
            {
                BDOP4q bdop4q;
                for( int i=0; i<4; i++ ) bdop4q[i].SetMinMax( m_vecBDOP[i][0], m_vecBDOP[i][1] );
                return bdop4q;
            }
#  endif
        //@}
#else //__ENABLE_TETSS_DCR_PATCH_BDT
        enum EConstants { cInvalidTID = 0xFFFFFFFF };
        inline uint32 GetVID( int i ) const { return m_vecVID[i]; }
        inline uint32 GetNTID( int i ) const { return m_vecNTID[i]; }
        inline void SetVID( int i, uint32 vid ) { m_vecVID[i] = vid; }
        inline void SetNTID( int i, uint32 ntid ) { m_vecNTID[i] = ntid; }
#endif
    private:
        uint32 m_vecVID[3]; //\todo COULD be uint16 if indices local to Element m_FirstVID
        uint32 m_vecNTID[3]; //\todo COULD be uint16 if indices local to Element m_FirstTID and cross-patch topology stored separately
#  ifdef __ENABLE_TETSS_DCR_PATCH_BDT_UNPACKED
        uint8 m_vecBDOP[4][2];
#  endif
    };

    uint32 m_NumElements;
    uint32 m_NumPatches;
    uint32 m_NumTriangles;
    uint32 m_NumVertices;
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
    uint32 m_NumNPID;
#endif
    Element* m_vecE; //Indices parallel to TetSS elements
    Patch* m_vecP; //stores all patches in per-element contiguous subarrays
    Triangle* m_vecT; //
    Vec3* m_vecV; // Relative to TetSS tr/dof
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
    uint32* m_vecNPID;
#endif
    //\todo Store as uint8* m_pBuffer if shared is supported

    inline DCR_TetSolidShape3()
    : m_NumElements(0), m_NumPatches(0), m_NumTriangles(0), m_NumVertices(0)
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
    , m_NumNPID(0)
#endif
    , m_vecE(0), m_vecP(0), m_vecT(0), m_vecV(0)
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
    , m_vecNPID(0)
#endif
        {}
    inline ~DCR_TetSolidShape3()
        {
            if(m_vecE) delete[] m_vecE; m_vecE = 0;
            if(m_vecP) delete[] m_vecP; m_vecP = 0;
            if(m_vecT) delete[] m_vecT; m_vecT = 0;
            if(m_vecV) delete[] m_vecV; m_vecV = 0;
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
            if(m_vecNPID) delete[] m_vecNPID; m_vecNPID = 0;
#endif
        }

    inline void SetBakedData( bool b_shared,
                              uint32 num_e, uint32 num_p, uint32 num_t, uint32 num_v, uint32 num_npid,
                              const Element* vec_e, const Patch* vec_p, const Triangle* vec_t, const Vec3* vec_v, const uint32* vec_npid )
        {
            if( b_shared ) GEO_LOG_WARNING("DCR_TetSolidShape3::SetBakedData ignores b_shared by now...");
            m_NumElements = num_e;
            m_NumPatches = num_p;
            m_NumTriangles = num_t;
            m_NumVertices = num_v;
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
            m_NumNPID = num_npid;
#endif
            m_vecE = new Element[num_e];
            m_vecP = new Patch[num_p];
            m_vecT = new Triangle[num_t];
            m_vecV = new Vec3[num_v];
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
            m_vecNPID = new uint32[num_npid];
#endif
            memcpy( m_vecE, vec_e, sizeof(Element)*num_e );
            memcpy( m_vecP, vec_p, sizeof(Patch)*num_p );
            memcpy( m_vecT, vec_t, sizeof(Triangle)*num_t );
            memcpy( m_vecV, vec_v, sizeof(Vec3)*num_v );
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
            memcpy( m_vecNPID, vec_npid, sizeof(uint32)*num_npid );
#endif
        }

        // Find neighbour-edge-in-triangle
        inline uint32 Find_NEINT( uint32 tid, uint32 ntid ) const
        {
            for( uint32 it_neint=0; it_neint<3; it_neint++ )
                if( m_vecT[ntid].GetNTID(it_neint) == tid )
                    return it_neint;
            return 0xFFFFFFFF;
        }

        inline bool Is_TID_In_PID( uint32 tid, uint32 pid ) const
        {
            return tid >= m_vecP[pid].m_FirstTID
                && tid < m_vecP[pid].m_FirstTID + m_vecP[pid].m_NumTriangles;
        }

        /* Find PID that contains a TID, globally:
           IMPORTANT: EXTREMELY SLOW
        */
        inline uint32 Find_PID_From_TID( uint32 tid ) const
        {
            for( uint32 it_pid=0; it_pid<m_NumPatches; it_pid++ )
                if( Is_TID_In_PID( tid, it_pid ) )
                    return it_pid;
            return 0xFFFFFFFF;
        }

        /* This should be fast enough For FloodT/FloodP, where a
           neighbour tid is always in the same P or in one of its
           neighbours.
           \note Works for ANY tid, regardless of hint_pid
        */
        inline uint32 Find_PID_From_TID_Hint( uint32 tid, uint32 hint_pid ) const
        {
            if( Is_TID_In_PID( tid, hint_pid ) )
                return hint_pid;
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
            const Patch& pd( m_vecP[hint_pid] );
            for( uint32 it_npip=0; it_npip<pd.m_NumNeighbours; it_npip++ )
            {
                uint32 npid( m_vecNPID[pd.m_FirstIndexNPID+it_npip] );
                if( Is_TID_In_PID( tid, npid ) )
                    return npid;
            }
#endif
            // Just in case tid is not in hint_pid or one of its neighbours...
            for( uint32 it_pid=0; it_pid<m_NumPatches; it_pid++ )
                if( Is_TID_In_PID( tid, it_pid ) )
                    return it_pid;
            return 0xFFFFFFFF;
        }
};

//----------------------------------------------------------------
// TetSolidShape3
//----------------------------------------------------------------

/* Tetrahedron Solid 3D
   Stores:
   - Vertices
   - Tetrahedrons
   - Boundary Triangles
*/
class TetSolidShape3: public IShape3
{
public:
    enum EConstants
    {
        cVersion = 3 //Added layers
        //cVersion = 2 //Removed unused vertex attribs
        //cVersion = 1 //Initial version
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
    struct vertex_type: public tetrahedron_vertex_topology {};
    struct tetrahedron_type: public tetrahedron_topology {};
    struct boundary_face_type: public tetrahedron_boundary_face_topology {};
    struct tetrahedron_layer_type
    {
        tetsolid3_feature_index_type m_FirstTID;
        tetsolid3_feature_count_type m_NumTetrahedrons;
    };

public:
    TetSolidShape3();
    ~TetSolidShape3();

    //!\name Dynamic IShape implementation
    //@{
    EShapeType GetType() const { return eShape_TetSolid3; }
    unsigned int GetNumSDOF() const { return m_NumV; }
    const sdof_type* GetVecDefaultSDOF() const { return m_vecPoints; }
    const Real* GetVecDefaultDOF() const { return reinterpret_cast<const Real*>(m_vecPoints); }
    void ComputeBVD( bv::BoundingVolume3& bv, const transform_type& transform, const sdof_type* vec_sdof ) const;

    IDomainSampler3* CreateDomainSampler() const { return 0; } //\todo
    //@}

    //! Init from external arrays
    void SetBakedData( bool b_shared,
                       uint32 num_v, uint32 num_t, uint32 num_bf, uint32 num_ll,
                       const Vec3* vec_points, const vertex_type* vec_v, const tetrahedron_type* vec_t, const boundary_face_type* vec_bf, const tetrahedron_layer_type* vec_l );
    //TEMPORAL: hacks to add annotations to non-editable TetSS
    inline void SetBakedDCR_StrictlyNonshared_UglyHack( DCR_TetSolidShape3* p_nonshared_dcr ) { m_pDCR = p_nonshared_dcr; } //never shared, breaks SetBakedData() strict const enforcement
#ifdef __GEO_TETSS_ENABLE_BVH
    inline void SetBakedBVH_StrictlyNonshared_UglyHack( BVH_TetSolidShape3* p_nonshared_bvh ) { m_pBVH = p_nonshared_bvh; } //never shared, breaks SetBakedData() strict const enforcement
#endif

    inline uint32 GetNumV() const { return m_NumV; }
    inline uint32 GetNumT() const { return m_NumT; }
    inline uint32 GetNumBF() const { return m_NumBF; }
    inline uint32 GetNumL() const { return m_NumL; }

    inline const Vec3* GetVecPoints() const { return m_vecPoints; }
    inline const vertex_type* GetVecV() const { return m_vecV; }
    inline const tetrahedron_type* GetVecT() const { return m_vecT; }
    inline const boundary_face_type* GetVecBF() const { return m_vecBF; }
    inline const tetrahedron_layer_type* GetVecL() const { return m_vecL; }

    // \todo
    inline bool IsClosed() const { return true; } //!< Must be closed by construction
    inline bool IsSimplyConnected() const { return true; } //!< Check! false if there are holes
    inline bool IsConnected() const { return true; } //!< Check!
    inline bool IsConvex() const { return false; } //!< Not so useful for solid meshes...

    // Vertex access
    inline Vec3 V_Pos_0( uint32 vid ) const { return m_vecPoints[vid]; }
    inline Vec3 V_Pos( uint32 vid, const sdof_type* vec_sdof ) const { return vec_sdof[vid]; }

    // Tetrahedron access
    inline uint32 T_VID( uint32 tid, uint32 vit ) const { return m_vecT[tid].m_vecVID[vit]; }
    inline void T_VecVID( uint32 tid, uint32* vec_vid ) const { for( int i=0;i<4;i++) vec_vid[i] = m_vecT[tid].m_vecVID[i]; }
    inline uint32 T_NTID( uint32 tid, uint32 ntit ) const { return m_vecT[tid].m_vecNTID[ntit]; }
    Vec3 T_Barycenter_0( uint32 tid ) const;
    Vec3 T_Barycenter( uint32 tid, const sdof_type* vec_sdof ) const;
    void T_Edge_0( uint32 tid, uint32 eit, Vec3& p0, Vec3& p1 ) const;

    // Boundary Face access
    inline uint32 BF_VID( uint32 bfid, uint32 vibf ) const { return m_vecBF[bfid].m_vecVID[vibf]; }
    Vec3 BF_Barycenter( uint32 bfid, const sdof_type* vec_sdof ) const;

    // Layers access
    inline uint32 L_FirstTID( unsigned int lid ) const { return m_vecL[lid].m_FirstTID; }
    inline uint32 L_NumT( unsigned int lid ) const { return m_vecL[lid].m_NumTetrahedrons; }

    // Get annotations or 0 if none
    /*\todo Annotations are built for specific shape params and
     default DOF, which they can refer implicitly assuming constant
     indices/values, if annotation-specific data is required to deal
     with runtime DOF, this annotation-DOF should be allocated by
     IObject. Ex: Transformed/Deformable vertex positions in a DCR
    */
    inline const DCR_TetSolidShape3* GetDCR() const { return m_pDCR; }
#ifdef __GEO_TETSS_ENABLE_BVH
    inline BVH_TetSolidShape3* GetBVH() const { return m_pBVH; } //\note non-const until topology/geometry are split, because we need to refit it!
#endif

protected:
    void ClearBakedData();

protected:

    uint32 m_NumV;
    uint32 m_NumT;
    uint32 m_NumBF;
    uint32 m_NumL;

    const Vec3* m_vecPoints;
    const vertex_type* m_vecV;
    const tetrahedron_type* m_vecT;
    const boundary_face_type* m_vecBF;
    const tetrahedron_layer_type* m_vecL;

    uint32* m_pBuffer; //!< If exists, the shape is NOT SHARED

    DCR_TetSolidShape3* m_pDCR; //strictly nonshared, by now...
#ifdef __GEO_TETSS_ENABLE_BVH
    BVH_TetSolidShape3* m_pBVH; //strictly nonshared, by now...
#endif

    /*\todo
private:
    class DomainSampler;
    */
};

/* Editable Tetrahedron Solid Shape 3D
*/
class EditableTetSolidShape3: public TetSolidShape3
{
public:
    EditableTetSolidShape3();
    ~EditableTetSolidShape3();

    void Clear();
    void Set( const TetSolidShape3& tss3 );

    //\name Begin/End edition
    //@{
    void BeginEdition(); //\note Incremental, preserves BakedData if any, call Clear() to reset completely
    bool EndEdition();
    inline bool IsBeingEdited() const { return m_IsBeingEdited; }

    tetsolid3_feature_index_type AddVertex( const Vec3& point );
    tetsolid3_feature_index_type AddTetrahedron( uint32 vid0, uint32 vid1, uint32 vid2, uint32 vid3 );
    void SetVertex( uint32 vid, const Vec3& point );
    void Transform( const Transform3& tr );
    //@}

    //\name Global edition (outside Begin/End)
    //@{
    bool RemoveTetrahedrons( uint32 num_tid, const uint32* vec_tid );
    bool FixNonManifoldFeatures( const std::vector<tetsolid3_face_id_type>& vec_split_fid );
    //@}

    // Add annotations
    bool AddDCR( const IShape3* p_shape, const Transform3& tr_s2tss );//, EDCRtype dcrt ); //TriSurface3, eDCR_Partition );
#ifdef __GEO_TETSS_ENABLE_BVH
    bool AddBVH(); //bvh type!?
#endif

private:

    bool RebuildAllFromEditData();
    void ClearEditData();
    void RebuildBakedData();
    bool FixDegeneracies();
    void RebuildTopology();
    void RebuildLayers();

    bool FixDegenerateVertices();
    bool FixDegenerateTetrahedrons();

private:
    struct editable_vertex_type: public vertex_type
    {
        inline editable_vertex_type() : m_Pos(0,0,0) {} //\note Required by m_addV.resize()
        inline editable_vertex_type( const Vec3& pos )
        : m_Pos(pos) {}
        Vec3 m_Pos;
    };
    struct editable_tetrahedron_type: public tetrahedron_type
    {
        inline editable_tetrahedron_type( tetsolid3_feature_index_type vid0, tetsolid3_feature_index_type vid1, tetsolid3_feature_index_type vid2, tetsolid3_feature_index_type vid3 )
            {
                m_vecVID[0] = vid0; m_vecVID[1] = vid1; m_vecVID[2] = vid2; m_vecVID[3] = vid3;
                m_vecNTID[0] = cTetSolid3_InvalidFeatureIndex; m_vecNTID[1] = cTetSolid3_InvalidFeatureIndex; m_vecNTID[2] = cTetSolid3_InvalidFeatureIndex; m_vecNTID[3] = cTetSolid3_InvalidFeatureIndex;
            }
    };
    struct editable_boundary_face_type: public boundary_face_type
    {
        inline editable_boundary_face_type( tetsolid3_feature_index_type vid0, tetsolid3_feature_index_type vid1, tetsolid3_feature_index_type vid2 )
            { m_vecVID[0] = vid0; m_vecVID[1] = vid1; m_vecVID[2] = vid2; }
    };
    struct editable_tetrahedron_layer_type: public tetrahedron_layer_type
    {
        editable_tetrahedron_layer_type()
            { m_FirstTID = 0; m_NumTetrahedrons = 0; } //\note Required by m_addL.resize()
        editable_tetrahedron_layer_type( tetsolid3_feature_index_type first_tid, tetsolid3_feature_count_type num_tetrahedrons )
            { m_FirstTID = first_tid; m_NumTetrahedrons = num_tetrahedrons; }
    };

    bool m_IsBeingEdited;
    std::vector<editable_vertex_type> m_addV;
    std::vector<editable_tetrahedron_type> m_addT;
    std::vector<editable_boundary_face_type> m_addBF;
    std::vector<editable_tetrahedron_layer_type> m_addL;
};

//\name Free functions
//@{
Real ComputeVolume( const TetSolidShape3& tss3, const Transform3& transform, const Vec3* vec_sdof );
//@}

//!\name Make methods
//@{
void Make_TetSolidShape3_Box( EditableTetSolidShape3& etss3, const Vec3& half_sizes );
void Make_TetSolidShape3_Box( EditableTetSolidShape3& etss3, const Vec3& half_sizes, unsigned int num_x, unsigned int num_y, unsigned int num_z );
//@}

DCR_TetSolidShape3* Create_DCR_TetSolidShape3_From_TriSurfaceShape3( const TetSolidShape3& solid, const Transform3& solid_tr, const Vec3* vec_solid_sdof,
                                                                     const TriSurfaceShape3& surface, const Transform3& surface_tr, const Vec3* vec_surface_sdof );

bool Check_DCR_TetSolidShape3( const DCR_TetSolidShape3& dcr );

#ifdef __GEO_TETSS_ENABLE_BVH
/*! BV(Element)
  \note generic implementation for any supported BV
  \todo Add pre-transformed vec_node_pos version
  IMPORTANT: Both std::bind and boost::bind REQUIRE a copy-constructor
             for const TriSurfaceShape3& surface param, which causes
             *LOTS* of trouble, therefore we pass a const
             TriSurfaceShape3* p_shape to avoid them.
*/
template <typename EIT, typename BVT>
inline void GEBV_TetSolidShape3_E( const TetSolidShape3* p_solid, const TetSolidShape3::transform_type& tr, const TetSolidShape3::sdof_type* vec_node_pos,
                                   EIT eid, BVT& bv )
{
    GEO_ASSERT( 0 != vec_node_pos );

    bv = BVT( tr * p_solid->V_Pos( p_solid->T_VID(eid,0), vec_node_pos ) );
    bv.Merge( tr * p_solid->V_Pos( p_solid->T_VID(eid,1), vec_node_pos ) );
    bv.Merge( tr * p_solid->V_Pos( p_solid->T_VID(eid,2), vec_node_pos ) );
    bv.Merge( tr * p_solid->V_Pos( p_solid->T_VID(eid,3), vec_node_pos ) );
};

/*! BV(BDOP)
  \note Has AABB3, DOP3_K6, DOP3_K14 specializations, and falls back to BV(E) for Sphere3
  \todo Add pre-transformed vec_node_pos version
  \pre 0 <= eid << p_dcr->m_NumElements
*/
template <typename EIT, typename BVT>
void GEBV_TetSolidShape3_BDOP( const TetSolidShape3* p_solid, const DCR_TetSolidShape3* p_dcr, const Transform3& tr, const TetSolidShape3::sdof_type* vec_node_pos,
                               EIT eid, BVT& bv );
#endif

} //namespace geo

#endif // GEO_SHAPE_TETSOLIDSHAPE3_H
