#ifndef GEO_BV_GBOUNDING_VOLUME_HIERARCHY_H
#define GEO_BV_GBOUNDING_VOLUME_HIERARCHY_H

#include <Geo/Config.h>
#include <Geo/bv/BoundingVolume.h>
#include <vector>
#include <boost/function.hpp>

//For specific instantiations
#include <Geo/bv/GAABB.h>
#include <Geo/bv/GKDOP.h>
#include <Geo/bv/GSphere.h>

namespace geo {
namespace bv {

/*! Simple implementation of a generic BVH
  - Entries are not stored, only their indices
  - Leaf-nodes store per-entry BV
  - Internal-nodes store the BV that bounds their children BVs

  Assumes:
  - bv_type.Init( entry_type )
  Requires:
  - bv_type.Merge( bv_type )
  - bv_type.TestOverlap( bv_type )

  Design decisions:
  - We assume a reasonable number of entries O(10^3) (eg: Elements in
    a FEM mesh)
  - We DO NOT store entry geometry, the GBVH only needs the BV of
    externally allocated entries through the GEBVD delegate.
  - We DO NOT try to support lazy refit or track any kind of per-entry
    modification status/timestamp, assume ALL entries have changed
    on Refit().
    - Externally, global Refit() can be called lazily, though.
  - We assume a SINGLE global Refit() call per frame followed by
    SEVERAL lazy Test() calls.
  - We use recursive calls for some ops, instead of faster iterative
    stack-based algorithms to keep code SIMPLE

  \todo Use Scratchpad-allocated stacks OR, reserve a member
        m_stackNID with known m_Height and reuse if for all hierarchy
        traversal methods.

  \todo Wrapped Vs Layered: Layered is easier to program and requires
        NO entry classification, but is much less tight. 
        method?) and layered internal nodes (bottom-up, O(#nodes)

  \todo For EntryT = Segment/Triangle/Tetrahedron, non-leaf BVH nodes
        do NOT NEED to eval the actual entries, only the actual vertices
        that compose them...

  \todo Fully _Static and _Dynamic policies are also possible, but by
        now we'll only support _StaticTopology_DynamicGeometry
        (_ST_DG) as required by deformable solids.

  \todo Geometry could be in specific reference system.
        - BVH1.Test(BVH2) will need to transform BV1 into BV2 reference system.
          - Specialize for Translation-only, and ignore Rotation, so that
            AABB/KDOP remain so, otherwise OBB/OKDOP would result,
            complicating TestOverlap(BV1,BV2)
        => Provide a TestBoundingVolumePairDelegate (TBVPD) that
           receives both BV1 and BV2, transforms/converts them at
           will, tests for ovrlap and returns true/false. The BVH then
           would be decoupled from the actual BV reference systems.
           => Default to bv1.TestOverlap(bv2)
           => Better, call it TEPD TestEntryPairDelegate and receive
              entry indexes AND their BV, so that user-specific per-entry
              filtering is possible if, for example, BDOP is available...
              => TEPD does not NEED to work for non-leaf nodes

   \todo Could accept different refits:
   - Layered: Refit-Merge-BV: using BV.Merge(BV), so that nodes store BV(left_node,right_node) => CHEAPER
     - \todo Per-entry BV do NOT NEED to be stored anywhere, they can
       be used to refit their parent node BV incrementally and
       discarded afterwards.
   - Wrapped: Refit-Merge-Entry: using BV.Merge(entry), so that nodes store BV(left_entries, right_entries) => TIGHTER
     => Could be done bottom-up propagating current leaf-list and computing BV-entries.
   => Per AABB/KDOP, Layered == Wrapped!?!?!?! (per Spheres no, paper
      Guibas "2002 Collision Detection for Deforming Necklaces" ho
      explica amb detall.

   \todo Forcing leaf nodes to contain exactly 1 entry simplifies the
         code considerably, but also increases tree height and requires
         storing per-entry stuff. However, NOT doing it also complicates
         data structures (non-leaf nodes need to store the entry range in a
         specific entry index array).
         => Given a "full tree" with leaf=entry, a "partial tree" can
            be built by condensing whole subtrees into leaves.

   \todo The GBoundingVolumeHierarchy_ST_DG structure could be split
         into an editable ST part (serializable, shareable, does not
         depend on runtime DOF) and a non-editable DG part (exclusive,
         runtime-instance/DOF dependent)

   \todo boost::function Vs std::function Vs template param functor,
         which *COULD* be a lambda
         => boost::function is faster than std::function, but a
         lambda/functor template param would be WAY WAY faster,
*/
template <typename BoundingVolumeT, typename EntryIndexT>
class GBoundingVolumeHierarchy_ST_DG //_Static | _Dynamic
{
public:
    typedef BoundingVolumeT bv_type;
    typedef EntryIndexT entry_index_type;
    typedef entry_index_type entry_count_type;
    typedef entry_index_type node_index_type; //\todo Could be smaller if leaf >> entry, but by now leaf == entry, and num_nodes = 2*num_entries - 1
    typedef GBoundingVolumeHierarchy_ST_DG<BoundingVolumeT,EntryIndexT> bvh_type;
    typedef mal::GVec<Real,bv_type::cDimension> vec_type;
    enum EConstants { cDimension = bv_type::cDimension, cInvalidNodeId = node_index_type(~0), cMaxNodes = node_index_type(~0) - 1 };

    struct Node;

    typedef boost::function<void (entry_index_type eid, bv_type& bv)> GEBVD; //GetEntryBoundingVolumeDelegate
    typedef boost::function<bool (const bvh_type& bvh, const Node& node, uint32 level)> MND; //MapNodesDelegate
    typedef boost::function<bool (entry_index_type eid1, entry_index_type eid2, const bv_type& bv1, const bv_type& bv2)> TEPD; //TestEntryPairDelegate
    //typedef boost::function<bool (const bv_type& bv1, const bv_type& bv2)> TBVPD; //TestBoundingVolumePairDelegate

public:
    inline GBoundingVolumeHierarchy_ST_DG() : m_NumEntries(0), m_Height(0), m_GeometryTimeStamp(0) {}

    // Editable/DynamicTopology API
    //@{
    // Rebuild methods: Selects best available implementation
    void Rebuild( entry_count_type num_entries, const GEBVD& gebvd ) { Rebuild_BottomUp(num_entries,gebvd); } //\todo Choose BU/TD depending on num_entries!?!?
    void Rebuild_BottomUp( entry_count_type num_entries, const GEBVD& gebvd ) { Rebuild_BottomUp_N3(num_entries,gebvd); }
    void Rebuild_TopDown( entry_count_type num_entries, const GEBVD& gebvd ) { Rebuild_TopDown_NLogN(num_entries,gebvd); }
    //@}

    // Runtime/DynamicGeometry API
    //@{
    void Refit( const GEBVD& gebvd )
        {
            if( IsEmpty() ) return;
            // Update leaf BV from their single-entries (we assume ALL entries have changed)
            // \note weak-post-order node layout allows all leaves to be stored together (first sub-array or even separate array)
            // \note IMPORTANT: leaves/entries are NOT NECESSARILY in the same order (Rebuild_TopDown changes it!)
            for( unsigned int it_leaf=0; it_leaf<m_NumEntries; it_leaf++ ) //#leaves == #entries
                gebvd( m_vecNode[it_leaf].m_Topology.m_FirstEID, m_vecNode[it_leaf].m_Geometry.m_BV );
            // Bottom-Up refit
            // \note weak-post-order node layout ensures all children are up-to-date before their parent is processed
            for( unsigned int it_node=m_NumEntries; it_node<m_vecNode.size(); it_node++ )
            {
                Node& node( m_vecNode[it_node] );
                node.m_Geometry.m_BV = bv_type( m_vecNode[node.m_Topology.m_Left].m_Geometry.m_BV ).Merge( m_vecNode[node.m_Topology.m_Right].m_Geometry.m_BV );
            }
        }

    /*\todo
    void Refit_Lazy( const GEBVD& gebvd )
        {
            if( IsEmpty() ) return;
            // Update leaf BV from their single-entries (we assume ALL entries have changed)
            // \note weak-post-order node layout allows all leaves to be stored together (first sub-array or even separate array)
            // \note IMPORTANT: leaves/entries are NOT NECESSARILY in the same order (Rebuild_TopDown changes it!)
            for( unsigned int it_entry=0; it_entry<m_NumEntries; it_entry++ )
                gebvd( m_vecNode[it_leaf].m_Topology.m_FirstEID, m_vecNode[it_entry].m_Geometry.m_BV );
            //\todo stack-based refit where upwards propagation stops when already contained in previous BV
            //std::vector<uint32> stackNID; //\todo USE Scratchpad-Stack to avoid new/delete!
        }
    */
    //@}

    /*\todo Consider bv1_type != bv2_type, just requires adding template params.
      template<BoundingVolumeT2>
      bool Test( BoundingVolumeT2& bv,
                 boost::function<bool (const bv_type& bv1, const BoundingVolumeT2& bv2)> tbvpd,
                 std::vector<entry_index_type>& vec_entry_id ) const
    */
    bool Test( const bv_type& bv2, std::vector<entry_index_type>& vec_entry_id ) const
        {
            vec_entry_id.clear();
            if( IsEmpty() ) return false;
            m_stackNID.clear(); //\todo Consider Scratchpad-Stack to avoid new/delete!
            m_stackNID.push_back( GetRootIndex() );
            while( !m_stackNID.empty() )
            {
                const Node& node( m_vecNode[m_stackNID.back()] );
                m_stackNID.pop_back();
                if( TestOverlap( node.m_Geometry.m_BV, bv2 ) ) //\todo if( tbvpd( node.m_Geometry.m_BV, bv2 ) )
                {
                    if( node.IsLeaf() )
                    {
                        //\todo If ever 1 leaf > 1 entry, add all leaf entries to vec_entry_id
                        vec_entry_id.push_back( node.m_Topology.m_FirstEID );
                    }
                    else
                    {
                        m_stackNID.push_back( node.m_Topology.m_Left );
                        m_stackNID.push_back( node.m_Topology.m_Right );
                    }
                }
            }
            return vec_entry_id.size() > 0;
        }

    /*\todo By now, simple recursive call on bvh root nodes.
      \todo DECIDE IF TBVPD or TEPD is used!!!!
      \todo TBVPD parameter could support:
            - bv1_type != bv2_type
              - requires adding additional template params
              - entry_pair_container could also support different index types
            - bv1 reference system != bv2 reference system
              - TBVPD/TEPD could implicitly convert before testing
      \todo TEPD could support:
            - User-specific testing for entry1 vs entry2 (eg: exact
              B-Slab test after detecting BV(B-Slab) overlap)
      \todo Consider returning leaf pairs, instead of entry pairs?
            => Better, return a "container" that offers an entry-pair
               iterator, which generates on-line from the stored
               leaf-pairs.
               => This avoids storing all entry-pairs (potentially a
                  lot) in an array just to iterate them afterwards,
                  encapsulating the "foreach e1 in l1{ foreach e2 in l2 {
                  do X } } double loop
    */
    // SIMPLE version
    bool Test( const bvh_type& bvh2, std::vector< std::pair<entry_index_type,entry_index_type> >& vec_pair_entry_id ) const
        {
            vec_pair_entry_id.clear();
            if( IsEmpty() ) return false;
            TestRecursive( bvh2, GetRoot(), bvh2.GetRoot(), vec_pair_entry_id );
            return vec_pair_entry_id.size() > 0;
        }

    // TEPD version \note TEPD is only invoked for leaf-nodes/entries, NOT for internal-nodes.
    bool Test( const bvh_type& bvh2, const TEPD& tepd, std::vector< std::pair<entry_index_type,entry_index_type> >& vec_pair_entry_id ) const
        {
            vec_pair_entry_id.clear();
            if( IsEmpty() ) return false;
            TestRecursive( bvh2, GetRoot(), bvh2.GetRoot(), tepd, vec_pair_entry_id );
            return vec_pair_entry_id.size() > 0;
        }

    // Self-collision test \todo Consider ESCC, BN-Tree or some other optimization
    bool TestSelf( std::vector< std::pair<entry_index_type,entry_index_type> >& vec_pair_entry_id ) const
        {
            if( IsEmpty() ) return false;
            return false;
        }

    /*\todo Consider _Preorder, _Postorder, etc... versions
      \todo Use node-stack to iterate tree in given order while Delegate returns true
    */
    void Map( MND mnd ) const
        {
            if( IsEmpty() ) return;
            m_stackPairNID.clear(); //\todo USE Scratchpad-Stack to avoid new/delete!
            m_stackPairNID.push_back( std::make_pair( GetRootIndex(), 0 ) );
            while( !m_stackPairNID.empty() )
            {
                const Node& node( m_vecNode[m_stackPairNID.back().first] );
                uint32 level( m_stackPairNID.back().second );
                m_stackPairNID.pop_back();
                if( mnd( *this, node, level )
                    && !node.IsLeaf() )
                {
                    m_stackPairNID.push_back( std::make_pair( node.m_Topology.m_Left, level+1 ) );
                    m_stackPairNID.push_back( std::make_pair( node.m_Topology.m_Right, level+1 ) );
                }
            }
      }


    /* SIMPLE raycast on BVH that returns ALL pierced entries.
       \note NO early out on ray clipping
       \note NO exact ray hits on UNKWNOWN entry geometry
      \todo CONSIDER using Map() with a proper delegate instead

      \todo CONSIDER TRED-version test-ray-entry-delegate, given an
            entry, perform an exact ray-cast on its contents and build
            a valid GRayHit<D>, allowing for ray-clipping/early-out
            with ray_dir ordered L/R traversal
    */
    bool RayCast( const vec_type& ray_pos, const vec_type& ray_dir, const Interval& interval,
                  std::vector<entry_index_type>& vec_entry_id ) const
        {
            vec_entry_id.clear();
            if( IsEmpty() ) return false;
            //\todo Recursive, use (non-precomputed) Node-Sorting to prioritize "closest" child L/R
            m_stackNID.clear(); //\todo Consider Scratchpad-Stack to avoid new/delete!
            m_stackNID.push_back( GetRootIndex() );
            while( !m_stackNID.empty() )
            {
                const Node& node( m_vecNode[m_stackNID.back()] );
                m_stackNID.pop_back();
                if( TestRay<cDimension>( node.m_Geometry.m_BV, ray_pos, ray_dir, interval ) ) //\todo <cDimension> required due to template param MISMATCH between GAABB<unsigned> and GVec<int>... MUST FIX THIS AT ONCE, IT HAPPENS EVERYWHERE!
                {
                    if( node.IsLeaf() )
                    {
                        //\todo If ever 1 leaf > 1 entry, add all leaf entries to vec_entry_id
                        vec_entry_id.push_back( node.m_Topology.m_FirstEID );
                    }
                    else
                    {
                        m_stackNID.push_back( node.m_Topology.m_Left );
                        m_stackNID.push_back( node.m_Topology.m_Right );
                    }
                }
            }
            return vec_entry_id.size() > 0;
        }

    finline uint32 GetHeight() const { return m_Height; }
    finline bool IsEmpty() const { return m_NumEntries == 0; }

public: //TEMP: to drawdebug...
    struct node_topology
    {
        node_index_type m_Left, m_Right; //Required for top-down Test() \todo If preorder, Left = Parent+1 and only Right needs to be saved!
        entry_count_type m_FirstEID, m_NumEntries; //|todo If only leafs need this, reuse m_Left/m_Right instead! Even better, first/num can be completely implicit if leaf == entry and leafs stored contiguously, entry_index = &node[entry_index] - &node[0]!!
        //node_index_type m_Parent; //\todo Required for bottom-up refit, but can be avoided if nodes in post-order so that node < parent(node)
    };
    struct leaf_entries
    {
        entry_count_type m_FirstEID, m_NumEntries;
    };
    struct node_geometry
    {
        bv_type m_BV;
    };
    struct Node
    {
        node_topology m_Topology; //\todo union node_topology, leaf_entries?
        node_geometry m_Geometry;
        //finline bool IsLeaf() const { return m_Topology.m_NumEntries > 0; }
        finline bool IsLeaf() const { return m_Topology.m_Left == m_Topology.m_Right; } //\todo BETTER, does not require num_entries to exist
    };

    uint32 m_NumEntries;
    uint32 m_Height;
    std::vector<Node> m_vecNode;

    uint32 m_GeometryTimeStamp; //!< Last known *external* geometry timestamp, so that we can *externally* know if the BVH is up-to-date or not

    //\todo Consider Scratchpad-Stack to avoid new/delete!
    mutable std::vector< uint32 > m_stackNID;
    mutable std::vector< std::pair<uint32,uint32> > m_stackPairNID;
    //

    //std::vector<entry_index_type> m_vecEntryID; // All entry id partitioned into node subarrays \todo If leaf==entry => NO NEED to store indices explicitly...

    /*\todo Consider this, where Leaves and Internal nodes may have DIFFERENT bv_type and topology
    struct leaf_type
    {
        node_index_type m_Parent; \todo Could be ignored if BV are stored, but required if not to refit parent BV otherwise
        entry_index_type m_EntryId; \todo Could be implicit (the same) in leaf_type index in m_vecLeaf array
        leaf_bv_type m_BV;
    };
    struct node_type
    {
        entry_count_type m_FirstEntry, m_NumEntries;
        node_index_type m_Left, m_Right; //Required for top-down Test()
        node_bv_type m_BV;
    };
    std::vector<leaf_type> m_vecLeaf;
    std::vector<node_type> m_vecNode;
    */

public:
    uint32 GetRootIndex() const { return m_vecNode.size()-1; } //\todo THIS MAY CHANGE
    const bvh_type::Node& GetRoot() const { return m_vecNode[ GetRootIndex() ]; }
    const bvh_type::Node& GetNode( uint32 nid ) const { return m_vecNode[nid]; }

private:
    struct bottom_up_node
    {
        entry_index_type m_FirstEID;
        entry_count_type m_NumEntries;
        node_index_type m_Left, m_Right;
        bv_type m_BV;
        inline bottom_up_node() {}
    };

    /* BottomUp brute-force implementation inspired in RTCD pg 246,
       but avoiding explicit pointers/new/delete.
       \note Cost is O(N^3) in num_entries
       \note Nodes are merged greedily merge according to a given heuristic
       \note Merged nodes are added at the end and therefore the
             sequential node layout guarantees that all children appear
             before their parent, but is NOT strict post-order.
    */
    void Rebuild_BottomUp_N3( entry_count_type num_entries, const GEBVD& gebvd )
        {
            m_NumEntries = num_entries;
            GEO_ASSERT( num_entries > 0 && 2*num_entries-1 <= cMaxNodes );
            // Add all leaf-nodes and compute their BV
            std::vector<bottom_up_node> vecBUN(num_entries);
            std::vector<node_index_type> vecOpen(num_entries);
            for( unsigned int it_entry=0; it_entry<num_entries; it_entry++ )
            {
                bottom_up_node& bun( vecBUN[it_entry] );
                bun.m_FirstEID = it_entry;
                bun.m_NumEntries = 1;
                bun.m_Left = bun.m_Right = cInvalidNodeId;
                gebvd( it_entry, bun.m_BV );
                vecOpen[it_entry] = it_entry;
            }
            // While open nodes, merge best open node pair according to heuristic
            unsigned int num_merges(0);
            while( vecOpen.size() > 1 )
            {
                // Find best open node pair to merge
                std::pair<node_index_type,node_index_type> best_oid(cInvalidNodeId,cInvalidNodeId);
                // = FindBestMerge( vecBUN, vecOpen ); \todo Consider different strategies
                {
                    Real best_h( mal::Infinity<Real>() );
                    for( uint32 i=0; i<vecOpen.size()-1; i++ )
                    {
                        for( uint32 j=i+1; j<vecOpen.size(); j++ )
                        {
                            Real h;
                            // = Heuristic( vecBUN[vecOpen[i]].m_BV, vecBUN[vecOpen[j]].m_BV ); \todo Consider different heuristics
                            {
                                h = bv::ComputeVolume( bv_type( vecBUN[vecOpen[i]].m_BV ).Merge( vecBUN[vecOpen[j]].m_BV ) );
                            }
                            if( h < best_h )
                            {
                                best_h = h;
                                best_oid = std::make_pair(i,j);
                            }
                        }
                    }
                }
                // Add merge node at the end
                vecBUN.push_back( bottom_up_node() );
                bottom_up_node& bun( vecBUN.back() );
                bun.m_Left = vecOpen[best_oid.first];
                bun.m_Right = vecOpen[best_oid.second];
                bun.m_FirstEID = cInvalidNodeId;
                bun.m_NumEntries = 0;
                bun.m_BV = bv_type( vecBUN[bun.m_Left].m_BV ).Merge( vecBUN[bun.m_Right].m_BV );
                // Close child nodes (\note As second > first, MUST reuse slots in this order, otherwise second may be overwritten by first after slot reuse
                vecOpen[ best_oid.second ] = vecOpen.back(); vecOpen.pop_back();
                vecOpen[ best_oid.first ] = vecOpen.back(); vecOpen.pop_back();
                // Open merge node
                vecOpen.push_back( vecBUN.size()-1 );
                num_merges++;
            }
            // GEO_LOG("Rebuild_BottomUp_N3() Built with %d entries, %d nodes, %d merges", (int)num_entries, (int)vecBUN.size(), (int)num_merges );

            // Bake ST into the GBoundingVolumeHierarchy_ST_DG
            //GEO_LOG("Baking...");
            m_vecNode.resize( vecBUN.size() );
            for( unsigned int it_node=0; it_node<vecBUN.size(); it_node++ )
            {
                Node& node( m_vecNode[it_node] );
                bottom_up_node& bun( vecBUN[it_node] );
                node.m_Topology.m_FirstEID = bun.m_FirstEID;
                node.m_Topology.m_NumEntries = bun.m_NumEntries;
                node.m_Topology.m_Left = bun.m_Left;
                node.m_Topology.m_Right = bun.m_Right;
                node.m_Geometry.m_BV = bun.m_BV;
            }

            //---- Common part
            // Compute level
            uint32 max_level(0);
            Map( [&max_level]( const bvh_type& bvh, const bvh_type::Node& node, uint32 level )
                 {
                     if( max_level < level ) max_level = level;
                     return true;
                 }
            );
            m_Height = max_level+1;
            // GEO_LOG("Rebuild_BottomUp_N3() Tree height = %d", (int)m_Height );
            // Reserve stacks
            m_stackNID.reserve(m_Height);
            m_stackPairNID.reserve(2*m_Height);
        }

    /* Similar to BF, but instead of searching all open node pairs
       O(N^2), pre-generate all topologic pairs O(N) and update them
       lazily/amortized using an MF-set scheme when pairs are merged.
       \note Cost is O(N^2) in num_entries
       \note When an open pair (e_i,e_j) is selected and merged, all
             remaining open pairs that involve *any* node in the same set as
             e_i or e_j become out of date.
       \todo REQUIRES topology!!
       \todo See if a heap-like scheme is possible to find best open
             pair in O(log(N)) instead of O(N) WITHOUT having to
             update all of them when pairs are merged
             => Actually, out-of-date pairs can only have an heuristic
                too small, so when the best pair is extracted from the
                heap *and* updated, if it was out-of-date, it should
                be reinserted into the heap and a different best pair
                extracted until an up-to-date one is found
    */
    void Rebuild_BottomUp_N2( entry_count_type num_entries, const GEBVD& gebvd )
        {
            //\todo
        }

    struct top_down_node
    {
        entry_index_type m_FirstEID;
        entry_count_type m_NumEntries;
        node_index_type m_Left, m_Right;
        bv_type m_BV;
        //inline top_down_node() {}
    };

    /* Very basic algorithm:
       - First pass: Sort entries/leafnodes recursively along alternating axis 0,1,2...
         - \todo TRY using std::nth_element O(N) instead of std::sort O(N log N)
       - Second pass: Build actual internal nodes recursively, without reordering leaves or child nodes AT ALL
       - This 2-pass alg avoids the need for an additional array of entry indices that can be reordered during internal node construction
       - Produces non-strict post-order layout, as required by Refit()
    */
    void Rebuild_TopDown_NLogN( entry_count_type num_entries, const GEBVD& gebvd )
        {
            m_NumEntries = num_entries;
            GEO_ASSERT( num_entries > 0 && 2*num_entries-1 <= cMaxNodes );
            // GEO_LOG("Rebuild_TopDown_NLogN() computing BV...");
            // Add all leaf-nodes and compute their BV
            std::vector<top_down_node> vecTDN(num_entries);
            for( unsigned int it_entry=0; it_entry<num_entries; it_entry++ )
            {
                top_down_node& tdn( vecTDN[it_entry] );
                tdn.m_FirstEID = it_entry;
                tdn.m_NumEntries = 1;
                tdn.m_Left = tdn.m_Right = cInvalidNodeId;
                gebvd( it_entry, tdn.m_BV );
            }
            // Recursive call
            // GEO_LOG("Rebuild_TopDown_NLogN() Recursive...");
#define __ENABLE_TOPDOWN_SORT_HACK
#ifdef __ENABLE_TOPDOWN_SORT_HACK
            Rebuild_TopDown_NLogN_Recursive_SortLeaves( vecTDN, 0, vecTDN.size(), 0 );
#endif
            Rebuild_TopDown_NLogN_Recursive_CreateInternal( vecTDN, 0, vecTDN.size() );
            // GEO_LOG("Rebuild_TopDown_NLogN_Recursive() Built with %d entries, %d nodes", (int)num_entries, (int)vecTDN.size() );

            // Bake ST into the GBoundingVolumeHierarchy_ST_DG
            // GEO_LOG("Rebuild_TopDown_NLogN() Baking...");
            m_vecNode.resize( vecTDN.size() );
            for( unsigned int it_node=0; it_node<vecTDN.size(); it_node++ )
            {
                Node& node( m_vecNode[it_node] );
                top_down_node& tdn( vecTDN[it_node] );
                node.m_Topology.m_FirstEID = tdn.m_FirstEID;
                node.m_Topology.m_NumEntries = tdn.m_NumEntries;
                node.m_Topology.m_Left = tdn.m_Left;
                node.m_Topology.m_Right = tdn.m_Right;
                node.m_Geometry.m_BV = tdn.m_BV;
            }

#define __ENABLE_GBVH_DEBUG_TOPDOWN
#ifdef __ENABLE_GBVH_DEBUG_TOPDOWN //TEMP: This helped catching a bug, keep it for a while
            std::vector<uint32> vec_num_parents( m_vecNode.size(), 0 );
            for( unsigned int it_node=0; it_node<m_vecNode.size(); it_node++ )
            {
                GEO_ASSERT( m_vecNode[it_node].IsLeaf() || m_vecNode[it_node].m_Topology.m_Left < it_node );
                GEO_ASSERT( m_vecNode[it_node].IsLeaf() || m_vecNode[it_node].m_Topology.m_Right < it_node );
                GEO_ASSERT( m_vecNode[it_node].IsLeaf() || m_vecNode[it_node].m_Topology.m_FirstEID == cInvalidNodeId );
                if( !m_vecNode[it_node].IsLeaf() )
                {
                    vec_num_parents[m_vecNode[it_node].m_Topology.m_Left]++;
                    vec_num_parents[m_vecNode[it_node].m_Topology.m_Right]++;
                }
            }
            for( unsigned int it_node=0; it_node<vec_num_parents.size(); it_node++ )
            {
                if( it_node == vec_num_parents.size()-1 )
                    GEO_ASSERT( vec_num_parents[it_node] == 0 );
                else
                    GEO_ASSERT( vec_num_parents[it_node] == 1 );
            }
#endif

            //---- Common part
            // Compute level
            uint32 max_level(0);
            Map( [&max_level]( const bvh_type& bvh, const bvh_type::Node& node, uint32 level )
                 {
                     if( max_level < level ) max_level = level;
                     return true;
                 }
            );
            m_Height = max_level+1;
            // GEO_LOG("Rebuild_TopDown_NLogN() Tree height = %d", (int)m_Height );

            // Reserve stacks
            m_stackNID.reserve(m_Height);
            m_stackPairNID.reserve(2*m_Height);
        }

    /* Lazy programmer recursive sort along alternating axis.
       \note Currently limited to 3 coordinate axis, but could accept
             num_axis and a lambda <= predicate as params to support
             arbitrary axis (KDOP)

       IMPORTANT: MUST use STRICTLY WEAK ordering <. If <= is used,
                  std::sort() CAUSES MEMORY CORRUPTION, see
                  https://schneide.wordpress.com/2010/11/01/bug-hunting-fun-with-stdsort/
                  This caused mem corruption and segfaults AFTER the
                  actual sort, in "unrelated"
                  Rebuild_TopDown_NLogN_Recursive_CreateInternal()
                  new() calls on push_back()
     */
    void Rebuild_TopDown_NLogN_Recursive_SortLeaves( std::vector<top_down_node>& vec_tdn, uint32 left_entry, uint32 right_entry, uint32 axis )
        {
            uint32 count_entries = right_entry-left_entry;
            if( count_entries > 1 )
            {
                std::sort( vec_tdn.begin()+left_entry, vec_tdn.begin()+right_entry,
                           [axis]
                           (const top_down_node& n1, const top_down_node& n2)
                           { return n1.m_BV.GetPos()[axis] < n2.m_BV.GetPos()[axis]; } ); //IMPORTANT: MUST be a STRICTLY WEAK ordering <
                uint32 mid_entry = left_entry + count_entries / 2;
                Rebuild_TopDown_NLogN_Recursive_SortLeaves( vec_tdn, left_entry, mid_entry, (axis+1)%cDimension );
                Rebuild_TopDown_NLogN_Recursive_SortLeaves( vec_tdn, mid_entry, right_entry, (axis+1)%cDimension );
            }
        }

    /* Lazy programmer recursive midpoint splitting of vec_tdn.
       \note This is NOT an auto-balancing top-down classification, as
             it does not use ANY geometric heuristic on the nodes/leaves.
             - The reason is that the entries are STRICTLY leaf nodes,
               and there is no "array of entry ids" that we can
               reorder/split as we descend. We only have the vec_tdn
               where new nodes can be added safely, but NOT
               moved/sorted, as previous nodes store indices that
               would be silently invalidated.
             - However, if the initial entry/leaves ARE pre-classified
               (as done in
               Rebuild_TopDown_NLogN_Recursive_SortLeaves()), the tree
               is auto-balanced.
    */
    uint32 Rebuild_TopDown_NLogN_Recursive_CreateInternal( std::vector<top_down_node>& vec_tdn, uint32 left_entry, uint32 right_entry )
        {
            uint32 count_entries = right_entry-left_entry;
            if( count_entries > 1 )
            {
                // Split entry/leaves sub-array
                uint32 mid_entry = left_entry + count_entries / 2;
                uint32 left_nid = Rebuild_TopDown_NLogN_Recursive_CreateInternal( vec_tdn, left_entry, mid_entry );
                uint32 right_nid = Rebuild_TopDown_NLogN_Recursive_CreateInternal( vec_tdn, mid_entry, right_entry );
                // Add parent node at the end to ensure post-order
                vec_tdn.push_back( top_down_node() );
                top_down_node& tdn( vec_tdn.back() );
                tdn.m_Left = left_nid;
                tdn.m_Right = right_nid;
                tdn.m_FirstEID = cInvalidNodeId;
                tdn.m_NumEntries = count_entries;
                tdn.m_BV = bv_type( vec_tdn[tdn.m_Left].m_BV ).Merge( vec_tdn[tdn.m_Right].m_BV );
                return vec_tdn.size()-1;
            }
            else // if( count_entries == 1 )
            {
                // This must be an entry/leaf, just return it's index
                return left_entry;
            }
        }

private:
    /*! We use the "generic informed simpled depth-first traversal" strategy from RTCD pg.256
      \todo See stack-based iterative version also in RTCD
     */
    void TestRecursive( const bvh_type& bvh2, const Node& n1, const Node& n2, std::vector< std::pair<entry_index_type,entry_index_type> >& vec_pair_entry_id ) const
        {
            GEO_ASSERT(!IsEmpty());
            if( TestOverlap( n1.m_Geometry.m_BV, n2.m_Geometry.m_BV ) )
            {
                if( n1.IsLeaf() && n2.IsLeaf() )
                {
                    //\todo If ever 1 leaf > 1 entry, add all leaf entries to vec_entry_id
                    vec_pair_entry_id.push_back( std::make_pair( n1.m_Topology.m_FirstEID, n2.m_Topology.m_FirstEID ) );
                }
                else if( n2.IsLeaf()  //\note this is the DescendLarger(n1,n2) strategy in RTCD pg.257, we use the _Fast() methods because exact KDOP volume is VERY expensive to compute
                         || (!n1.IsLeaf() && bv::ComputeVolume_Fast(n1.m_Geometry.m_BV) >= bv::ComputeVolume_Fast(n2.m_Geometry.m_BV)) )
                {
                    TestRecursive( bvh2, m_vecNode[n1.m_Topology.m_Left], n2, vec_pair_entry_id );
                    TestRecursive( bvh2, m_vecNode[n1.m_Topology.m_Right], n2, vec_pair_entry_id );
                }
                else
                {
                    TestRecursive( bvh2, n1, bvh2.m_vecNode[n2.m_Topology.m_Left], vec_pair_entry_id );
                    TestRecursive( bvh2, n1, bvh2.m_vecNode[n2.m_Topology.m_Right], vec_pair_entry_id );
                }
            }
        }
    /*! We use the "generic informed simpled depth-first traversal" strategy from RTCD pg.256
      \todo See stack-based iterative version also in RTCD
    */
    void TestRecursive( const bvh_type& bvh2, const Node& n1, const Node& n2, const TEPD& tepd, std::vector< std::pair<entry_index_type,entry_index_type> >& vec_pair_entry_id ) const
        {
            GEO_ASSERT(!IsEmpty());
            if( n1.IsLeaf() && n2.IsLeaf() )
            {
                //\todo If ever 1 leaf > 1 entry, add all leaf entries to vec_entry_id
                if( tepd( n1.m_Topology.m_FirstEID, n2.m_Topology.m_FirstEID, n1.m_Geometry.m_BV, n2.m_Geometry.m_BV ) )
                    vec_pair_entry_id.push_back( std::make_pair( n1.m_Topology.m_FirstEID, n2.m_Topology.m_FirstEID ) );
            }
            else if( TestOverlap( n1.m_Geometry.m_BV, n2.m_Geometry.m_BV ) )
            {
                if( n2.IsLeaf()  //\note this is the DescendLarger(n1,n2) strategy in RTCD pg.257, we use the _Fast() methods because exact KDOP volume is VERY expensive to compute
                    || (!n1.IsLeaf() && bv::ComputeVolume_Fast(n1.m_Geometry.m_BV) >= bv::ComputeVolume_Fast(n2.m_Geometry.m_BV)) )
                {
                    TestRecursive( bvh2, m_vecNode[n1.m_Topology.m_Left], n2, tepd, vec_pair_entry_id );
                    TestRecursive( bvh2, m_vecNode[n1.m_Topology.m_Right], n2, tepd, vec_pair_entry_id );
                }
                else
                {
                    TestRecursive( bvh2, n1, bvh2.m_vecNode[n2.m_Topology.m_Left], tepd, vec_pair_entry_id );
                    TestRecursive( bvh2, n1, bvh2.m_vecNode[n2.m_Topology.m_Right], tepd, vec_pair_entry_id );
                }
            }
        }
};

/* \todo Consider renaming it "query" instead of iterator...
    template <typename BVHT1, typename BVHT2>
    class bvh2bvh_overlap_iterator
    {
        typedef std::pair< bvh_type1::entry_index_type, bvh_type2::entry_index_type > result_type;
        bvh2bvh_overlap_iterator( T1& bvh1, T2& bvh2, TBVPD tbvpd, TEPD tbvpd )
        : m_rBVH1(bvh1), m_BVH2(bvh2)
        , m_TBVPD(tbvpd)
        , m_Overlap(bvh_type1::cInvalidEntryIndex,bvh_type2::cInvalidEntryIndex)
            {}
        void operator++() { find next overlap or finish; }
        bool operator bool() const { return !m_Stack.empty(); }
        const result_type& operator*() const { return m_Overlap; }
    private:
        const bvh_type1& m_rBVH1;
        const bvh_type2& m_rBVH2;
        TBVPD/TEPD m_TBVPD;
        \todo m_stackPairNID;
        result_type m_Overlap;
    };

\todo Fully generic test supports different BV and entry types, as well as specific TEPD and a TBVPD
//\todo TBVPD: supports bv1!=bv2 and different ref sys
//\todo TEPD: leaf-only, supports all TBVPD cases as well as entry-specific info not available to the BVH (eg: BSlab/BDOP)
template< typename BVHT1, typename BVHT2 >
bool Test( const BVHT1& bvh1, const BVHT1& bvh2,
           const boost::function<bool (const BVHT1::bv_type& bv1, const BVHT2::bv_type& bv2)>& tbvpd,
           const boost::function<bool (BVHT1::entry_index_type eid1, BVHT2::entry_index_type eid2, const BVHT1::bv_type& bv1, const BVHT2::bv_type& bv2)>& tepd,
           std::vector< std::pair<typename BVHT1::entry_index_type,typename BVHT1::entry_index_type> > vec_overlaps ) const
{
}
*/

//----------------------------------------------------------------
// Specific instantiations
//----------------------------------------------------------------
typedef GBoundingVolumeHierarchy_ST_DG<geo::bv::Sphere2,geo::feature_index_type> BVH_ST_DG_Sphere2;
typedef GBoundingVolumeHierarchy_ST_DG<geo::bv::AABB2,geo::feature_index_type> BVH_ST_DG_AABB2;
typedef GBoundingVolumeHierarchy_ST_DG<geo::bv::DOP2_K8,geo::feature_index_type> BVH_ST_DG_DOP2_K8;
typedef GBoundingVolumeHierarchy_ST_DG<geo::bv::Sphere3,geo::feature_index_type> BVH_ST_DG_Sphere3;
typedef GBoundingVolumeHierarchy_ST_DG<geo::bv::AABB3,geo::feature_index_type> BVH_ST_DG_AABB3;

}} //namespace geo::bv

#endif // GEO_BV_GBOUNDING_VOLUME_HIERARCHY_H
