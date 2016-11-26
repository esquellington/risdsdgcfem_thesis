#ifndef GEO_BP_G_BIPARTITION_TREE_H
#define GEO_BP_G_BIPARTITION_TREE_H

#include <algorithm>
#include <vector>
#include <util/GPoolDA.h>

namespace geo {
namespace bp {

enum ETraverseFunctorResult {
    eBPT_TraverseStop  = 0,
    eBPT_TraverseLeft  = 1,
    eBPT_TraverseRight = 2,
    eBPT_TraverseBoth  = (eBPT_TraverseLeft | eBPT_TraverseRight)
};

/*! Generic Bipartition Tree with configurable:
  - Entry type
  - Node base type
  - Classification criterion
  - Splitting criterion

  UserNodeT must fulfill the following non-virtual interface:
  - bool User_NeedsSplit( const EntryT *vec_entry, uint32 size, int level, void *user_data )
    - return true if must split, false otherwise
  - bool User_Classify( const EntryT &entry, int level, void *user_data );
    - return true if Left, false if Right
  - void User_Init( UserNodeT *p_parent, const EntryT *vec_entry, uint32 size, int level, void *user_data );
    - allows user-init when a node is created
  - void User_SetChildren( UserNodeT *p_left, UserNodeT *p_right, int level, void *user_data );
    - notifies user of newly created child nodes (during a split)

  \note While this specific tree implementation may not be the most
  efficient for all cases (kd-tree,k-DOP tree,BSP,etc...), it's good
  enough as a reference implementation that can be instantiated,
  derived or just copied and modified. Therefore, we *must* keep it as
  SIMPLE and GENERIC as possible.
*/
template <typename EntryT, typename UserNodeT>
class GBipartitionTree
{
public:
    typedef EntryT entry_type;
    typedef UserNodeT user_node_type;

    //! BPT node
    struct node_type: public user_node_type
    {
        uint32 m_FirstIdx;
        uint32 m_Size;
        uint32 m_Free;
        node_type *m_pParent;
        node_type *m_pLeft;
        node_type *m_pRight;
            
        inline void Init( node_type *p_parent, uint32 first_idx, uint32 size, uint32 capacity )
            {
                m_FirstIdx = first_idx; m_Size = size; m_Free = capacity-m_Size;
                m_pParent = p_parent; m_pLeft = 0; m_pRight = 0;
            }
        inline uint32 Size() const { return m_Size; }
        inline uint32 Capacity() const { return m_Size+m_Free; }
        inline bool HasFreeSlots() const { return (m_Free > 0); }
        inline bool IsLeaf() const { return 0 == m_pLeft && 0 == m_pRight; } //\todo Could opt as 0 == m_pLeft + m_pRight
        inline void Clear()
            {
                m_FirstIdx = m_Size = m_Free = 0;
                m_pParent = m_pLeft = m_pRight = 0;
            }
    };

    /*! Forward entry iterator
      It's just a pointer to the current element but stores also the end pointer for safety.    
    */
    class entry_iterator
    {
    public:
        inline entry_iterator() : m_pCurrent(0), m_pEnd(0) {}
        inline entry_iterator( entry_type *p_begin, entry_type *p_end ) : m_pCurrent(p_begin), m_pEnd(p_end) {}
        inline const entry_type &operator*() const { UTIL_ASSERT(IsValid()); return *m_pCurrent; } 
        inline const entry_type *operator->() const { UTIL_ASSERT(IsValid()); return m_pCurrent; }
        inline entry_type &operator*() { UTIL_ASSERT(IsValid()); return *m_pCurrent; } 
        inline entry_type *operator->() { UTIL_ASSERT(IsValid()); return m_pCurrent; }
        //! Simple increment. May return INVALID entries when free slots are allowed
        inline void operator++() { UTIL_ASSERT(IsValid()); m_pCurrent++; }
        //! Fast-skip free slots. MAY FAIL if right subtree is empty but not null \todo Maybe we SHOULD avoid splitting if a child is empty!
        //  inline void operator++() { UTIL_ASSERT(IsValid()); do { m_pCurrent++; } while(0==*m_pCurrent); }
        //! Expensive-skip free slots with bound-check. Handles empty non-null right subtrees correctly
        //  inline void operator++() { UTIL_ASSERT(IsValid()); do { m_pCurrent++; } while(0==*m_pCurrent && IsValid()); }
        inline bool operator==( const entry_iterator &other ) const { return (m_pCurrent == other.m_pCurrent); }
        inline bool IsValid() const { return (m_pCurrent < m_pEnd); }
    private:
        entry_type *m_pCurrent;
        entry_type *m_pEnd;
    };

    //\todo class node_iterator
  
public:
    GBipartitionTree( unsigned int capacity, unsigned int max_nodes )
    : m_Capacity(capacity)
    , m_poolNodes(max_nodes)
    , m_UserData(0)
        {
            m_vecEntries = new entry_type[m_Capacity];
            m_RootNode.Clear();
        }
    ~GBipartitionTree()
        {
            if( m_vecEntries ) delete m_vecEntries;            
        }
    
    //! Sets user-provided stuff
    void SetUserData( void *user_data ) { m_UserData = user_data; }
    
    //!\name Entry-related methods
    //@{
    //! Add raw entry, no classification/splitting
    void Add( const entry_type &entry )
        {
            if( m_Size < m_Capacity ) m_vecEntries[m_Size++] = entry;
            else ;//Realloc
        }
        
    //! Try to add without rebuilding
    bool Insert( entry_type& entry )
        {
            /*
              FindNode();
              node.Insert( entry );
              node.UserInsert( entry );
            */
            return true;
        }

    //! remove without rebuilding
    bool Remove( entry_type& entry )
        {
            /*
              FindNode();
              node.Remove( entry );
              node.UserRemove( entry );
            */
            return true;
        }

    //! find entry
    /*
      node_type *Find( entry_type& entry )
      {
            
      }
    */
    //@}

    //!\name Node access methods
    //@{
    const user_node_type *GetRoot() const { return &m_RootNode; }
    const user_node_type *GetLeft( const user_node_type &unode ) const
        { return static_cast<const node_type&>(unode).m_pLeft; }
    const user_node_type *GetRight( const user_node_type &unode ) const
        { return static_cast<const node_type&>(unode).m_pRight; }
    
    user_node_type *GetRoot() { return &m_RootNode; }
    user_node_type *GetLeft( const user_node_type &unode )
        { return static_cast<const node_type&>(unode).m_pLeft; }
    user_node_type *GetRight( const user_node_type &unode )
        { return static_cast<const node_type&>(unode).m_pRight; }
    
    entry_iterator GetEntryIterator( const user_node_type &unode ) const
        {
            const node_type node( static_cast<const node_type&>(unode) );
            return entry_iterator( &m_vecEntries[node.m_FirstIdx],
                                   &m_vecEntries[node.m_FirstIdx+node.m_Size] );
        }  
    //@}
  
    //!\name Whole-Tree related methods
    //@{
    //! Rebuild whole tree
    void Rebuild()
        {
            // Init root
            m_RootNode.Init( 0, 0, m_Size, m_Capacity );
            m_RootNode.User_Init( 0, &m_vecEntries[0], m_Size, 0, m_UserData );

            // Split recursively according to user-provided method
            //\todo use static alloc stack or at least reuse across different calls!
            std::vector< std::pair<node_type*,int> > Stack;
            Stack.push_back( std::make_pair(&m_RootNode,int(0)) );
            while( !Stack.empty() )
            {
                node_type &node( *Stack.back().first );
                int level( Stack.back().second );
                Stack.pop_back();
                if( node.User_NeedsSplit( &m_vecEntries[node.m_FirstIdx], node.m_Size, level, m_UserData ) )
                {
                    Split( node, level );
                    if( 0 != node.m_pLeft ) Stack.push_back( std::make_pair(node.m_pLeft,level+1) );
                    if( 0 != node.m_pRight ) Stack.push_back( std::make_pair(node.m_pRight,level+1) );
                }
            }
        }        
  
    /*! Maps a functor to tree nodes in preorder until the functor
      returns false or the whole tree is visited. Returns number of
      visited nodes.
      
      This can be used for any classification-independent tree-traversal
    */
    template <typename FunctorT>
    unsigned int TraversePreorder( FunctorT &functor )
        {
            //\todo use static alloc stack or at least reuse across different calls!
            std::vector< std::pair<node_type*,int> > Stack;
            Stack.push_back( std::make_pair(&m_RootNode,int(0)) );
            unsigned int num_visited(0);
            while( !Stack.empty() )                
            {
                node_type &node( *Stack.back().first );
                int level( Stack.back().second );
                Stack.pop_back();
                num_visited++;
                int32 result = functor( static_cast<const user_node_type&>(node),
                                        &m_vecEntries[node.m_FirstIdx], node.m_Size,
                                        level, m_UserData );
                if( eBPT_TraverseStop != result )
                {
                    if( 0 != node.m_pLeft && 0 != (eBPT_TraverseLeft & result) )
                        Stack.push_back( std::make_pair(node.m_pLeft,level+1) );
                    if( 0 != node.m_pRight && 0 != (eBPT_TraverseRight & result) )
                        Stack.push_back( std::make_pair(node.m_pRight,level+1) );                    
                }
                else
                    return num_visited;
            }
            return num_visited;
        }

    void Clear()
        {
            m_Size = 0;
            m_poolNodes.Clear();
            m_RootNode.Clear();
        }

    //! Resets with a new capacity
    void Reset( unsigned int capacity )
        {
            if( capacity > m_Capacity )
            {
                delete m_vecEntries;
                m_Capacity = capacity;
                m_vecEntries = new entry_type[m_Capacity];
            }
            Clear();
        }    
    //@}
  
private:
    /*! Single-pass split of:
      node = [ first..last : free ]
      into:
      [ left : free/2 : right ]
      where
      left = [ first_l..last_l : free/4 ]
      right = [ first_r..last_r : free/4 ]
    */
    void Split( node_type &node, int level )
        {
            GEO_ASSERT( node.IsLeaf() );
            
            // Split entry sub-array
            int32 last_left( node.m_FirstIdx );
            int32 first_right( node.m_FirstIdx + node.m_Size );
            while( last_left < first_right )
            {
                // Classify according to user-provided callback
                if( node.User_Classify( m_vecEntries[last_left], level, m_UserData ) ) last_left++;
                else std::swap( m_vecEntries[last_left], m_vecEntries[--first_right] );
            }

            // Create Left and/or Right childs
            int32 size_left( last_left - node.m_FirstIdx );
            if( size_left > 0 )
            {
                node.m_pLeft = m_poolNodes.New();
                node.m_pLeft->Init( &node, node.m_FirstIdx, size_left, size_left );
                node.m_pLeft->User_Init( &node, &m_vecEntries[node.m_FirstIdx], size_left, level+1, m_UserData );
            }

            int32 size_right( node.m_FirstIdx + node.m_Size - first_right );
            if( size_right > 0 )
            {
                node.m_pRight = m_poolNodes.New();
                node.m_pRight->Init( &node, first_right, size_right, size_right );
                node.m_pRight->User_Init( &node, &m_vecEntries[first_right], size_right, level+1, m_UserData );
            }

            // Set user children
            node.User_SetChildren( node.m_pLeft, node.m_pRight, level, m_UserData );
        }
        
private:        
    //!\name Buffer management
    //@{
    uint32 m_Size;
    uint32 m_Capacity;
    entry_type *m_vecEntries;
    //@}

    //!\name Tree management
    //@{
    node_type m_RootNode;
    util::GPoolDA<node_type> m_poolNodes;
    //@}

    //! \name User stuff
    //@{
    void *m_UserData;
    //@}

    // Iteration \todo IMPORTANT: Distingir entre leaf_iterator i
    //tree_iterator. El primer es SEQUENCIAL, mentre que el segon
    //ha de saltar-se els espais FREE intermitjos, o bé a saco, o
    //bé recursivament!!
};

}} //namespace geo::bp

#endif // GEO_BP_G_BIPARTITION_TREE_H
