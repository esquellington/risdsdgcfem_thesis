#ifndef UTIL_GSIMPLETREE_H
#define UTIL_GSIMPLETREE_H

#include <util/Config.h>
#include <list>

namespace util {

/*! Simple Tree container with *aggregation*-based use.

  \note Client ItemT class does NOT need to *inherit* from
  GSimpleTree<ItemT>, just *aggregate* it as an attribute and set
  themselves as the node Item. Therefore, a single Item can be
  referenced in several trees!
*/
template <typename ItemT>
class GSimpleTree
{
public:
    class iterator; //!\todo?
    typedef ItemT item_type;
    typedef GSimpleTree<ItemT> node_type;

private:
    item_type *m_pItem;
    node_type *m_pParent;
    node_type *m_pFirstChild;
    node_type *m_pNextSibling;

public:
    GSimpleTree() : m_pItem(0), m_pParent(0), m_pFirstChild(0), m_pNextSibling(0) {}

    inline bool IsRoot() const { return 0 == m_pParent; }
    inline bool IsLeaf() const { return 0 == m_pFirstChild; }
    
    inline void Reset() { m_pItem = 0; m_pParent = 0; m_pFirstChild = 0; m_pNextSibling = 0; }
    inline void SetItem( item_type *p_item ) { m_pItem = p_item; }
    inline void SetParent( node_type *p_parent ) { m_pParent = p_parent; }
    inline void SetNextSibling( node_type *p_sibling ) { m_pNextSibling = p_sibling; }    
    inline void AddChild( node_type *p_child )
        {
            UTIL_ASSERT( p_child->GetParent() == 0 );
            // Add at the end of sibling list
            if( 0 == m_pFirstChild )        
                m_pFirstChild = p_child;
            else
            {
                node_type *pNode(m_pFirstChild);
                while( 0 != pNode->GetNextSibling() )
                    pNode = pNode->GetNextSibling();
                pNode->SetNextSibling(p_child);
            }            
            p_child->SetParent(this);
        }
    inline void RemoveChild( node_type *p_child )
        {
            RGDK_ASSERT( !IsLeaf() && p_child->GetParent() == this );
            if( p_child == m_pFirstChild )
                m_pFirstChild = p_child->GetNextSibling();
            else
            {
                // Find node prev to p_child, and relink it to node succ to p_child
                node_type *pPred(m_pFirstChild);
                while( pPred->GetNextSibling() != p_child ) pPred = pPred->GetNextSibling();
                pPred->SetNextSibling( p_child->GetNextSibling() );
            }
        }

    inline item_type *GetItem() const { return m_pItem; }
    inline node_type *GetParent() const { return m_pParent; }
    inline node_type *GetFirstChild() const { return m_pFirstChild; }
    inline node_type *GetNextSibling() const { return m_pNextSibling; }    
};
    
} // namespace util

#endif // UTIL_GSIMPLETREE_H
