#ifndef UTIL_GPOOL_SA_H
#define UTIL_GPOOL_SA_H

#include <util/Config.h>

#define __ENABLE_ISVALID_TRACKING //!< Performs extra checks on IsValid(), but is slower and yelds bigger sizeof(Entry)

namespace util
{

/*! Iterable pool of T items with strict capacity max_items (<2^16)
    and O(1) runtime operations.

  Memory is statically allocated on construction.
    
  Pool functionality:
  - Creation is O(N)
  - New is O(1)
  - Delete is O(1)
  - Size() is O(1)

  Iterator functionality:
  - Begin() is O(1)
  - *Iterator is O(1)
  - ++Iterator is O(1)
*/
template <typename T>
class GPoolSA
{
public:
    class iterator;
    enum EConstants { cInvalidIdx = uint16(0xFFFF) };
    
public:
    GPoolSA( unsigned int max_items )
    : m_Buffer(0)
    , m_MaxItems(uint16(max_items))
    , m_Size(0)
    {
        UTIL_ASSERT( m_MaxItems > 1 && m_MaxItems < cInvalidIdx );
        m_Buffer = new Entry[max_items];
        Clear();
    }
    ~GPoolSA() { if(m_Buffer) delete[] m_Buffer; }

    void Clear()
    {
        m_Size = 0;
        m_FirstValid = cInvalidIdx;
        m_FirstFree = 0;
        for( uint16 i=0; i<m_MaxItems; i++ )
        {
            m_Buffer[i].m_Prev = (i-1+m_MaxItems) % m_MaxItems;
            m_Buffer[i].m_Next = (i+1) % m_MaxItems;
#ifdef __ENABLE_CHECK_ISVALID            
            m_Buffer[i].m_IsValid = false;
#endif            
        }
    }
    
    //! \name Object creation/destruction
    //@{
    T *New()
    {
        UTIL_ASSERT(m_Size < m_MaxItems);

        uint16 idx = m_FirstFree;
#ifdef __ENABLE_CHECK_ISVALID
        m_Buffer[idx].m_IsValid = true;
#endif
            
        // Remove from free list
        if( m_Buffer[m_FirstFree].m_Next == m_FirstFree )
            m_FirstFree = cInvalidIdx;
        else
        {
            // Recompute first free and relink
            m_FirstFree = m_Buffer[idx].m_Next;
            m_Buffer[m_Buffer[idx].m_Prev].m_Next = m_FirstFree;
            m_Buffer[m_FirstFree].m_Prev = m_Buffer[idx].m_Prev;
        }

        // Add to valid list
        if( cInvalidIdx == m_FirstValid )
        {
            // Set as single valid
            m_FirstValid = idx;
            m_Buffer[idx].m_Prev = idx;
            m_Buffer[idx].m_Next = idx;
        }
        else
        {
            // Link to neighbours
            m_Buffer[idx].m_Prev = m_Buffer[m_FirstValid].m_Prev;
            m_Buffer[idx].m_Next = m_FirstValid;
            // Link from neighbours
            m_Buffer[ m_Buffer[idx].m_Prev ].m_Next = idx;
            m_Buffer[ m_Buffer[idx].m_Next ].m_Prev = idx;            
        }
        m_Size++;
        return &(m_Buffer[idx].m_Item);
    }
    
    inline void Delete( T *p_item )
    {
        UTIL_ASSERT( IsValid(p_item) );
        Entry *lEntry = (Entry*)p_item;
        uint16 idx = uint16(lEntry - m_Buffer);
#ifdef __ENABLE_CHECK_ISVALID
        m_Buffer[idx].m_IsValid = false;
#endif

        // Remove from valid list
        if( m_Buffer[idx].m_Next == idx )
            m_FirstValid = cInvalidIdx;
        else
        {
            // Recompute first valid and relink
            m_FirstValid = m_Buffer[idx].m_Next;
            m_Buffer[m_Buffer[idx].m_Prev].m_Next = m_FirstValid;
            m_Buffer[m_FirstValid].m_Prev = m_Buffer[idx].m_Prev;
        }
        
        // Add to free list
        if( cInvalidIdx == m_FirstFree )
        {
            // Set as single free
            m_FirstFree = idx;
            m_Buffer[m_FirstFree].m_Prev = idx;
            m_Buffer[m_FirstFree].m_Next = idx;
        }
        else
        {
            // Link to neighbours
            m_Buffer[idx].m_Prev = m_Buffer[m_FirstFree].m_Prev;
            m_Buffer[idx].m_Next = m_FirstFree;
            // Link from neighbours
            m_Buffer[ m_Buffer[idx].m_Prev ].m_Next = idx;
            m_Buffer[ m_Buffer[idx].m_Next ].m_Prev = idx;
        }
        
        m_Size--;
    }

    inline bool IsValid( const T* p_item ) const
    {
        Entry *p_entry = (Entry*)p_item;
        bool bIsValid( p_entry >= m_Buffer
                       && p_entry < m_Buffer + m_MaxItems
#ifdef __ENABLE_CHECK_ISVALID
                       && m_Buffer[ p_entry - m_Buffer ].m_IsValid );
#else
                       );
#endif
        return bIsValid;
    }
    //@}
    
    //! \name Object iteration (valid objects only)
    //@{
    inline bool IsEmpty() const { return 0 == m_Size; }
    inline unsigned int Size() const { return m_Size; }
    inline unsigned int Capacity() const { return m_MaxItems; }
    inline iterator Begin() const { return iterator(*this); }
    inline iterator Find( const T *p_item ) const { return iterator(*this,p_item); }
    //@}

public:
    //! Simple forward iterator
    class iterator
    {
    private:
        iterator( const GPoolSA &pool ) : m_pPool(&pool), m_Index(pool.m_FirstValid) {} //!< Sets to first valid idx (may be -1)
        iterator( const GPoolSA &pool, const T* p_item ) //!< Create iterator pointing to a specific item
        : m_pPool(&pool)
        {
            UTIL_ASSERT( pool.IsValid(p_item) );
            const Entry *pEntry = reinterpret_cast<const Entry*>(p_item);
            m_Index = uint16(pEntry - pool.m_Buffer);
        }        
        friend class GPoolSA;
        
    public:
        iterator() : m_pPool(0), m_Index(cInvalidIdx) {}
        iterator( const iterator &it ) : m_pPool(it.m_pPool), m_Index(it.m_Index) {}
        
        inline const T* operator->() const { UTIL_ASSERT(IsValid()); return &m_pPool->m_Buffer[m_Index].m_Item; }
        inline const T& operator*() const { UTIL_ASSERT(IsValid()); return m_pPool->m_Buffer[m_Index].m_Item; }
        inline const T* GetPtr() const { UTIL_ASSERT(IsValid()); return &m_pPool->m_Buffer[m_Index].m_Item; }
        
        inline T* operator->() { UTIL_ASSERT(IsValid()); return &m_pPool->m_Buffer[m_Index].m_Item; }
        inline T& operator*() { UTIL_ASSERT(IsValid()); return m_pPool->m_Buffer[m_Index].m_Item; }
        inline T* GetPtr() { UTIL_ASSERT(IsValid()); return &m_pPool->m_Buffer[m_Index].m_Item; }

        inline const GPoolSA *GetPool() const { return m_pPool; }
        
        inline void operator++()
        {
            UTIL_ASSERT(IsValid());
            
            // If its the last element the iterator becomes invalid, otherwise advance
            if( m_pPool->m_Buffer[m_Index].m_Next == m_pPool->m_FirstValid )
                m_Index = cInvalidIdx;
            else
                m_Index = m_pPool->m_Buffer[m_Index].m_Next;
        }
        
        inline bool IsValid() const { return m_Index != cInvalidIdx; }
        
    private:
        const GPoolSA *m_pPool;
        uint16 m_Index;
    };
    
private:
    struct Entry
    {
        T m_Item;
        uint16 m_Prev, m_Next;
#ifdef __ENABLE_CHECK_ISVALID
        bool m_IsValid;
#endif
    };
    Entry *m_Buffer;
    uint16 m_MaxItems;
    uint16 m_Size;
    uint16 m_FirstValid;
    uint16 m_FirstFree;

    friend class iterator;
};

} // namespace util

#endif // UTIL_GPOOL_SA_H
