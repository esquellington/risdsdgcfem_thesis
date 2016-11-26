#ifndef UTIL_GINDEXED_POOL_DA_H
#define UTIL_GINDEXED_POOL_DA_H

#include <util/Config.h>
#include <vector>
#include <algorithm>

namespace util
{

/*! Non-Iterable pool of T items with indexed access, dynamic
    allocation and O(1) runtime operations.
  
  Pool functionality:
  - Creation is O(1)
  - New is O(1)
  - Delete is O(1)
  - Size() is O(1)
*/
template < typename T, typename IndexT = unsigned int >
class GIndexedPoolDA
{
public:
    typedef T entry_type;
    typedef IndexT index_type;
    
public:
    GIndexedPoolDA( size_t hint_max_size = 0 )
    : m_HintMaxSize( index_type(hint_max_size) )
    , m_Size(0)
    {
        // Compute max_index according to index_type size, \todo WE ASSUME unsigned, SHOULD use some kind of numeric_limits<index_type>::max()!!
        switch( sizeof(index_type) )
        {
        case 1: m_MaxIndex = index_type(0xFF); break;
        case 2: m_MaxIndex = index_type(0xFFFF); break;
        case 4: m_MaxIndex = index_type(0xFFFFFFFF); break;
        case 8: m_MaxIndex = index_type(0xFFFFFFFFFFFFFFFF); break;
        default: GEO_ASSERT(false); break;
        }
        UTIL_ASSERT( hint_max_size <= m_MaxIndex );
        Clear();
        if( m_HintMaxSize > 0 ) m_vecEntries.reserve( m_HintMaxSize );
    }
    ~GIndexedPoolDA() { Clear(); }
    
    const T &operator[]( index_type index ) const { UTIL_ASSERT( IsValid(index) ); return m_vecEntries[index]; }
    T &operator[]( index_type index ) { UTIL_ASSERT( IsValid(index) ); return m_vecEntries[index]; }
    
    inline void Clear()
    {
        m_vecEntries.clear();
        m_vecFree.clear();
        m_Size = 0;
    }
    
    //! \name Object creation/destruction
    //@{
    index_type New( const T &value = T() )
    {
        /*
        UTIL_LOG_WARNING_IF( m_Size > m_HintMaxSize,
                        "GIndexedPoolDA: hint_max_size = %d surpassed, a GPoolSA would have crashed here",
                        m_HintMaxSize );
        */
        UTIL_ASSERT( m_Size < m_MaxIndex );
        m_Size++;
        if( !m_vecFree.empty() )
        {
            index_type index = m_vecFree.back();
            m_vecFree.pop_back();
            m_vecEntries[index] = value;
            return index;
        }
        else
        {
            index_type index = m_vecEntries.size();
            m_vecEntries.push_back( value );
            return index;
        }
    }
    
    inline void Delete( index_type index )
    {
        UTIL_ASSERT( IsValid(index) );
        m_Size--;
        m_vecFree.push_back( index );
    }
    
    inline bool IsValid( index_type index ) const
    {
        // In range and not free (SLOW! just for debug, CONSIDER adding a bitfield to flag entry validity)
        return index >= 0
            && index < (index_type)m_vecEntries.size()
            && m_vecFree.end() == std::find( m_vecFree.begin(), m_vecFree.end(), index );
    }
    //@}
    
    //! \name Sizes
    //@{
    inline bool IsEmpty() const { return 0 == m_Size; }
    inline unsigned int Size() const { return m_Size; }
    inline unsigned int Capacity() const { return m_MaxIndex; }
    inline unsigned int NumAllocated() const { return m_vecEntries.size(); }
    //@}
    
    //!\pre Requires T::operator<()
    //\todo void Sort()
    
private:
    index_type m_HintMaxSize;
    index_type m_Size;
    index_type m_MaxIndex;
    std::vector<T> m_vecEntries;
    std::vector<index_type> m_vecFree;
};

} // namespace util

#endif // UTIL_GINDEXED_POOL_DA_H
