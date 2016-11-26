#ifndef UTIL_GINDEXED_REFCOUNTED_POOL_DA_H
#define UTIL_GINDEXED_REFCOUNTED_POOL_DA_H

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
template < typename T, typename IndexT = unsigned int, typename RefCountT = uint8 >
class GIndexedRefCountedPoolDA
{
public:
    typedef T entry_type;
    typedef IndexT index_type;
    typedef RefCountT refcount_type;
    
public:
    GIndexedRefCountedPoolDA( size_t hint_max_size = 0 )
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
        default: GEO_LOG_ASSERT(false, "Unsupported index_type size %d", (int)sizeof(index_type) ); break;
        }
        UTIL_ASSERT( hint_max_size <= m_MaxIndex );
        switch( sizeof(refcount_type) )
        {
        case 1: m_MaxRC = refcount_type(0xFF); break;
        case 2: m_MaxRC = refcount_type(0xFFFF); break;
        case 4: m_MaxRC = refcount_type(0xFFFFFFFF); break;
        case 8: m_MaxRC = refcount_type(0xFFFFFFFFFFFFFFFF); break;
        default: GEO_LOG_ASSERT(false, "Unsupported refcount_type size %d", (int)sizeof(refcount_type) ); break;
        }        
        Clear();
        if( m_HintMaxSize > 0 )
        {
            m_vecEntries.reserve( m_HintMaxSize );
            m_vecRC.reserve( m_HintMaxSize );
        }
    }
    ~GIndexedRefCountedPoolDA() { Clear(); }

    void Trace() const
    {
        UTIL_LOG( "GIRCPDA<>::Trace() S = %d, E = %d, RC = %d, F = %d",
                  m_Size, (int)m_vecEntries.size(), (int)m_vecRC.size(), (int)m_vecFree.size() );
    }
    
    const T &operator[]( index_type index ) const { UTIL_ASSERT( IsValid(index) ); return m_vecEntries[index]; }
    T &operator[]( index_type index ) { UTIL_ASSERT( IsValid(index) ); return m_vecEntries[index]; }
    
    inline void Clear()
    {
#ifdef PROFILE_DEBUG
        //\todo Check if all entries have RC <= 1, otherwise we're deleting referenced data!
        unsigned int acc_ref_count(0);
        for( unsigned int i=0; i<m_vecRC.size(); i++ ) acc_ref_count += m_vecRC[i];
        UTIL_LOG_ERROR_IF( acc_ref_count != m_Size,
                           "GIndexedRefCountedPoolDA::Clear() %d != %d, probably deleting referenced data!", acc_ref_count, m_Size );        
        /*
        UTIL_LOG_ASSERT( acc_ref_count == m_Size,
                         "GIndexedRefCountedPoolDA::Clear() %d != %d, probably deleting referenced data!", acc_ref_count, m_Size );
        */
#endif
        m_vecEntries.clear();
        m_vecFree.clear();        
        m_vecRC.clear();
        m_Size = 0;
    }
    
    //! \name Object creation/destruction
    //@{
    index_type New( const T &value = T() )
    {
        /*
        UTIL_LOG_WARNING_IF( m_Size > m_HintMaxSize,
                        "GIndexedRefCountedPoolDA: hint_max_size = %d surpassed, a GPoolSA would have crashed here",
                        m_HintMaxSize );
        */
        UTIL_ASSERT( m_Size < m_MaxIndex );
        m_Size++;
        if( !m_vecFree.empty() )
        {
            index_type index = m_vecFree.back();
            m_vecFree.pop_back();
            m_vecEntries[index] = value;
            m_vecRC[index] = refcount_type(1);
            //UTIL_LOG("New %d",index);
            return index;
        }
        else
        {
            index_type index = m_vecEntries.size();
            m_vecEntries.push_back( value );
            m_vecRC.push_back( refcount_type(1) );
            //UTIL_LOG("New %d",index);
            return index;
        }
    }
    
    inline void Delete( index_type index )
    {
        //UTIL_LOG("Delete %d",index);        
        UTIL_ASSERT( IsValid(index) );
        UTIL_ASSERT( m_vecRC[index] == refcount_type(1) );
        m_Size--;
        m_vecRC[index] = refcount_type(0);
        m_vecFree.push_back( index );
    }

    inline void IncRef( index_type index )
    {
        //UTIL_LOG("IncRef %d",index);        
        UTIL_ASSERT( IsValid(index) );
        UTIL_ASSERT( m_vecRC[index] < m_MaxRC );
        m_vecRC[index]++;
    }

    inline void DecRef( index_type index )
    {
        //UTIL_LOG("DecRef %d",index);
        UTIL_ASSERT( IsValid(index) );
        UTIL_ASSERT( m_vecRC[index] > refcount_type(0) );
        m_vecRC[index]--;
        if( m_vecRC[index] == refcount_type(0) )
        {
            m_Size--;
            m_vecFree.push_back( index );
        }
    }
    
    inline bool IsValid( index_type index ) const
    {
        // In range and not free (=> refcount > 0)
        return index >= 0
            && index < (index_type)m_vecRC.size()
            && m_vecRC[index] > refcount_type(0);
        //\note SLOW and unnecessary && m_vecFree.end() == std::find( m_vecFree.begin(), m_vecFree.end(), index );
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
    refcount_type m_MaxRC;
    std::vector<T> m_vecEntries;
    std::vector<refcount_type> m_vecRC;
    std::vector<index_type> m_vecFree;
};

} // namespace util

#endif // UTIL_GINDEXED_REFCOUNTED_POOL_DA_H
