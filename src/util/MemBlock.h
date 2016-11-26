#ifndef UTIL_MEM_BLOCK_H
#define UTIL_MEM_BLOCK_H

#include <util/Config.h>
#include <memory.h>

namespace util {

/*! Simple memory block with fixed max size.
  Memory can be either internally or externally allocated:
  - Internal alloc: Destructor and Clear will deallocate it. 4b-aligned.
  - External alloc: No allocation nor destruction.
*/
class MemBlock
{
public:
    MemBlock( uint32 size_in_bytes = 0, int8 *p_data = 0 ) : m_Size(0), m_pData(0), m_bHasInternalAllocation(false) { Init( size_in_bytes, p_data ); }
    virtual ~MemBlock() { Clear(); }

    inline void Clear() { if( HasInternalAllocation() ) delete[] m_pData; m_Size = 0; m_pData = 0; }
    inline void Init( uint32 size_in_bytes = 0, int8 *p_data = 0 )
        {
            Clear();
            if( p_data == 0 && size_in_bytes > 0 ) { m_pData = (int8*)new int32[(size_in_bytes+3)/4]; m_Size = size_in_bytes; m_bHasInternalAllocation = true; }
            else if( p_data != 0 && size_in_bytes > 0 ) { m_pData = p_data; m_Size = size_in_bytes; m_bHasInternalAllocation = false; }
            else { m_pData = 0; m_Size = 0; m_bHasInternalAllocation = false; }
        }
    inline bool Realloc( uint32 size_in_bytes )
        {
            UTIL_ASSERT( size_in_bytes >= m_Size
                         && (0 == m_Size || m_bHasInternalAllocation) );
            int8* pNewData = (int8*)new int32[(size_in_bytes+3)/4]; //\todo this should be a PROPER realloc... use malloc()/realloc() to improve performance!?
            if( 0 != pNewData )
            {
                if( 0 != m_pData )
                {
                    memcpy( pNewData, m_pData, m_Size );
                    delete[] m_pData;
                }
                m_Size = size_in_bytes;
                m_pData = pNewData;
                m_bHasInternalAllocation = true;
                return true;
            }
            else
                return false;
        }

    inline uint32 GetSize() const { return m_Size; }
    inline int8 *GetData() { return m_pData; }
    inline const int8 *GetData() const { return m_pData; }
    inline bool HasInternalAllocation() const { return m_bHasInternalAllocation; }

    inline int8 &operator[]( uint32 offset ) { UTIL_ASSERT( offset >= 0 && offset < (uint32)m_Size ); return m_pData[offset]; }
    inline const int8 &operator[]( uint32 offset ) const { UTIL_ASSERT( offset >= 0 && offset < (uint32)m_Size ); return m_pData[offset]; }

private:
    uint32 m_Size;
    int8* m_pData;
    bool m_bHasInternalAllocation;
};

} //util

#endif //UTIL_MEM_BLOCK_H
