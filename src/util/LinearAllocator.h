#ifndef UTIL_LINEAR_ALLOCATOR_H
#define UTIL_LINEAR_ALLOCATOR_H

#include <util/Config.h>
#include <util/MemBlock.h>

//#define __UTIL_ENABLE_LINEAR_ALLOCATOR_TRACE

namespace util {

class LinearAllocator : public MemBlock
{
public:
    //enum EConstants { cAlignment = 16 };
    enum EConstants { cAlignment = 4 };
public:
    inline LinearAllocator( uint32 size_in_bytes, int8* p_data = 0 )
    : m_MB(size_in_bytes,p_data), m_pEnd(m_MB.GetData()) {}
    inline ~LinearAllocator() {}

    inline void Clear() { m_MB.Clear(); m_pEnd = 0; }

    inline void* Malloc( uint32 size )
        {
            if( size == 0 )
            {
#ifdef __UTIL_ENABLE_LINEAR_ALLOCATOR_TRACE
                UTIL_LOG_ERROR( "LinearAllocator::Malloc(0) works but dangerous!" );
#endif
                return m_pEnd;
            }
            uint32 aligned_size = ComputeAlignedSize(size);
            if( aligned_size <= GetAvailableBytes() )
            {
#ifdef __UTIL_ENABLE_LINEAR_ALLOCATOR_TRACE
                UTIL_LOG( "LinearAllocator::Malloc(%u) => %d in [%p,%p]", size, aligned_size, m_pEnd, m_pEnd+aligned_size );
#endif
                void* pBegin = m_pEnd;
                m_pEnd += aligned_size;
                return pBegin;
            }
            else
            {
                UTIL_LOG_ERROR( "LinearAllocator::Malloc(%u) out of memory!", size );
                return 0;
            }
        }

    inline void Rewind( void* p )
        {
            //\todo IsAligned() may fail for the FIRST frame, as m_MB.m_pData alignment is NOT FORCED to cAlignment...
            UTIL_ASSERT( IsAligned(p) );
            UTIL_ASSERT( p >= m_MB.GetData() );
            UTIL_ASSERT( p <= m_pEnd ); //\note Equality required to support 0-size allocations, otherwise it's an error
#ifdef __UTIL_ENABLE_LINEAR_ALLOCATOR_TRACE
            UTIL_LOG( "LinearAllocator::Rewind(%p)", p );
#endif
            m_pEnd = reinterpret_cast<int8*>(p);
        }

    // POD and POD[N] do not call ctor/dtor
    template <typename T>
    inline T* NewPOD() { return reinterpret_cast<T*>(Malloc(sizeof(T))); }
    template <typename T>
    inline T* NewArrayPOD( uint32 count ) { return reinterpret_cast<T*>(Malloc(count*sizeof(T))); }

    inline uint32 GetAvailableBytes() const { return m_MB.GetData() + m_MB.GetSize() - m_pEnd; }
    inline void* GetBeginPtr() { return m_MB.GetData(); }
    inline void* GetEndPtr() { return m_pEnd; }

protected:
    inline uint32 ComputeAlignedSize( uint32 size ) const { return (size + (cAlignment - 1)) & ~(cAlignment - 1); }
    inline bool IsAligned( void* p ) const { return 0 == size_t(p) % cAlignment; }

private:
    MemBlock m_MB;
    int8* m_pEnd;
};

} //util

#endif //UTIL_LINEAR_ALLOCATOR_H
