#ifndef UTIL_SCOPED_ALLOCATOR_H
#define UTIL_SCOPED_ALLOCATOR_H

#include <util/Config.h>
#include <util/LinearAllocator.h>

//#define __UTIL_ENABLE_SCOPED_ALLOCATOR_TRACE

namespace util {

/*Scoped Allocator (aka Frame Allocator)

  Receives an external LinearAllocator/Scratchpad and defines a
  self-rewinding scope inside it.

  \todo Non-POD data is correctly created/destroyed using placement
        new and a finalizer stack.

  \note See doc at http://dice.se/publications/scope-stack-allocation/
  and code at http://pastebin.com/h7nU8JE2 (\note there's a POTENTIAL
  BUG in that code at line 24: u8 *result = m_ptr += size)

  \note Se also http://www.codingwisdom.com/codingwisdom/2012/09/you-need-a-frame-allocator.html
*/
class ScopedAllocator
{
public:
    inline ScopedAllocator( LinearAllocator& la, const char* name )
    : m_rLinearAllocator(la)
    , m_pBeginScope(la.GetEndPtr())
#ifdef __UTIL_ENABLE_SCOPED_ALLOCATOR_TRACE
    , m_Name(name)
#endif
        {
#ifdef __UTIL_ENABLE_SCOPED_ALLOCATOR_TRACE
            UTIL_LOG( "%s Push at %p", m_Name.GetStr(), m_pBeginScope );
#endif
        }
    inline ~ScopedAllocator()
        {
#ifdef __UTIL_ENABLE_SCOPED_ALLOCATOR_TRACE
            UTIL_LOG( "%s Pop to %p", m_Name.GetStr(), m_pBeginScope );
#endif
            m_rLinearAllocator.Rewind(m_pBeginScope);
        }

    // POD and POD[N] do not call ctor/dtor
    template <typename T>
    inline T* NewPOD() { return m_rLinearAllocator.NewPOD<T>(); }
    template <typename T>
    inline T* NewArrayPOD( uint32 count ) { return m_rLinearAllocator.NewArrayPOD<T>(count); }

    /* \todo New Object > POD
       - Call default (placement) constructors
       - Save destructor register, call destructors on pop
      template <typename T>
      inline T* NewObject() {  }
    */
    /* \todo New Array Object > POD
       - Call default (placement) constructors
       - Save destructor register
      template <typename T>
      inline T* NewArrayObject() {}
    */

private:
    LinearAllocator& m_rLinearAllocator;
    void* m_pBeginScope;
#ifdef __UTIL_ENABLE_SCOPED_ALLOCATOR_TRACE
    String64 m_Name;
#endif
};

} //util

#endif //UTIL_SCOPED_ALLOCATOR_H
