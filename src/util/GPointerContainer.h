#ifndef UTIL_GPOINTERCONTAINER_H
#define UTIL_GPOINTERCONTAINER_H

#include <util/Config.h>

#define __USE_TR1_UNORDERED_SET
#ifdef __USE_TR1_UNORDERED_SET
#  include <tr1/unordered_set>
#else
#  include <list>
#endif

namespace util {

/*! Simple Add/Remove/Find container for T* elements.
  - Implements a nice forward iterator that does NOT require comparison with .end()
  - Due to unordered_set, all ops are O(1) or amortized O(1)
  
  \note I'm not sure if unordered_set iterators are invalidated after any
  insertion/removal... msdn says yes.
*/
template <typename T>
class GPointerContainer
{
public:
    class iterator;
#ifdef __USE_TR1_UNORDERED_SET
    typedef std::tr1::unordered_set<T*> internal_container_type;
#else
    typedef std::list<T*> internal_container_type;
#endif
    
public:
    GPointerContainer() {}

#ifdef __USE_TR1_UNORDERED_SET
    inline void Add( T *p_item ) { m_IC.insert(p_item); }
    inline void Remove( T *p_item ) { m_IC.erase(p_item); }
    inline iterator Find( const T *p_item ) const { return iterator(*this,m_IC.find(const_cast<T*>(p_item))); } //\todo const_cast is a HACK req to compile
#else
    inline void Add( T *p_item ) { m_IC.push_back(p_item); }
    inline void Remove( T *p_item ) { m_IC.remove(p_item); }
    inline iterator Find( const T *p_item ) const { return iterator(*this,std::find(m_IC,p_item)); }
#endif
    inline void Clear() { return m_IC.clear(); }
    inline iterator Begin() const { return iterator(*this); }
    inline unsigned int Size() const { return m_IC.size(); }
    inline bool IsEmpty() const { return m_IC.empty(); }
    
public:
    //! Simple forward iterator
    class iterator
    {
    public:
        typedef typename internal_container_type::const_iterator internal_iterator_type;
        
    private:
        iterator( const GPointerContainer &container ) : m_pContainer(&container), m_Iterator(container.m_IC.begin()) {}
        iterator( const GPointerContainer &container, const internal_iterator_type &iterator )
        : m_pContainer(&container), m_Iterator(iterator) {}
        friend class GPointerContainer;
        
    public:
        iterator() : m_pContainer(0) {}
        iterator( const iterator &it ) : m_pContainer(it.m_pContainer), m_Iterator(it.m_Iterator) {}

        inline const T* operator->() const { UTIL_ASSERT(IsValid()); return *m_Iterator; }
        inline const T& operator*() const { UTIL_ASSERT(IsValid()); return *(*m_Iterator); }
        inline operator const T*() const { UTIL_ASSERT(IsValid()); return *m_Iterator; }
        
        inline T* operator->() { UTIL_ASSERT(IsValid()); return *m_Iterator; }
        inline T& operator*() { UTIL_ASSERT(IsValid()); return *(*m_Iterator); }
        inline operator T*() { UTIL_ASSERT(IsValid()); return *m_Iterator; }
        
        inline void operator ++() { ++m_Iterator; }
        inline bool IsValid() const { return 0 != m_pContainer && m_Iterator != m_pContainer->m_IC.end(); }
        
    private:
        const GPointerContainer *m_pContainer;
        internal_iterator_type m_Iterator;
    };
    
private:
    internal_container_type m_IC;
    friend class iterator;    
};
    
} // namespace util

#endif // UTIL_GPOINTERCONTAINER_H
