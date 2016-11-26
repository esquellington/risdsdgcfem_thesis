#ifndef UTIL_GPOINTERLIST_H
#define UTIL_GPOINTERLIST_H

#include <util/Config.h>
#include <list>

namespace util {

/*! Simple Add/Remove/Find list for T* elements.
  - Implements a nice forward iterator that does NOT require comparison with .end()
*/
template <typename T>
class GPointerList
{
public:
    class iterator;
    typedef std::list<T*> internal_list_type;
    
public:
    GPointerList() {}

    inline void Add( T *p_item ) { m_IC.push_back(p_item); }
    inline void Remove( T *p_item ) { m_IC.remove(p_item); }
    inline iterator Find( T *p_item ) const { return iterator(*this,std::find(m_IC,p_item)); }
    inline void Clear() { return m_IC.clear(); }
    inline iterator Begin() const { return iterator(*this); }
    inline unsigned int Size() const { return m_IC.size(); }
    
public:
    //! Simple forward iterator
    class iterator
    {
    public:
        typedef typename internal_list_type::const_iterator internal_iterator_type;
        
    private:
        iterator( const GPointerList &list ) : m_pList(&list), m_Iterator(list.m_IC.begin()) {}
        iterator( const GPointerList &list, const internal_iterator_type &iterator )
        : m_pList(&list), m_Iterator(iterator) {}
        friend class GPointerList;
        
    public:
        iterator() : m_pList(0) {}
        iterator( const iterator &it ) : m_pList(it.m_pList), m_Iterator(it.m_Iterator) {}

        inline const T* operator->() const { UTIL_ASSERT(IsValid()); return *m_Iterator; }
        inline const T& operator*() const { UTIL_ASSERT(IsValid()); return *(*m_Iterator); }
        inline operator const T*() const { UTIL_ASSERT(IsValid()); return *m_Iterator; }
        
        inline T* operator->() { UTIL_ASSERT(IsValid()); return *m_Iterator; }
        inline T& operator*() { UTIL_ASSERT(IsValid()); return *(*m_Iterator); }
        inline operator T*() { UTIL_ASSERT(IsValid()); return *m_Iterator; }
        
        inline void operator ++() { ++m_Iterator; }
        inline bool IsValid() const { return 0 != m_pList && m_Iterator != m_pList->m_IC.end(); }
        
    private:
        const GPointerList *m_pList;
        internal_iterator_type m_Iterator;
    };
    
private:
    internal_list_type m_IC;
    friend class iterator;
};
    
} // namespace util

#endif // UTIL_GPOINTERLIST_H
