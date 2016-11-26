#ifndef UTIL_TRIAD_H
#define UTIL_TRIAD_H

#include <util/Config.h>

namespace util {

/*! std::pair + 1
  Must fulfill:
  - Assignable (default op== and op(const &)
  Optionally:
  - DefaultConstructible
  - EqualityComparable
  - LessThanComparable
*/
template <typename T1, typename T2, typename T3>
struct triad
{
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;

    triad() {}
    triad( const T1 &a, const T2 &b, const T3 &c ): first(a), second(b), third(c) {}

    T1 first;
    T2 second;
    T3 third;
};

template <typename T1, typename T2, typename T3>
triad<T1,T2,T3> make_triad( const T1 &a, const T2 &b, const T3 &c ) { return triad<T1,T2,T3>(a,b,c); }

template <typename T1, typename T2, typename T3>
bool operator==( const triad<T1,T2,T3>& a, const triad<T1,T2,T3>& b )
{
    return ( a.first == b.first && a.second == b.second && a.third == b.third );
}

template <typename T1, typename T2, typename T3>
bool operator!=( const triad<T1,T2,T3>& a, const triad<T1,T2,T3>& b )
{
    return !( a.first == b.first && a.second == b.second && a.third == b.third );
}

template <typename T1, typename T2, typename T3>
bool operator<(const triad<T1,T2,T3>& a, const triad<T1,T2,T3>& b )
{
    if( a.first  < b.first  ) return true;
    if( b.first  < a.first  ) return false;
    if( a.second < b.second ) return true;
    if( b.second < a.second ) return false;
    return a.third < b.third;
}

} //util

#endif // UTIL_TRIAD_H
