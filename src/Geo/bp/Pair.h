#ifndef GEO_BP_PAIR_H
#define GEO_BP_PAIR_H

#include <Geo/Config.h>
#include <Geo/bp/Proxy.h>

namespace geo { namespace bp
{

/*! BP overlapping pair */
struct Pair
{
    enum EState { eInvalid = 0, eNew, ePersistent, eVanished };
    
    //\name Attribs managed by IBroadPhase
    //@{    
    Proxy *m_pProxy1;
    Proxy *m_pProxy2;
    mutable void *m_UserData;
    //@}
    
    //\name Attribs managed by PairContainer
    //@{
    uint32 m_TimeStamp; //!< Timestamp and m_State could be merged into 32 bit only, yielding 16 byte structure alignment
    EState m_State;
    //\todo USe type TimeStamp24 that leaves 8 higher bits free for State, for example
    //@

    //! \name Invalidation of proxy/object references (when Obj is deleted externally)
    //@{
    inline void Invalidate1() { m_pProxy1 = 0; }
    inline void Invalidate2() { m_pProxy2 = 0; }
    inline bool IsValid1() const { return 0 != m_pProxy1; }
    inline bool IsValid2() const { return 0 != m_pProxy2; }
    //@}
            
    //! \name Access to per-object data (hidden in Proxy)
    //@{
    template <typename T> inline T GetObj1() const { return reinterpret_cast<T>(m_pProxy1->m_pObj); }
    template <typename T> inline T GetObj2() const { return reinterpret_cast<T>(m_pProxy2->m_pObj); }
    //@}

    //! \name Access user-data
    //@{
    template <typename T> inline void SetUserData( T user_data ) const { m_UserData = reinterpret_cast<void*>(user_data); }
    template <typename T> inline T GetUserData() const { return reinterpret_cast<T>(m_UserData); }
    //@}
};

}} //namespace geo::bp

#endif // GEO_BP_PAIR_H
