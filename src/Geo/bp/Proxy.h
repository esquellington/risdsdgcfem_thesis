#ifndef GEO_BP_PROXY_H
#define GEO_BP_PROXY_H

#include <Geo/Config.h>
#include <Geo/bv/BoundingVolume.h>

namespace geo { namespace bp
{

/*! BP object proxy

  Contains:
  - A void* to the user-object it represents
  - Flags that define the proxy nature
  - Generic Bounding Volume with internal Timestamp that changes on each modification
  - Timestamp of the last time it accessed the BV

  \note Specific BP may derive Proxy should NOT expose the derived type to BP users.
*/
struct Proxy
{
    //!\name Attribs (keep it small)
    //@{
    void *m_pObj;              //!< User-defined object
    uint16 m_PartId;           //!< part id for multi-part objects 
    Flags8 m_Flags;            //!< Internal flags
    Flags8 m_UserFlags;        //!< User-defined flags
    uint32 m_TimeStamp;        //!< Last change TS (Flags+TS could be easily merged into a single uint32(8+24)
    const bv::IBoundingVolume *m_pBV; //!< Generic BV, allocated elsewhere
    //\todo Use type TimeStamp24 that leaves 8 higher bits free for Flags, for example, or even TimeStamp16 with 8+8 bits free for internal+user flags
    //@}
    
    //! BP proxy flags
    enum EFlags {
        // Generic-BP flags
        eDefault  = 0,       //!< !Boundary && !Static
        eBoundary = (1<<0),  //!< Boundary proxy (no B2B tests)
        eStatic   = (1<<1),  //!< Static objects allow some optimizations in IBroadPhase::Update()
        eNestedBP = (1<<2),  //!< Nested BP proxy, whose m_pObj is an IBroadPhase.
        // Specific-BP flags
        // eBig, eHuge, eFast, ....
    };
    
    Proxy() : m_pObj(0), m_PartId(0), m_Flags(eDefault), m_UserFlags(0), m_TimeStamp(0), m_pBV(0) {}
    ~Proxy() {}
    
    //! \name Access to typed void* data
    //@{
    template <typename T> inline T GetObj() const { return reinterpret_cast<T>(m_pObj); }
    //@}

    //! \name Flags helpers
    //@{
    inline bool IsBoundary() const { return m_Flags.Test(eBoundary); }
    inline bool IsStatic() const { return m_Flags.Test(eStatic); }
    inline bool IsNestedBP() const { return m_Flags.Test(eNestedBP); }
    //@}
};

}} //namespace geo::bp

#endif // GEO_BP_PROXY_H
