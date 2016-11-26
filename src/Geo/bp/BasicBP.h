#ifndef GEO_BP_BASIC_BP_H
#define GEO_BP_BASIC_BP_H

#include "IBroadPhase.h"
#include <util/GPoolDA.h>

namespace geo {
namespace bp {

/*! Simplest broadphase implementation:
  - TestAll is O(n^2)
  - TestSingle is O(n)
*/
class BasicBP: public IBroadPhase
{
public:
    BasicBP( unsigned int max_objects );
    ~BasicBP();

    EBroadPhaseType GetType() const { return eBP_Basic; }

    //!\name Proxy management
    //@{
    const Proxy *Add( void *p_obj, uint16 part_id,
                      const bv::IBoundingVolume *p_bv,
                      Flags8 proxy_flags,
                      Flags8 user_flags );
    void Remove( const Proxy *p_proxy );
    const Proxy *FindProxy( void *p_obj, uint16 part_id = 0 ) const;
    void Clear();
    unsigned int RemoveProxies( Flags32 proxy_flags ); //! Remove proxies matching flags, returns num removed
    // Proxy iteration, do not Add/Remove proxies while iterating
    const Proxy *FirstProxy() const;
    const Proxy *NextProxy() const;
    //@}

    //!\name Overlap Queries
    //@{
    void Update();
    unsigned int TestAll( PairContainer &pair_container,
                          Flags32 test_flags,
                          CB_PairFilter cb_pair_filter = CBPF_True ) const;
    unsigned int TestSingle( const Proxy *p_proxy,
                             PairContainer &pair_container,
                             Flags32 test_flags,
                             CB_PairFilter cb_pair_filter = CBPF_True ) const;
    unsigned int TestRay2( const np::Ray2 &ray,
                           std::vector<Proxy*> &vec_hit_proxies,
                           Flags32 test_flags,
                           CB_RayFilter cb_ray_filter = CBRF_True ) const;
    unsigned int TestRay3( const np::Ray3 &ray,
                           std::vector<Proxy*> &vec_hit_proxies,
                           Flags32 test_flags,
                           CB_RayFilter cb_ray_filter = CBRF_True ) const;
    //@}

private:
    typedef util::GPoolDA<Proxy> PoolProxy;
    PoolProxy m_poolProxiesDefault;
    PoolProxy m_poolProxiesBoundary;
    mutable PoolProxy::iterator m_ProxyIterator;
};

}} //namespace geo::bp

#endif // GEO_BP_BASIC_BP_H
