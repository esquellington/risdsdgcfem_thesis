#ifndef GEO_BP_SAP_BP_H
#define GEO_BP_SAP_BP_H

#include "IBroadPhase.h"
#include <util/GPoolDA.h>
#include <vector>

namespace geo {
namespace bp {

/*! Sweep and Prune broadphase implementation:
  - Update is O(n log n)
  - TestAll is O(n log n)
  - TestSingle is O(n)

  \sa in-place mergesort: http://en.literateprograms.org/Merge_sort_%28C_Plus_Plus%29
*/
class SweepAndPruneBP: public IBroadPhase
{
public:
    enum EConstants {
        cMaxDirtyProxiesBeforeTotalResort = 0,  //!< \todo partial-sort disabled by now
        cReservedHeapSize = 128
    };

public:
    SweepAndPruneBP( unsigned int max_objects );
    ~SweepAndPruneBP();

    EBroadPhaseType GetType() const { return eBP_SweepAndPrune; }

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
                           CB_RayFilter cb_ray_filter ) const;
    unsigned int TestRay3( const np::Ray3 &ray,
                           std::vector<Proxy*> &vec_hit_proxies,
                           Flags32 test_flags,
                           CB_RayFilter cb_ray_filter ) const;
    //@}

public:
    //! Specific proxy with SaP axis interval, public required
    struct ProxySAP : public Proxy
    {
        Interval m_Interval;
        finline bool operator<( const ProxySAP &proxy_sap ) const { return m_Interval.Min() < proxy_sap.m_Interval.Min(); }
    };

private:
    typedef util::GPoolDA<ProxySAP> PoolProxy;
    PoolProxy m_poolProxiesDefault;
    PoolProxy m_poolProxiesBoundary;
    PoolProxy m_poolProxiesStaticBoundary;   //!< Specific pool to avoid re-sorting them
    mutable PoolProxy::iterator m_ProxyIterator;

    typedef std::vector<ProxySAP*> HeapProxy;
    mutable HeapProxy m_heapOpen1; //mutable because we alloc it once but use it in TestAll() CONST method
    mutable HeapProxy m_heapOpen2; //mutable because we alloc it once but use it in TestAll() CONST method

private:
    void SortProxies( PoolProxy &proxies );
    unsigned int Sweep1( const PoolProxy &proxies,
                         PairContainer &pair_container,
                         CB_PairFilter cb_pair_filter,
                         OrderedPTF OPTF ) const;
    unsigned int Sweep2( const PoolProxy &proxies1, const PoolProxy &proxies2,
                         PairContainer &pair_container,
                         CB_PairFilter cb_pair_filter,
                         OrderedPTF OPTF ) const;
private:
    struct Stats
    {
        finline void Reset() { m_NumPairTests = 0; m_MaxPairTests = 0; m_NumPairOverlaps = 0; m_MaxHeapSize = 0; }
        uint32 m_NumPairTests;
        uint32 m_MaxPairTests; //BF
        uint32 m_NumPairOverlaps;
        uint32 m_MaxHeapSize;
    };
    mutable Stats m_Stats;
};

}} //namespace geo::bp

#endif // GEO_BP_SAP_BP_H
