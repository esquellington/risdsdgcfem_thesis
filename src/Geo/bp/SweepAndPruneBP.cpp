#include "SweepAndPruneBP.h"
#include <Geo/np/Projection.h>
#include <Geo/bv/bv.h>

// STL includes required by heap
#include <vector>
#include <algorithm>
#include <functional>

#define __USE_SAP_SWEEP1 //!< If disabled, we use O(n^2) brute-force algorithm
#define __USE_SAP_SWEEP2 //!< If disabled, we use O(n^2) brute-force algorithm

//define __ENABLE_SORT_TRACKING
//#define __ENABLE_EVENT_TRACKING
//#define __ENABLE_STATS

#ifdef __ENABLE_STATS
#  define GEO_STATS_IF(x) { x; }
#  define GEO_STATS_SET(x,y) { x = (y); }
#  define GEO_STATS_INC(x) { x++; }
#  define GEO_STATS_DEC(x) { x--; }
#  define GEO_STATS_ADD(x,y) { x += (y); }
#else
#  define GEO_STATS_IF(x)
#  define GEO_STATS_SET(x,y)
#  define GEO_STATS_INC(x)
#  define GEO_STATS_DEC(x)
#  define GEO_STATS_ADD(x,y)
#endif

#if defined(__ENABLE_SORT_TRACKING) || defined(__ENABLE_EVENT_TRACKING) || defined(__ENABLE_STATS)
#  define __ENABLE_TRACKING
#  include <iostream> //TEMPORAL!!
#endif

namespace geo {
namespace bp {

SweepAndPruneBP::SweepAndPruneBP( unsigned int max_objects )
: m_poolProxiesDefault(max_objects)
, m_poolProxiesBoundary(max_objects)
, m_poolProxiesStaticBoundary(max_objects)
{
    m_heapOpen1.reserve( cReservedHeapSize );
    GEO_STATS_IF( m_Stats.Reset() );
}

SweepAndPruneBP::~SweepAndPruneBP()
{
}

//---- Proxy management
const Proxy *SweepAndPruneBP::Add( void *p_obj, uint16 part_id,
                                   const bv::IBoundingVolume *p_bv,
                                   Flags8 proxy_flags,
                                   Flags8 user_flags )
{
    GEO_ASSERT( 0 == FindProxy( p_obj, part_id ) );
    ProxySAP *pProxy = ( proxy_flags.Test(Proxy::eBoundary) )
                       ? ( proxy_flags.Test(Proxy::eStatic) ) ? m_poolProxiesStaticBoundary.New() : m_poolProxiesBoundary.New()
                       : m_poolProxiesDefault.New();
    pProxy->m_pObj = p_obj;
    pProxy->m_PartId = part_id;
    int32 lAutoFlags = 0;
    pProxy->m_Flags = proxy_flags | lAutoFlags;
    pProxy->m_UserFlags = user_flags;
    pProxy->m_pBV = p_bv;
    if( geo::bv::eBV_LSS2 != p_bv->GetType() ) GEO_LOG_ERROR("BV type %d unsupported", p_bv->GetType() );
    pProxy->m_TimeStamp = 0; //\note On addition it's out-of-wync wrt p_bv->GetTimeStamp() to enforce Sort() in the next Update()
    return pProxy;
}

// IsValid() already tested in GPoolDA
void SweepAndPruneBP::Remove( const Proxy *p_proxy )
{
    ProxySAP *pProxy( static_cast<ProxySAP*>( const_cast<Proxy*>(p_proxy) ) );
    if( p_proxy->IsBoundary() )
    {
        if( p_proxy->IsStatic() ) m_poolProxiesStaticBoundary.Delete( pProxy );
        else m_poolProxiesBoundary.Delete( pProxy);
    }
    else
        m_poolProxiesDefault.Delete( pProxy );
    //\note All pairs involved will vanish on the next Test
}

const Proxy *SweepAndPruneBP::FindProxy( void *p_obj, uint16 part_id ) const
{
    for( PoolProxy::iterator it_proxy=m_poolProxiesDefault.Begin(); it_proxy.IsValid(); ++it_proxy )
        if( it_proxy->m_pObj == p_obj && it_proxy->m_PartId == part_id )
            return it_proxy.GetPtr();
    for( PoolProxy::iterator it_proxy=m_poolProxiesBoundary.Begin(); it_proxy.IsValid(); ++it_proxy )
        if( it_proxy->m_pObj == p_obj && it_proxy->m_PartId == part_id )
            return it_proxy.GetPtr();
    for( PoolProxy::iterator it_proxy=m_poolProxiesStaticBoundary.Begin(); it_proxy.IsValid(); ++it_proxy )
        if( it_proxy->m_pObj == p_obj && it_proxy->m_PartId == part_id )
            return it_proxy.GetPtr();
    return 0;
}

void SweepAndPruneBP::Clear()
{
    m_poolProxiesDefault.Clear();
    m_poolProxiesBoundary.Clear();
    m_poolProxiesStaticBoundary.Clear();
}

unsigned int SweepAndPruneBP::RemoveProxies( Flags32 proxy_flags )
{
    // Iterate over all proxies and remove matching ones
    unsigned int num_removed(0);
    // Default proxies
    PoolProxy::iterator it_proxy=m_poolProxiesDefault.Begin();
    while( it_proxy.IsValid() )
    {
        if( it_proxy->m_Flags.Test(proxy_flags) )
        {
            ProxySAP *pDeleted(it_proxy.GetPtr());
            ++it_proxy;
            m_poolProxiesDefault.Delete(pDeleted);
            num_removed++;
        }
        else
            ++it_proxy;
    }
    // Boundary proxies
    it_proxy=m_poolProxiesBoundary.Begin();
    while( it_proxy.IsValid() )
    {
        if( it_proxy->m_Flags.Test(proxy_flags) )
        {
            ProxySAP *pDeleted(it_proxy.GetPtr());
            ++it_proxy;
            m_poolProxiesBoundary.Delete(pDeleted);
            num_removed++;
        }
        else
            ++it_proxy;
    }
    // Static Boundary proxies
    it_proxy=m_poolProxiesStaticBoundary.Begin();
    while( it_proxy.IsValid() )
    {
        if( it_proxy->m_Flags.Test(proxy_flags) )
        {
            ProxySAP *pDeleted(it_proxy.GetPtr());
            ++it_proxy;
            m_poolProxiesStaticBoundary.Delete(pDeleted);
            num_removed++;
        }
        else
            ++it_proxy;
    }
    return num_removed;
}

// Proxy iteration, do not Add/Remove proxies while iterating
const Proxy *SweepAndPruneBP::FirstProxy() const
{
    m_ProxyIterator = m_poolProxiesDefault.Begin();
    if( !m_ProxyIterator.IsValid() ) m_ProxyIterator = m_poolProxiesBoundary.Begin();
    if( !m_ProxyIterator.IsValid() ) m_ProxyIterator = m_poolProxiesStaticBoundary.Begin();
    if( m_ProxyIterator.IsValid() ) return m_ProxyIterator.GetPtr();
    else return 0;
}

const Proxy *SweepAndPruneBP::NextProxy() const
{
    GEO_ASSERT( m_ProxyIterator.IsValid() );
    ++m_ProxyIterator;
    if( m_ProxyIterator.IsValid() ) return m_ProxyIterator.GetPtr();
    else if( m_ProxyIterator.GetPool() == &m_poolProxiesDefault ) m_ProxyIterator = m_poolProxiesBoundary.Begin();
    else if( m_ProxyIterator.GetPool() == &m_poolProxiesBoundary ) m_ProxyIterator = m_poolProxiesStaticBoundary.Begin();

    if( m_ProxyIterator.IsValid() ) return m_ProxyIterator.GetPtr();
    else return 0;
}

//---- Overlap Queries

/*! \todo Allow option to assume always ALL default proxies dirty, thus
  directly re-sorting them without previous check.
*/
void SweepAndPruneBP::Update()
{
#ifdef __USE_SAP_SWEEP1
    SortProxies( m_poolProxiesDefault );
#endif
#ifdef __USE_SAP_SWEEP2
    SortProxies( m_poolProxiesBoundary );
#endif
#ifdef  __ENABLE_SORT_TRACKING
    // Do nothing with staticboundary, maybe assert that they haven't changed
    for( PoolProxy::iterator it_proxy = m_poolProxiesStaticBoundary.Begin(); it_proxy.IsValid(); ++it_proxy )
        GEO_ASSERT( !it_proxy->m_pBV->IsDirty(it_proxy->m_TimeStamp) );
#endif
}

/*! Strict Weak Ordering binary predicate to be used with
  heap_open. As we want the heap to report the *smallest*
  element, we invert the comparison.
*/
class HeapProxy_Comparison: public std::binary_function< const SweepAndPruneBP::ProxySAP *, const SweepAndPruneBP::ProxySAP*, bool >
{
public:
    HeapProxy_Comparison() {}
    inline bool operator()( const SweepAndPruneBP::ProxySAP *p_proxy1, const SweepAndPruneBP::ProxySAP *p_proxy2 )
        {
            return p_proxy1->m_Interval.Max() > p_proxy2->m_Interval.Max();
        }
};

/*\todo Iterate over default, boundary and static proxies, keeping a
  stack-based list of end-event sorted open proxies. Iterate in
  sweep-axis-order, advancing to the first event across all 4 open
  lists at each step.

  \todo Try to AVOID handling all coincident event-type permutations
  specifically, just ARRANGE comparisons so that a specific "priority"
  is considered... ex: "in case of coincidence of begin and end, always
  report overlap".

  \todo We could SaP along N > 1 axis in several passes using the
  persistent pair cache to perform the accumulated AND of all axis
  overlaps:
  - pc.update()
  - Save vanished from last bp::update()
  - axis1
  - pc.update()
  - axis2
  ...
  - Restore vanished from last bp::update()

  \note std::heap methods from <algorithm> are a bit confusing:
  - A vector<> structure is required.
  - New heap elements are inserted in 2 steps: vector.push_back(), push_heap(begin,end,comparison)
  - The top element is vector.first()
  - The top element is removed in 2 steps: pop_heap(begin,end,comparison), vector.pop_back()
*/
unsigned int SweepAndPruneBP::TestAll( PairContainer &pair_container,
                                       Flags32 test_flags,
                                       CB_PairFilter cb_pair_filter ) const
{
    // Clear pairs if non-persistent, update (purge/vanish) them otherwise
    if( !test_flags.Test(eTest_Persistent) ) pair_container.Clear();
    else pair_container.Update();

    // Select persistent or non-persistent OrderedPairTestFunction
    OrderedPTF OPTF = (test_flags.Test(eTest_Persistent)) ? OPTF_Persistent : OPTF_NotPersistent;

    unsigned int num_overlaps(0);

    GEO_STATS_IF( m_Stats.Reset() );

#ifdef __ENABLE_TRACKING
    static int sCount(0);
    if( ++sCount > 100 ) sCount = 0;
#endif

    // Internal * (Boundary + StaticBoundary)
    if( test_flags.Test(eTest_Boundary) )
    {
        num_overlaps += Sweep2( m_poolProxiesDefault, m_poolProxiesBoundary, pair_container, cb_pair_filter, OPTF );
        num_overlaps += Sweep2( m_poolProxiesDefault, m_poolProxiesStaticBoundary, pair_container, cb_pair_filter, OPTF );
    }

    // Internal * Internal
    if( test_flags.Test(eTest_Internal) )
        num_overlaps += Sweep1( m_poolProxiesDefault, pair_container, cb_pair_filter, OPTF );

    GEO_STATS_SET( m_Stats.m_NumPairOverlaps, num_overlaps );

#ifdef __ENABLE_STATS
    if( sCount == 0 ) std::cout << "Stats:" << std::endl
                                << "#Proxies = " << m_poolProxiesDefault.Size()
                                << ", " << m_poolProxiesBoundary.Size()
                                << ", " << m_poolProxiesStaticBoundary.Size() << std::endl
                                << "#Tests = " << m_Stats.m_NumPairTests << std::endl
                                << "#Overlaps = " << m_Stats.m_NumPairOverlaps << std::endl
                                << "#Open = " << m_Stats.m_MaxHeapSize << std::endl
                                << "%Useful = " << ((m_Stats.m_NumPairTests>0)
                                                    ? float(m_Stats.m_NumPairOverlaps)/m_Stats.m_NumPairTests
                                                    : -1.0f) << std::endl
                                << "%SaP/Basic = " << float(m_Stats.m_NumPairTests) / m_Stats.m_MaxPairTests << std::endl;
#endif

    return num_overlaps;
}

unsigned int SweepAndPruneBP::TestSingle( const Proxy *p_proxy,
                                          PairContainer &pair_container,
                                          Flags32 test_flags,
                                          CB_PairFilter cb_pair_filter ) const
{
    GEO_LOG_ERROR("bp::SweepAndPruneBP::TestSingle not yet implemented");
    return 0;
}

unsigned int SweepAndPruneBP::TestRay2( const np::Ray2 &ray,
                                        std::vector<Proxy*> &vec_hit_proxies,
                                        Flags32 test_flags,
                                        CB_RayFilter cb_ray_filter ) const
{
    GEO_LOG_WARNING("bp::SweepAndPruneBP::TestRay2 not yet implemented" );
    return 0;
    /* Dummy implementation that does NOT use S&P, enable if required
    np::RayHit2 rh;
    // Internal
    if( test_flags.Test(eTest_Internal) )
        for( PoolProxy::iterator it_int = m_poolProxiesDefault.Begin(); it_int.IsValid(); ++it_int )
            if( it_int->m_pBV->GetDimension() == 2 && cb_ray_filter( it_int.GetPtr() ) && bv::TestRay2(it_int->m_pBV,ray,rh) )
                vec_hit_proxies.push_back( it_int.GetPtr() );
    // Boundary
    if( test_flags.Test(eTest_Boundary) )
        for( PoolProxy::iterator it_bnd = m_poolProxiesBoundary.Begin(); it_bnd.IsValid(); ++it_bnd )
            if( it_bnd->m_pBV->GetDimension() == 2 && cb_ray_filter( it_bnd.GetPtr() ) && bv::TestRay2(it_bnd->m_pBV,ray,rh) )
                vec_hit_proxies.push_back( it_bnd.GetPtr() );
    // Static Boundary
    if( test_flags.Test(eTest_Boundary) )
        for( PoolProxy::iterator it_bnd = m_poolProxiesStaticBoundary.Begin(); it_bnd.IsValid(); ++it_bnd )
            if( it_bnd->m_pBV->GetDimension() == 2 && cb_ray_filter( it_bnd.GetPtr() ) && bv::TestRay2(it_bnd->m_pBV,ray,rh) )
                vec_hit_proxies.push_back( it_bnd.GetPtr() );
    return vec_hit_proxies.size();
    */
}

unsigned int SweepAndPruneBP::TestRay3( const np::Ray3 &ray,
                                        std::vector<Proxy*> &vec_hit_proxies,
                                        Flags32 test_flags,
                                        CB_RayFilter cb_ray_filter ) const
{
    GEO_LOG_WARNING("bp::SweepAndPruneBP::TestRay2 not yet implemented" );
    return 0;
}

//---- Internal

/*! For a given proxy pool, pull per-proxy BV and check timestamp. If
  it's out-of date, save the proxy in a temporary list. If there's a
  few changed proxies, re-insert them, otherwise, perform a complete
  sort of the pool.
*/
void SweepAndPruneBP::SortProxies( PoolProxy &proxies )
{
    const Proxy *list_dirty_proxies[cMaxDirtyProxiesBeforeTotalResort];
    int num_dirty_proxies(0);
    for( PoolProxy::iterator it_proxy = proxies.Begin(); it_proxy.IsValid(); ++it_proxy )
        if( it_proxy->m_pBV->TestDirtyAndUpdate(it_proxy->m_TimeStamp) )
        {
            if( num_dirty_proxies < cMaxDirtyProxiesBeforeTotalResort )
                list_dirty_proxies[num_dirty_proxies] = it_proxy.GetPtr(); // Save in the to-insert list

            //TEMPORAL!!
            if( geo::bv::eBV_LSS2 != it_proxy->m_pBV->GetType() )
                GEO_LOG_ERROR("BV type %d unsupported", it_proxy->m_pBV->GetType() );

            //\todo bv-type-polimorphic Projection with generic dimension too
            const bv::LSS2 &lss2( it_proxy->m_pBV->As<bv::LSS2>() );
            np::GProjection_LSS<2,0>( lss2.GetPos0(), lss2.GetPos1(), lss2.GetRadius(), it_proxy->m_Interval );
            num_dirty_proxies++;
        }
    if( num_dirty_proxies > 0 )//cMaxDirtyProxiesBeforeTotalResort )
        proxies.Sort();
    /*
    else
    {
        GEO_LOG_ERROR("SweepAndPruneBP::SortProxies() Partial sorting not yet implemented");
        GEO_ASSERT(false);
        ;//reinsert. It can be efficiently done if the elements were
         //removed from the sorted pool list, sorted on their own, and
         //then merged into the sorted list again.
    }
    */
#ifdef __ENABLE_SORT_TRACKING
    // Log and Check correct sort
    static int sCount(0);
    if( ++sCount > 100 )
    {
        std::cout << "<SaP> " << std::endl;
        std::cout << "%Dirty = " << float(num_dirty_proxies)/proxies.Size() << std::endl;
        Real last_elem = proxies.IsEmpty() ? Real(0) : proxies.Begin()->m_Interval.Min();
        for( PoolProxy::iterator it_proxy = proxies.Begin(); it_proxy.IsValid(); ++it_proxy )
        {
            Real elem( it_proxy->m_Interval.Min() );
            if( last_elem > elem )
                std::cout << "SHIIIT " << last_elem << " > " << elem << std::endl;
            GEO_ASSERT( last_elem <= elem );
            last_elem = elem;
            std::cout << elem << std::endl;
        }
        std::cout << "<\\SaP>" << std::endl;
        sCount = 0;
    }
#endif
}

unsigned int SweepAndPruneBP::Sweep1( const PoolProxy &proxies,
                                      PairContainer &pair_container,
                                      CB_PairFilter cb_pair_filter,
                                      OrderedPTF OPTF ) const
{
    if( proxies.Size() < 2 ) return 0;
    unsigned int num_overlaps(0);

#ifdef __USE_SAP_SWEEP1
    m_heapOpen1.clear();
    PoolProxy::iterator it_opening( proxies.Begin() );
    m_heapOpen1.push_back( it_opening.GetPtr() );
    std::push_heap( m_heapOpen1.begin(), m_heapOpen1.end(), HeapProxy_Comparison() );
    Real next_closing( it_opening->m_Interval.Max() );
    ++it_opening;
    Real next_opening( it_opening->m_Interval.Min() );
#ifdef __ENABLE_EVENT_TRACKING
    if( sCount == 0 )
    {
        std::cout << "<i2i>" << std::endl;
        std::cout << "First Opening [ " << (int)proxies.Begin().GetPtr() << " ] = "
                  << proxies.Begin()->m_Interval.Min() << std::endl;
        std::cout << "First Events O = " << next_opening << ", C = " << next_closing << " ]" << std::endl;
        sCount = 0;
    }
#endif
    while( it_opening.IsValid() || !m_heapOpen1.empty() )
    {
        // In case of opening/closing coincidence (=), we
        // prioritize opening in order to consider the just-open
        // proxy(es) when closing the coincident proxy(es).
        if( it_opening.IsValid() && next_opening <= next_closing )
        {
#ifdef __ENABLE_EVENT_TRACKING
            if( sCount == 0 ) std::cout << "Opening [ " << (int)it_opening.GetPtr() << " ] = " << next_opening << std::endl;
#endif
            // Add proxy to open heap
            m_heapOpen1.push_back( it_opening.GetPtr() );
            std::push_heap( m_heapOpen1.begin(), m_heapOpen1.end(), HeapProxy_Comparison() );
            // Advance to next opening event if any
            ++it_opening;
            if( it_opening.IsValid() )
            {
                next_opening = it_opening->m_Interval.Min();
                next_closing = m_heapOpen1.front()->m_Interval.Max();
            }
            else next_opening = 1000001.0f;
        }
        else if( !m_heapOpen1.empty() && next_opening > next_closing )
        {
#ifdef __ENABLE_EVENT_TRACKING
            if( sCount == 0 ) std::cout << "Closing [ " << (int)m_heapOpen1.front() << " ] = " << next_closing << std::endl;
#endif
            // Test closing against all open, including other
            // proxies potentially closing at the same axis max
            // value.
            ProxySAP *pClosing = m_heapOpen1.front();
            std::pop_heap( m_heapOpen1.begin(), m_heapOpen1.end(), HeapProxy_Comparison() );
            m_heapOpen1.pop_back();
            if( !m_heapOpen1.empty() )
            {
                for( HeapProxy::const_iterator it_ho=m_heapOpen1.begin(); it_ho!=m_heapOpen1.end(); ++it_ho )
                {
#ifdef __ENABLE_EVENT_TRACKING
                    if( sCount == 0 ) std::cout << "Testing ( " << (int)pClosing << ", " << (int)*it_ho << " )"
                                                << " [ " << pClosing->m_Interval.Min() << "," << pClosing->m_Interval.Max()
                                                << "] vs [" << (*it_ho)->m_Interval.Min() << "," << (*it_ho)->m_Interval.Max() << "]"
                                                << std::endl;
#endif
                    num_overlaps += OPTF( pClosing, *it_ho, pair_container, cb_pair_filter );
                    GEO_STATS_INC( m_Stats.m_NumPairTests );
                }
                // Advance to next closing event
                next_closing = m_heapOpen1.front()->m_Interval.Max();
            }
            else
                next_closing = 1000000.0f;
        }
        GEO_STATS_SET( m_Stats.m_MaxHeapSize, mal::Max( m_Stats.m_MaxHeapSize, m_heapOpen1.size() ) );
    }

#else //__USE_SAP_SWEEP1

    //---- Perform O(n^2) tests
    for( PoolProxy::iterator it1 = proxies.Begin(); it1.IsValid(); ++it1 )
    {
        PoolProxy::iterator it2( it1 );
        for( ++it2; it2.IsValid(); ++it2 )
            num_overlaps += OPTF( it1.GetPtr(), it2.GetPtr(), pair_container, cb_pair_filter );
    }
    GEO_STATS_ADD( m_Stats.m_NumPairTests, proxies.Size() * (proxies.Size()-1) / 2 );

#endif //__USE_SAP_SWEEP1

    GEO_STATS_ADD( m_Stats.m_MaxPairTests, proxies.Size() * (proxies.Size()-1) / 2 );

    return num_overlaps;
}

unsigned int SweepAndPruneBP::Sweep2( const PoolProxy &proxies1, const PoolProxy &proxies2,
                                      PairContainer &pair_container,
                                      CB_PairFilter cb_pair_filter,
                                      OrderedPTF OPTF ) const
{
    if( proxies1.Size() == 0 || proxies2.Size() == 0 ) return 0;
    unsigned int num_overlaps(0);

#ifdef __USE_SAP_SWEEP2

    m_heapOpen1.clear();
    m_heapOpen2.clear();
    PoolProxy::iterator it_opening1( proxies1.Begin() );
    PoolProxy::iterator it_opening2( proxies2.Begin() );
    if( it_opening1->m_Interval.Min() < it_opening2->m_Interval.Min() )
    {
        m_heapOpen1.push_back( it_opening1.GetPtr() );
        //std::push_heap( m_heapOpen1.begin(), m_heapOpen1.end(), HeapProxy_Comparison() ); unnecessary, no need to sort single element
        ++it_opening1;
    }
    else
    {
        m_heapOpen2.push_back( it_opening2.GetPtr() );
        //std::push_heap( m_heapOpen2.begin(), m_heapOpen2.end(), HeapProxy_Comparison() ); unnecessary, no need to sort single element
        ++it_opening2;
    }
    // Continue while there's Pool1 proxies to open or to close
    while( it_opening1.IsValid() || !m_heapOpen1.empty() )
    {
        if( it_opening1.IsValid()
            && (!it_opening2.IsValid() || it_opening1->m_Interval.Min() < it_opening2->m_Interval.Min())
            && (m_heapOpen1.empty() || it_opening1->m_Interval.Min() <= m_heapOpen1.back()->m_Interval.Max())
            && (m_heapOpen2.empty() || it_opening1->m_Interval.Min() <= m_heapOpen2.back()->m_Interval.Max()) )
        {
            // Open next 1
            m_heapOpen1.push_back( it_opening1.GetPtr() );
            std::push_heap( m_heapOpen1.begin(), m_heapOpen1.end(), HeapProxy_Comparison() );
            ++it_opening1;
        }
        else if( it_opening2.IsValid()
                 //&& (!it_opening2.IsValid() || it_opening2->Min() < it_opening1->Min())
                 && (m_heapOpen1.empty() || it_opening2->m_Interval.Min() <= m_heapOpen1.back()->m_Interval.Max())
                 && (m_heapOpen2.empty() || it_opening2->m_Interval.Min() <= m_heapOpen2.back()->m_Interval.Max()) )
        {
            // Open next 2
            m_heapOpen2.push_back( it_opening2.GetPtr() );
            std::push_heap( m_heapOpen2.begin(), m_heapOpen2.end(), HeapProxy_Comparison() );
            ++it_opening2;
        }
        else if( !m_heapOpen1.empty()
                 && (m_heapOpen2.empty() || m_heapOpen1.back()->m_Interval.Max() <= m_heapOpen2.back()->m_Interval.Max() ) )
        {
            // Close 1, test against all open 2
            ProxySAP *pClosing1( m_heapOpen1.back() );
            m_heapOpen1.pop_back();
            if( !m_heapOpen2.empty() )
            {
                for( HeapProxy::const_iterator it2=m_heapOpen2.begin(); it2 != m_heapOpen2.end(); it2++ )
                    num_overlaps += OPTF( *it2, pClosing1, pair_container, cb_pair_filter );
                GEO_STATS_ADD( m_Stats.m_NumPairTests, m_heapOpen2.size() );
            }
        }
        else if( !m_heapOpen2.empty() )
            //&& (m_heapOpen1.empty() || m_heapOpen2.head()->Max() <= m_heapOpen1.head()->Max() )
        {
            // Close 2, test against all open 1
            ProxySAP *pClosing2( m_heapOpen2.back() );
            m_heapOpen2.pop_back();
            if( !m_heapOpen1.empty() )
            {
                for( HeapProxy::const_iterator it1=m_heapOpen1.begin(); it1 != m_heapOpen1.end(); it1++ )
                    num_overlaps += OPTF( pClosing2, *it1, pair_container, cb_pair_filter );
                GEO_STATS_ADD( m_Stats.m_NumPairTests, m_heapOpen1.size() );
            }
        }
        else
            GEO_LOG_ERROR( "SweepAndPruneBP: This should never happen, no progress has been made, infinite loop!" );
    }

#else //__USE_SAP_SWEEP2
    for( PoolProxy::iterator it1 = proxies1.Begin(); it1.IsValid(); ++it1 )
        for( PoolProxy::iterator it2 = proxies2.Begin(); it2.IsValid(); ++it2 )
            num_overlaps += OPTF( it1.GetPtr(), it2.GetPtr(), pair_container, cb_pair_filter );
    GEO_STATS_ADD( m_Stats.m_NumPairTests, proxies1.Size() * proxies2.Size() );
#endif //__USE_SAP_SWEEP2

    GEO_STATS_ADD( m_Stats.m_MaxPairTests, proxies1.Size() * proxies2.Size() );

    return num_overlaps;
}

}} //namespace geo::bp
