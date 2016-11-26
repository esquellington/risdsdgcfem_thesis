#include "BasicBP.h"
#include <Geo/bv/TestRay.h>

namespace geo {
namespace bp {

BasicBP::BasicBP( unsigned int max_objects )
: m_poolProxiesDefault(max_objects)
, m_poolProxiesBoundary(max_objects)
{
}

BasicBP::~BasicBP()
{
}

//---- Proxy management
const Proxy *BasicBP::Add( void *p_obj, uint16 part_id,
                           const bv::IBoundingVolume *p_bv,
                           Flags8 proxy_flags,
                           Flags8 user_flags )
{
    GEO_ASSERT( 0 == FindProxy( p_obj, part_id ) );
    Proxy *pProxy = ( proxy_flags.Test(Proxy::eBoundary) ) ? m_poolProxiesBoundary.New() : m_poolProxiesDefault.New();
    pProxy->m_pObj = p_obj;
    pProxy->m_PartId = part_id;
    int32 lAutoFlags = 0;
    pProxy->m_Flags = proxy_flags | lAutoFlags;
    pProxy->m_UserFlags = user_flags;
    pProxy->m_pBV = p_bv;
    pProxy->m_TimeStamp = p_bv->GetTimeStamp();
    return pProxy;
}

// IsValid() already tested in GPoolDA
void BasicBP::Remove( const Proxy *p_proxy )
{
    if( p_proxy->IsBoundary() ) m_poolProxiesBoundary.Delete( const_cast<Proxy*>(p_proxy) );
    else m_poolProxiesDefault.Delete( const_cast<Proxy*>(p_proxy) );
    //\note All pairs involved will vanish on the next Test
}

const Proxy *BasicBP::FindProxy( void *p_obj, uint16 part_id ) const
{
    for( PoolProxy::iterator it_proxy=m_poolProxiesDefault.Begin(); it_proxy.IsValid(); ++it_proxy )
        if( it_proxy->m_pObj == p_obj && it_proxy->m_PartId == part_id )
            return it_proxy.GetPtr();
    for( PoolProxy::iterator it_proxy=m_poolProxiesBoundary.Begin(); it_proxy.IsValid(); ++it_proxy )
        if( it_proxy->m_pObj == p_obj && it_proxy->m_PartId == part_id )
            return it_proxy.GetPtr();
    return 0;
}

void BasicBP::Clear()
{
    m_poolProxiesDefault.Clear();
    m_poolProxiesBoundary.Clear();
}

unsigned int BasicBP::RemoveProxies( Flags32 proxy_flags )
{
    // Iterate over all proxies and remove matching ones
    unsigned int num_removed(0);
    // Default proxies
    PoolProxy::iterator it_proxy=m_poolProxiesDefault.Begin();
    while( it_proxy.IsValid() )
    {
        if( it_proxy->m_Flags.Test(proxy_flags) )
        {
            Proxy *pDeleted(it_proxy.GetPtr());
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
            Proxy *pDeleted(it_proxy.GetPtr());
            ++it_proxy;
            m_poolProxiesBoundary.Delete(pDeleted);
            num_removed++;
        }
        else
            ++it_proxy;
    }
    return num_removed;
}

// Proxy iteration, do not Add/Remove proxies while iterating
const Proxy *BasicBP::FirstProxy() const
{
    m_ProxyIterator = m_poolProxiesDefault.Begin();
    if( !m_ProxyIterator.IsValid() ) m_ProxyIterator = m_poolProxiesBoundary.Begin();
    if( m_ProxyIterator.IsValid() ) return m_ProxyIterator.GetPtr();
    else return 0;
}

const Proxy *BasicBP::NextProxy() const
{
    GEO_ASSERT( m_ProxyIterator.IsValid() );
    ++m_ProxyIterator;
    if( m_ProxyIterator.IsValid() ) return m_ProxyIterator.GetPtr();
    else if( m_ProxyIterator.GetPool() == &m_poolProxiesDefault ) m_ProxyIterator = m_poolProxiesBoundary.Begin();

    if( m_ProxyIterator.IsValid() ) return m_ProxyIterator.GetPtr();
    else return 0;
}

//---- Overlap Queries
void BasicBP::Update()
{
    // Pull per-proxy BV and check timestamp
    /*
    for( PoolProxy::iterator it_proxy=m_poolProxies.Begin(); it_proxy.IsValid(); ++it_proxy )
        if( it_proxy->m_pBV->TestDirtyAndUpdate(it_proxy->m_TimeStamp) )
            Reclassify(it_proxy);
    */
}

unsigned int BasicBP::TestAll( PairContainer &pair_container,
                               Flags32 test_flags,
                               CB_PairFilter cb_pair_filter ) const
{
    // Clear pairs if non-persistent, update (purge/vanish) them otherwise
    if( !test_flags.Test(eTest_Persistent) ) pair_container.Clear();
    else pair_container.Update();

    // Select persistent or non-persistent OrderedPairTestFunction
    OrderedPTF OPTF = (test_flags.Test(eTest_Persistent)) ? OPTF_Persistent : OPTF_NotPersistent;

    //---- Perform O(n^2) tests
    unsigned int num_pairs(0);

    // Internal * Boundary
    if( test_flags.Test(eTest_Boundary) )
        for( PoolProxy::iterator it_int = m_poolProxiesDefault.Begin(); it_int.IsValid(); ++it_int )
            for( PoolProxy::iterator it_bnd = m_poolProxiesBoundary.Begin(); it_bnd.IsValid(); ++it_bnd )
                num_pairs += OPTF( it_int.GetPtr(), it_bnd.GetPtr(), pair_container, cb_pair_filter );

    // Internal * Internal
    if( test_flags.Test(eTest_Internal) )
        for( PoolProxy::iterator it_int1 = m_poolProxiesDefault.Begin(); it_int1.IsValid(); ++it_int1 )
        {
            PoolProxy::iterator it_int2( it_int1 );
            for( ++it_int2; it_int2.IsValid(); ++it_int2 )
                num_pairs += OPTF( it_int1.GetPtr(), it_int2.GetPtr(), pair_container, cb_pair_filter );
        }

    return num_pairs;
}

unsigned int BasicBP::TestSingle( const Proxy *p_proxy,
                                  PairContainer &pair_container,
                                  Flags32 test_flags,
                                  CB_PairFilter cb_pair_filter ) const
{
    GEO_LOG_ERROR("bp::BasicBP::TestSingle not yet implemented");
    return 0;
}

unsigned int BasicBP::TestRay2( const np::Ray2 &ray,
                                std::vector<Proxy*> &vec_hit_proxies,
                                Flags32 test_flags,
                                CB_RayFilter cb_ray_filter ) const
{
    // Internal
    if( test_flags.Test(eTest_Internal) )
        for( PoolProxy::iterator it_int = m_poolProxiesDefault.Begin(); it_int.IsValid(); ++it_int )
            if( it_int->m_pBV->GetDimension() == 2
                && cb_ray_filter( it_int.GetPtr() )
                && bv::TestRay(static_cast<const bv::BoundingVolume2*>(it_int->m_pBV),ray) )
                vec_hit_proxies.push_back( it_int.GetPtr() );
    // Boundary
    if( test_flags.Test(eTest_Boundary) )
        for( PoolProxy::iterator it_bnd = m_poolProxiesBoundary.Begin(); it_bnd.IsValid(); ++it_bnd )
            if( it_bnd->m_pBV->GetDimension() == 2
                && cb_ray_filter( it_bnd.GetPtr() )
                && bv::TestRay(static_cast<const bv::BoundingVolume2*>(it_bnd->m_pBV),ray) )
                vec_hit_proxies.push_back( it_bnd.GetPtr() );
    return vec_hit_proxies.size();
}

unsigned int BasicBP::TestRay3( const np::Ray3 &ray,
                                std::vector<Proxy*> &vec_hit_proxies,
                                Flags32 test_flags,
                                CB_RayFilter cb_ray_filter ) const
{
    // Internal
    if( test_flags.Test(eTest_Internal) )
        for( PoolProxy::iterator it_int = m_poolProxiesDefault.Begin(); it_int.IsValid(); ++it_int )
            if( it_int->m_pBV->GetDimension() == 3
                && cb_ray_filter( it_int.GetPtr() )
                && bv::TestRay(static_cast<const bv::BoundingVolume3*>(it_int->m_pBV),ray) )
                vec_hit_proxies.push_back( it_int.GetPtr() );
    // Boundary
    if( test_flags.Test(eTest_Boundary) )
        for( PoolProxy::iterator it_bnd = m_poolProxiesBoundary.Begin(); it_bnd.IsValid(); ++it_bnd )
            if( it_bnd->m_pBV->GetDimension() == 3
                && cb_ray_filter( it_bnd.GetPtr() )
                && bv::TestRay(static_cast<const bv::BoundingVolume3*>(it_bnd->m_pBV),ray) )
                vec_hit_proxies.push_back( it_bnd.GetPtr() );
    return vec_hit_proxies.size();
}

}} //namespace geo::bp
