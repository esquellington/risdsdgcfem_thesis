#include "QuadraticTimeBP.h"

namespace geo {
namespace bp {

QuadraticTimeBP::QuadraticTimeBP( unsigned int max_objects )
: m_poolProxies(max_objects)
, m_ProxyBPT(max_objects,1+2+4) //3-level tree
{
    GEO_LOG_ERROR("geo::bp::QuadraticTimeBP is DEPRECATED, use geo::bp::BasicBP instead!" );
}

QuadraticTimeBP::~QuadraticTimeBP()
{
}

//---- Proxy management
const Proxy *QuadraticTimeBP::Add( void *p_obj, const bv::IBoundingVolume *p_bv, Flags32 proxy_flags )
{
    GEO_ASSERT( 0 == FindProxy(p_obj) );
    int32 lAutoFlags = 0;
    Proxy *pProxy = m_poolProxies.New();
    pProxy->m_pObj = p_obj;
    pProxy->m_Flags = proxy_flags | lAutoFlags;
    pProxy->m_pBV = p_bv;
    pProxy->m_TimeStamp = p_bv->GetTimeStamp();
    m_IsDirtyClusters = true;
    return pProxy;
}

void QuadraticTimeBP::Remove( const Proxy *p_proxy )
{
    GEO_ASSERT( m_poolProxies.IsValid(p_proxy) );
    m_poolProxies.Delete( const_cast<Proxy*>(p_proxy) );
    m_IsDirtyClusters = true;
    //\note All pairs involved will vanish on the next Test
}

const Proxy *QuadraticTimeBP::FindProxy( void *p_obj ) const
{
    for( PoolProxy::iterator it_proxy=m_poolProxies.Begin(); it_proxy.IsValid(); ++it_proxy )
        if( it_proxy->m_pObj == p_obj )
            return it_proxy.GetPtr();
    return 0;
}

void QuadraticTimeBP::Clear()
{
    m_poolProxies.Clear();
    m_IsDirtyClusters = true;
}

unsigned int QuadraticTimeBP::RemoveProxies( Flags32 proxy_flags )
{
    // Iterate over all proxies and remove matching ones
    unsigned int num_removed(0);
    PoolProxy::iterator it_proxy=m_poolProxies.Begin();
    while( it_proxy.IsValid() )
    {
        if( 0 != (it_proxy->m_Flags & proxy_flags) )
        {
            Proxy *pDeleted(it_proxy.GetPtr());
            ++it_proxy;
            m_poolProxies.Delete(pDeleted);
            num_removed++;            
        }
        else
            ++it_proxy;
    }

    // \todo Set completely dirty if any removed, by now...
    if( num_removed > 0 )
        m_IsDirtyClusters = true;
    return num_removed;
}

//---- Overlap Queries
void QuadraticTimeBP::Update()
{
    // Pull per-proxy BV and check timestamp
    /*
    for( PoolProxy::iterator it_proxy=m_poolProxies.Begin(); it_proxy.IsValid(); ++it_proxy )
        if( it_proxy->m_pBV->TestDirtyAndUpdate(it_proxy->m_TimeStamp) )
            Reclassify(it_proxy);
    */
    // Reclassify Boundary/Internal if necessary
    if( m_IsDirtyClusters ) RecomputeClusters();
}

unsigned int QuadraticTimeBP::TestAll( PairContainer &pair_container,
                                       Flags32 test_flags,
                                       CB_PairFilter cb_pair_filter )
{
    if( m_IsDirtyClusters ) RecomputeClusters();

    // Clear pairs if non-persistent, update (purge/vanish) them otherwise
    if( !test_flags.Test(eTest_Persistent) ) pair_container.Clear();
    else pair_container.Update();

    // Select persistent or non-persistent OrderedPairTestFunction
    OrderedPTF OPTF = (test_flags.Test(eTest_Persistent)) ? OPTF_Persistent : OPTF_NotPersistent;
    
    //---- Perform O(n^2) tests
    unsigned int num_pairs(0);

    // Internal * Boundary
    if( test_flags.Test(eTest_Boundary)
        && 0 != m_pClusterInternal
        && 0 != m_pClusterBoundary )
        for( ProxyBPT::entry_iterator it_int = m_ProxyBPT.GetEntryIterator(*m_pClusterInternal);
             it_int.IsValid();
             ++it_int )
            for( ProxyBPT::entry_iterator it_bnd = m_ProxyBPT.GetEntryIterator(*m_pClusterBoundary);
                 it_bnd.IsValid();
                 ++it_bnd )
                num_pairs += OPTF( *it_int, *it_bnd, pair_container, cb_pair_filter );
    
    // Internal * Internal
    if( test_flags.Test(eTest_Internal)
        && 0 != m_pClusterInternal )
        for( ProxyBPT::entry_iterator it_int1 = m_ProxyBPT.GetEntryIterator(*m_pClusterInternal);
             it_int1.IsValid();
             ++it_int1 )
        {
            ProxyBPT::entry_iterator it_int2( it_int1 );
            for( ++it_int2; it_int2.IsValid(); ++it_int2 )
                num_pairs += OPTF( *it_int1, *it_int2, pair_container, cb_pair_filter );
        }
    
    return num_pairs;
}

unsigned int QuadraticTimeBP::TestSingle( const Proxy *p_proxy,
                                          PairContainer &pair_container,
                                          Flags32 test_flags,
                                          CB_PairFilter cb_pair_filter )
{
    return 0;
}

//TEMPORAL: TEST
struct Functor_ComputeSize
{
    uint32 m_AccSize;
    Functor_ComputeSize() : m_AccSize(0) {}
    inline int32 operator()( const QuadraticTimeBP::Cluster &unode,
                             Proxy **vec_entry, uint32 size, int level,
                             void *user_data )
        {
            m_AccSize += size;
            return eBPT_TraverseBoth;
        }
};

//---- Internal methods
void QuadraticTimeBP::RecomputeClusters()
{
    // Resize BPT (num proxies may have increased)
    m_ProxyBPT.Reset( m_poolProxies.Size() );
    // Rebuild BPT
    for( PoolProxy::iterator it_proxy=m_poolProxies.Begin(); it_proxy.IsValid(); ++it_proxy )
        m_ProxyBPT.Add( it_proxy.GetPtr() );
    m_ProxyBPT.Rebuild();
    m_IsDirtyClusters = false;
    
    m_pClusterAll = m_ProxyBPT.GetRoot();
    m_pClusterBoundary = m_ProxyBPT.GetLeft( *m_pClusterAll );
    m_pClusterInternal = m_ProxyBPT.GetRight( *m_pClusterAll );
    
    //TEMPORAL: TEST
    Functor_ComputeSize functor;
    m_ProxyBPT.TraversePreorder( functor );
}

}} //namespace geo::bp
