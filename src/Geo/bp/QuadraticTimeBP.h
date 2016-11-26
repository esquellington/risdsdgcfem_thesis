#ifndef GEO_BP_QUADRATIC_TIME_BP_H
#define GEO_BP_QUADRATIC_TIME_BP_H

#include "IBroadPhase.h"
#include <Geo/bp/GBipartitionTree.h>

namespace geo {
namespace bp {

/*! Simplest broadphase implementation:
  - TestAll is O(n^2)
  - TestSingle is O(n)

  IMPORTANT: This BP should NOT be used, it's only kept as an example
of GBipartitionTree usage. Default/Boundary proxy clustering is much
more efficiently implemented using separate pools, so use BasicBP
instead!
*/
class QuadraticTimeBP: public IBroadPhase
{
public:
    QuadraticTimeBP( unsigned int max_objects );
    ~QuadraticTimeBP();

    EBroadPhaseType GetType() const { return eBP_QuadraticTime; }
    
    //!\name Proxy management
    //@{
    const Proxy *Add( void *p_obj, const bv::IBoundingVolume *p_bv, Flags32 proxy_flags );
    void Remove( const Proxy *p_proxy );
    const Proxy *FindProxy( void *p_obj ) const;
    void Clear();
    unsigned int RemoveProxies( Flags32 proxy_flags ); //! Remove proxies matching flags, returns num removed
    const Proxy *FirstProxy() const { return 0; }
    const Proxy *NextProxy() const { return 0; }
    //@}

    //!\name Overlap Queries
    //@{
    void Update();
    unsigned int TestAll( PairContainer &pair_container,
                          Flags32 test_flags,
                          CB_PairFilter cb_pair_filter = CBPF_True );
    unsigned int TestSingle( const Proxy *p_proxy,
                             PairContainer &pair_container,
                             Flags32 test_flags,
                             CB_PairFilter cb_pair_filter = CBPF_True );
    //@}

private:
    //! \name Internal methods
    //@{o
    void RecomputeClusters();
    //@}

public:
    /*! Cluster of proxies according to their flags.      
      \note It's public because we need to use it the definition of
      tree traversal functors in QuadraticTimeBP.cpp
    */
    struct Cluster
    {
        typedef Proxy * Entry;
        Flags32 m_Flags;
        inline void User_Init( Cluster *p_parent, const Entry *vec_entry, uint32 size, int level, void *user_data )
            {
                ;
            }
        inline bool User_NeedsSplit( const Entry *vec_entry, uint32 size, int level, void *user_data )
            {
                return level < 2;
            }
        inline bool User_Classify( const Entry &entry, int level, void *user_data )
            {
                return entry->m_Flags.Test( 1<<level );
            }
        inline void User_SetChildren( Cluster *p_left, Cluster *p_right, int level, void *user_data )
            {
                if(p_left) p_left->m_Flags = m_Flags;
                if(p_right) p_right->m_Flags = m_Flags;
            }
    };
    typedef GBipartitionTree< Proxy*, Cluster > ProxyBPT;
    
    
private:
    typedef util::GPoolDA<Proxy> PoolProxy;
    PoolProxy m_poolProxies;
    
    //! \name Clustering
    //@{
    bool m_IsDirtyClusters;  //!< If objects added or removed, clusters must be recomputed
    ProxyBPT m_ProxyBPT;
    const Cluster *m_pClusterAll;
    const Cluster *m_pClusterBoundary;
    const Cluster *m_pClusterInternal;
    /* Unnecessary in QTBP
    const Cluster *m_pClusterBnd_Static;
    const Cluster *m_pClusterBnd_Moving;
    const Cluster *m_pClusterInt_Static;
    const Cluster *m_pClusterInt_Moving;
    */
    //@}
};

}} //namespace geo::bp

#endif // GEO_BP_QUADRATIC_TIME_BP_H
