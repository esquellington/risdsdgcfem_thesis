#ifndef GEO_BP_I_BROAD_PHASE_H
#define GEO_BP_I_BROAD_PHASE_H

#include <Geo/Config.h>
#include <Geo/bp/Proxy.h>
#include <Geo/bp/Pair.h>
#include <Geo/bp/GPairContainerHT.h>
#include <Geo/bv/BoundingVolume.h>
#include <Geo/bv/TestOverlap.h> //TEMPORAL: Req by OPTF implementations, but I don't like placing it here...

#include <Geo/np/RayCast.h> //Required for Ray and RayHit types...
#include <vector> //Required for TestRay results

namespace geo { namespace bp
{

//! BP test flags
enum ETestFlags {
    // TestAll()/TestSingle() flags
    eTest_Internal    = (1<<0), //!< Int2Int (unordered) or Single2Int overlaps
    eTest_Boundary    = (1<<1), //!< Int2Boundary (ordered) or Single2Bondary overlaps
    // Common flags
    eTest_NonPersistent = 0,
    eTest_Persistent    = (1<<2), //!< Pairs are kept persistent, not simply added as eNew
    eTest_Discrete      = 0,
    eTest_Continuous    = (1<<3), //!< CCD: Swept BV is tested, not only current BV (req specific Proxy/BV support...)
};


//! BP types
enum EBroadPhaseType {
    eBP_Basic = 0,
    eBP_SweepAndPrune,
    eBP_QuadraticTime //!< deprecated
};

/*! Dimension-Independent Broad Phase interface.

  Usage:
  1) Add/Remove: Proxy lifetime.
  2) Update: If proxies have changed, may rebuild internal structures
  that accelerate overlap detection.
  3) TestAll / TestSingle: Configurable all-pairs or one-versus-all
  tests.

  Overlapping pair persistence is guaranteed when eTest_Persistent is
  enabled in TestAll/TestSingle calls, so that Pair::m_UserData can
  be used to externally "annotate" overlapping pairs with specific
  info (such as a persistent contact manifold). However, the order in
  which each TPairBP appear in the returned array may change each
  Update().

  After each TestAll/TestSingle() call, overlapping pairs can be
  accessed and their individual state checked. Note that when a pair
  is in eVanished state, its objects and userdata MAY NOT EXIST
  ANYMORE (eg: they were removed by the user), thus, as a general
  rule, you cannot use objects from out-of-date overlapping pairs.
*/
class IBroadPhase
{
public:
    typedef bool (*CB_PairFilter)( const Proxy *p_proxy1, const Proxy *p_proxy2 );
    static bool CBPF_True( const Proxy *p_proxy1, const Proxy *p_proxy2 ) { return true; }

    typedef bool (*CB_RayFilter)( const Proxy *p_proxy );
    static bool CBRF_True( const Proxy *p_proxy ) { return true; }

public:
    IBroadPhase() {}
    virtual ~IBroadPhase() {}

    virtual EBroadPhaseType GetType() const = 0;

    //!\name Nested-BP interface
    //@{
    virtual bool IsNested() const { return false; }
    virtual const IBroadPhase *GetParentBP() const { return 0; }
    virtual IBroadPhase *GetParentBP() { return 0; }
    //virtual const bv::IBoundingVolume &GetBV() const = 0;
    //@}

    //!\name Proxy management
    //@{
    virtual const Proxy *Add( void *p_obj, uint16 part_id,
                              const bv::IBoundingVolume *p_bv,
                              Flags8 proxy_flags,
                              Flags8 user_flags ) = 0;
    virtual void Remove( const Proxy *p_proxy ) = 0;
    virtual const Proxy *FindProxy( void *p_obj, uint16 part_id = 0 ) const = 0;
    virtual void Clear() = 0;
    virtual unsigned int RemoveProxies( Flags32 proxy_flags ) = 0; //! Remove proxies matching flags, returns #removed
    // Proxy iteration, do not Add/Remove proxies while iterating
    virtual const Proxy *FirstProxy() const = 0;
    virtual const Proxy *NextProxy() const = 0;
    //@}

    //!\name Overlap Queries
    //@{
    virtual void Update() = 0;
    virtual unsigned int TestAll( PairContainer &pair_container,
                                  Flags32 test_flags,
                                  CB_PairFilter cb_pair_filter = CBPF_True ) const = 0;
    virtual unsigned int TestSingle( const Proxy *p_proxy,
                                     PairContainer &pair_container,
                                     Flags32 test_flags,
                                     CB_PairFilter cb_pair_filter = CBPF_True ) const = 0;
    /*\todo Ugly, raycast needs ray dimension, but all IBroadPhase
     * interface is dimensionless... adding an abstract IRay with
     * generic dimension would be just ugly too, especially for
     * RayHit2... maybe we should just have a SINGLE Ray class (3d)
     * that can work in 2d (using only 2 coords...) */
    virtual unsigned int TestRay2( const np::Ray2 &ray,
                                   std::vector<Proxy*> &vec_hit_proxies,
                                   Flags32 test_flags,
                                   CB_RayFilter cb_ray_filter = CBRF_True ) const = 0;
    /*\todo Ugly, raycast needs ray dimension, but all IBroadPhase
     * interface is dimensionless... adding an abstract IRay with
     * generic dimension would be just ugly too, especially for
     * RayHit2... maybe we should just have a SINGLE Ray class (3d)
     * that can work in 2d (using only 2 coords...) */
    virtual unsigned int TestRay3( const np::Ray3 &ray,
                                   std::vector<Proxy*> &vec_hit_proxies,
                                   Flags32 test_flags,
                                   CB_RayFilter cb_ray_filter = CBRF_True ) const = 0;
    //@}

    //!\name SLOW helper methods
    //@{
    //\todo disabled because a single p_obj can appear with multiple part_id!!
    //inline void Remove( void *p_obj ) { Remove( FindProxy(p_obj) ); }
    //@}
};


//! Prototype of "Ordered Pair Test Function" with flags. Pair (a,b) is added in explicit order.
typedef unsigned int (*OrderedPTF)( Proxy *p_proxy1, Proxy *p_proxy2,
                                    PairContainer &pair_container,
                                    IBroadPhase::CB_PairFilter cb_pair_filter );

//! Test pair with persistence. Retunrs 1 if overlap, 0 if don't.
inline unsigned int OPTF_Persistent( Proxy *p_proxy1, Proxy *p_proxy2,
                                     PairContainer &pair_container,
                                     IBroadPhase::CB_PairFilter cb_pair_filter )
{
    if( cb_pair_filter( p_proxy1, p_proxy2 )
        && bv::TestOverlap( p_proxy1->m_pBV, p_proxy2->m_pBV ) )
    {
        bool bFound(false);
        Pair *pPair( pair_container.FindOrNew( p_proxy1, p_proxy2, bFound ) );
        if( bFound )
            pair_container.Persist( pPair );
        else
            pPair->m_UserData = 0; //Clients req userdata=0 in New pairs
        return 1;
    }
    else
        return 0;
}

//! Test pair without persistence. Retunrs 1 if overlap, 0 if don't.
inline unsigned int OPTF_NotPersistent( Proxy *p_proxy1, Proxy *p_proxy2,
                                        PairContainer &pair_container,
                                        IBroadPhase::CB_PairFilter cb_pair_filter )
{
    if( cb_pair_filter( p_proxy1, p_proxy2 )
        && bv::TestOverlap( p_proxy1->m_pBV, p_proxy2->m_pBV ) )
    {
        Pair *pPair = pair_container.NewNP( p_proxy1, p_proxy2 );
        pPair->m_UserData = 0; //Clients req userdata=0 in New pairs
        return 1;
    }
    else
        return 0;
}

}} //namespace geo::bp

#endif // GEO_BP_I_BROAD_PHASE_H
