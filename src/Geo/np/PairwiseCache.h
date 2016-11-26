#ifndef GEO_NP_PAIRWISE_CACHE_H
#define GEO_NP_PAIRWISE_CACHE_H

#include <Geo/Config.h>

namespace geo {
namespace np {

/* Generic-Dimension Pairwise-Algorithm-Specific cache (includes
 var-sized caches, such as stochastic deformable CD... allocated and
 managed by EACH Contact_X_Y method requiring it.
*/
template <unsigned D>
class GISpecificPairwiseCache
{
public:
    GISpecificPairwiseCache() {}
    virtual ~GISpecificPairwiseCache()
        {
            /*\todo Delete MUST work and clear any ref regardless of
              where it's called. Specifically, ISpecificPairwiseCache
              will probably be allocated in typed pools, but delete
              will be called when deallocating it's SINGLE referencing
              PairwiseCache, thus, delete must remove "this" from the
              pool. If pools are global static members of each
              ISpecificPairwiseCache subclass, this is easily
              achievable.
            */
        }

    /*TEMP: Not really necessary by now... could be useful if some
      common method must handle different ISPC types, such as
      VizPairwiseCache(), for example
    enum EType { eSPC_Stochastic };
    virtual EType GetType() const = 0;
    */
};

typedef GISpecificPairwiseCache<2> ISpecificPairwiseCache2;
typedef GISpecificPairwiseCache<3> ISpecificPairwiseCache3;

/* Generic-Dimension Pairwise Cache

   The actual data stored in the cache depends on the exact
   Contact_X_Y method used, and can even change for the same pair (ex:
   In Convex Polyhedron pair tests, cache can be a Separating Plane
   from a face or a from an edge, which are different cache types.

   For common and fixed-size cached data, the basic and optional
   PairwiseCache fields are enough. Most Contact_X_Y methods for rigid
   shapes fall into this category. On the other hand, for deformable
   shape contact algorithms cached data can grow arbitrarily as
   contact complexity is unbounded.

   Pair-specific and/or variable-sized data that cannot be statically
   pre-allocated will be stored as a virtual ISpecificPairwiseCache,
   allocated and updated by each Contact_X_Y() that requires it.

   \note EACH Contact_X_Y method MUST only handle the cache types it may
   generate, and ignore the rest...
*/
template <unsigned D>
class GPairwiseCache
{
public:
    finline GPairwiseCache() : m_bIsValid(false), m_pSPC(0) {}
    ~GPairwiseCache() { Invalidate(); }

    finline bool IsValid() const { return m_bIsValid; }
    finline void Invalidate() { m_bIsValid = false; if(m_pSPC) delete m_pSPC; }
    finline GISpecificPairwiseCache<D> *GetSPC() const { return m_pSPC; }
    finline void SetSPC( GISpecificPairwiseCache<D> *p_spc ) { if(m_pSPC) delete m_pSPC; m_pSPC = p_spc; }

private:
    bool m_bIsValid;
    GISpecificPairwiseCache<D> *m_pSPC;
};

typedef GPairwiseCache<2> PairwiseCache2;
typedef GPairwiseCache<3> PairwiseCache3;

/* Convex-Convex SPC.

   Supports different convex2convex cached data types (separating
   plane, witness, etc...) enriched with optional fields (distance,
   transforms...)

   \todo Plane and FeaturePair data might be share the same memory if
   we used a boost variant of Plane and FeaturePair... A union would
   not work I'm afraid.

   \todo We have to be VERY CAREFUL to avoid interference in shared
   cached stuff across Overlap/Contact/Proximity queries between the
   same objects!!
*/
class Convex2ConvexSPC2: public ISpecificPairwiseCache2
{
public:
    enum EFlags {
        // Global
        ePC_IsOverlapping                = (1<<0), //!< Whether objects were overlapping when last tested
        // Separating Plane witness
        ePC_Has_Plane_Explicit           = (1<<1), //!< Explicit plane geometry (N,D) in Shape2 local coords
        ePC_Has_Plane_FeaturePair        = (1<<2), //!< Feature pair (FID1,FID2) yields plane witness
        // Closes points witness
        ePC_Has_Closest_FeaturePair      = (1<<3), //!< Feature pair (FID1,FID2) that yields closest points
        // Min Depth Axis witness
        ePC_Has_MinDepthAxis_FeaturePair = (1<<4), //!< Feature (FID1,FID2) that yields min depth axis
        // Distance witness
        ePC_Has_Distance                 = (1<<5), //!< Separating or penetration (negative) distance
        // Old transforms
        ePC_Has_Transform12              = (1<<6), //!< Useful to infer bounds on distance/depth change, useless for shapes with DOF, though...
    };

    Convex2ConvexSPC2() { m_Flags.Reset(); }
    ~Convex2ConvexSPC2()
        {
            //\todo DEALLOC from pool
        }

    finline Flags32 GetFlags() { return m_Flags; }

    finline void SetOverlaping( bool b_overlapping ) { if(b_overlapping) m_Flags.Enable(ePC_IsOverlapping); else m_Flags.Disable(ePC_IsOverlapping); }

    finline void SetPlane( const Vec2 &normal, Real coeff_d ) { m_Flags.Enable(ePC_Has_Plane_Explicit); m_Normal = normal; m_CoeffD = coeff_d; }
    finline void SetPlane( feature_id fid1, feature_id fid2 ) { m_Flags.Enable(ePC_Has_Plane_FeaturePair); m_FID1 = fid1; m_FID2 = fid2; }
    finline void SetClosest( feature_id fid1, feature_id fid2 ) { m_Flags.Enable(ePC_Has_Closest_FeaturePair); m_FID1 = fid1; m_FID2 = fid2; }
    finline void SetMinDepthAxis( feature_id fid1, feature_id fid2 ) { m_Flags.Enable(ePC_Has_MinDepthAxis_FeaturePair); m_FID1 = fid1; m_FID2 = fid2; }

    finline void SetDistance( Real distance ) { m_Flags.Enable(ePC_Has_Distance); m_Distance = distance; }

    finline void SetTransform12( const Transform2 &tr12 ) { m_Flags.Enable(ePC_Has_Transform12); m_Tr12 = tr12; }

public:
    Flags32 m_Flags;
    // Explicit Plane (\todo Could be struct GPlane<2>)
    Vec2 m_Normal;
    Real m_CoeffD;
    // Feature Pair (\todo Could be struct FeaturePair)
    feature_id m_FID1;
    feature_id m_FID2;
    // Distance/Depth
    Real m_Distance;
    // Relative transform 1 wrt 2 (Tr_1^2)
    Transform2 m_Tr12;
//\todo Alloc from global static POOL
};

//\todo By now, all pairwise test caches are the same basic one, but eventually they may extend it if required
template <unsigned D> struct GContactCache : public GPairwiseCache<D> {};
typedef GContactCache<2> ContactCache2;
typedef GContactCache<3> ContactCache3;

template <unsigned D> struct GOverlapCache : public GPairwiseCache<D> {};
typedef GOverlapCache<2> OverlapCache2;
typedef GOverlapCache<3> OverlapCache3;

template <unsigned D> struct GProximityCache : public GPairwiseCache<D> {};
typedef GProximityCache<2> ProximityCache2;
typedef GProximityCache<3> ProximityCache3;

}} //namespace geo::np

#endif // GEO_NP_PAIRWISE_CACHE_H
