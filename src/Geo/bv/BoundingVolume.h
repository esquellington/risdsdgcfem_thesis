#ifndef GEO_BV_BOUNDING_VOLUME_H
#define GEO_BV_BOUNDING_VOLUME_H

#include <Geo/Config.h>

namespace geo {
namespace bv {

//! BV types (<256)
enum EBVType {

    // Simple BV
    eBV_Point2    = 1,
    eBV_Point3,
    eBV_Sphere2,
    eBV_Sphere3,
    eBV_AABB2,
    eBV_AABB3,
    eBV_DOP2_K4,
    eBV_DOP2_K8,
    eBV_DOP3_K6,
    eBV_DOP3_K14,
    eBV_OBB2,
    eBV_OBB3,
    eBV_LSS2, // Line-Swept-Sphere (aka Capsule)
    eBV_LSS3,

    //\todo Considear adding eBV_HalfSpace2/eBV_HalfSpace3

    // Hierarchical BV \todo These may NOT be basic BV but a different BVH concept
    //eBVH_SphereTree2,
    //eBVH_SphereTree3,
    //...

    // Special BV
    eBV_Void      = 0,        //<! Never overlaps
    eBV_Void2     = eBV_Void,
    eBV_Void3     = eBV_Void,
    eBV_Infinite  = 0xFF,     //<! Overlaps everything except eBV_Void
    eBV_Infinite2 = eBV_Infinite,
    eBV_Infinite3 = eBV_Infinite
};

//! trait for type-reflection. \note Must be defined for all bv types in their .h file
template <typename T> struct bv_type_of { enum { value = eBV_Void }; };

/*! Bounding Volume Base.
  Handles:
  - Builtin timestamp: Allows lazy computations in some places. The BV
                       timestamp will be "touched" whenever it
                       changes. Any observer can compare the timestamp
                       with the last known one (locally cached).
  - Enumerated type: Simplifies reflection and double dispaching while
                     avoiding the need of a virtual interface.
  - Dimension: 2 or 3

  \note No virtual methods, no complex data. This class is just a
  "header" for all BV types. If generic ops are required on BV,
  implement them as free functions as in TestOverlap.

  \note Despite the "I" prefix, IBoundingVolume is not an abstract
  class, but works "mostly" while requiring less memory.
*/
class IBoundingVolume
{
public:
    inline IBoundingVolume( EBVType type, unsigned int dimension )
    : m_TimeStamp(0)
    , m_Type(type)
    , m_Dimension( (2==dimension) ? 0 : 1 )
        { GEO_ASSERT( uint32(type)<256 && (dimension == 2 || dimension == 3) ); }
    inline EBVType GetType() const { return EBVType(m_Type); }
    inline unsigned int GetDimension() const { return (0==m_Dimension) ? 2:3; }
    inline uint32 GetTimeStamp() const { return m_TimeStamp; }
    inline bool IsDirty( uint32 ts ) const { return m_TimeStamp != ts; }
    inline bool TestDirtyAndUpdate( uint32& ts ) const { bool bDirty( m_TimeStamp != ts ); ts = m_TimeStamp; return bDirty; }

    //! \name pseudo-polimorphism
    //@{
    template <typename T> inline bool IsType() const { return GetType() == EBVType( bv_type_of<T>::value ); }
    template <typename T> inline const T& As() const { GEO_ASSERT(IsType<T>()); return *static_cast<const T*>(this); }
    template <typename T> inline T& As() { GEO_ASSERT(IsType<T>()); return *static_cast<T*>(this); }
    //@}

protected:
    inline void Touch() { m_TimeStamp++; }

private:
    unsigned int m_TimeStamp :24; //IMPORTANT: If this is used with some other external TS with different bitcount, they may become out of phase for values in the incompatible range! Consider global cTimestampBits constant!
    unsigned int m_Type      :7;
    unsigned int m_Dimension :1;
};

//\todo HACK
template < unsigned D >
class GBoundingVolumeD: public IBoundingVolume
{
public:
    static const unsigned int cDimension = D;
public:
    inline GBoundingVolumeD( EBVType type ) : IBoundingVolume(type,D) {}
};

typedef GBoundingVolumeD<2> BoundingVolume2;
typedef GBoundingVolumeD<3> BoundingVolume3;

//!\name Special BV types
//@{
class Void2 : public BoundingVolume2 { public: inline Void2() : BoundingVolume2(eBV_Void2) {} };
class Void3 : public BoundingVolume3 { public: inline Void3() : BoundingVolume3(eBV_Void3) {} };
class Infinite2 : public BoundingVolume2 { public: inline Infinite2() : BoundingVolume2(eBV_Infinite2) {} };
class Infinite3 : public BoundingVolume3 { public: inline Infinite3() : BoundingVolume3(eBV_Infinite3) {} };
template <> struct bv_type_of<Void2> { enum { value = eBV_Void2 }; };
template <> struct bv_type_of<Void3> { enum { value = eBV_Void3 }; };
template <> struct bv_type_of<Infinite2> { enum { value = eBV_Infinite2 }; };
template <> struct bv_type_of<Infinite3> { enum { value = eBV_Infinite3 }; };
//@}

//!\name Static polymorphic methods
//@{
// bool TestOverlap( BV1, BV2 );
// bool TestRay( BV, ray_pos, ray_dir, ray_interval );
// bool TestThickRay( BV, ray_pos, ray_dir, ray_interval, ray_thickness );
//@}

}} //namespace geo::bv

#endif // GEO_BV_BOUNDING_VOLUME_H
