#ifndef S2_BS_RAYCASTQUERY_H
#define S2_BS_RAYCASTQUERY_H

#include <Saphyre2/bs/Config.h>
#include <Saphyre2/bs/BSG.h>
#include <Geo/np/RayCast.h> //TEMP: For RayHit2

namespace S2 {

class ISyncEntity;

/*! IQuery \todo By now, it's just an interface, similar to
  ISyncEntity, but for queries... not sure if it'll get fleshed out,
  but allows ISyncEntity to have a SyncMe(IQuery* p_query) INDEPENDENT
  from SyncMe(ISyncEntity* p_child), which is a GOOD idea
*/
class IQuery
{
public:
    inline IQuery() {}
    virtual ~IQuery() {}
    inline machine_uint_type GetQID() const { return reinterpret_cast<machine_uint_type>(this); }

protected:
    friend class Universe;
    virtual bool ProcessReturn( const ds::ReturnIt& rit ) = 0;
};

//\todo should be in Geo because it may be used internally in specific geo::np::RayCast queries!!
enum ERayCastQueryFlags
{
    eRCQF_ClosestOnly        = 1<<0,
    eRCQF_BoundingVolumeOnly = 1<<1, //!< If defined, only BV hits are returned, fast but inaccurate.
    eRCQF_ComputeNormal      = 1<<2,
    eRCQF_ComputeFeatureId   = 1<<3,
    eRCQF_SortedByLambda     = 1<<4,
    eRCQF_Default            = 0
};

/*! RayCast Query
  Computes ray hits on a Universe
  \todo Maybe could be performed on ANY ISyncEntity, from a whole Universe to single Particle
*/
class RayCastQuery2D: public IQuery
{
public:
    struct Hit2D
    {
        ISyncEntity* m_pEntity;
        //geo::feature_id m_FeatureId; //\todo this should accept ANY feature and internal position on IEntity
        geo::np::RayHit2 m_RayHit;
    };

public:
    void Run( const Vec2& pos, const Vec2& dir,
              const Interval& lambda_interval,
              Real thickness,
              Flags32 flags = eRCQF_Default /*\todo filters, closes only, etc...*/ );
    //! Computes closest hit
    bool Closest( Hit2D& h2d ) const;

    //iterator_ray_hit_2d RayCast2D( ... ); //\todo Better, iterate over all hits, by now

    /*\todo Could "Add" new rays to a single RayCastQuery2D and Sync them in a single call...
      void Begin()
         ray_id Add( const Vec2& pos, const Vec2& dir, Real lambda_max, Real thickness, Flags32 flags );
      void End(bool sync)
      iterator_ray_hit_2d Begin( ray_id ) const
    */

private:
    friend class Universe;
    RayCastQuery2D( ISyncEntity* p_entity, ds::Channel* p_channel ) : m_pEntity(p_entity), m_pChannel(p_channel), m_NumHits(0) {}
    bool ProcessReturn( const ds::ReturnIt& rit );

private:
    ISyncEntity* m_pEntity;
    ds::Channel* m_pChannel;
    uint32 m_NumHits;
    Hit2D m_ClosestH2D;
};

/*! RayCast Query
  Computes ray hits on a Universe
  \todo Maybe could be performed on ANY ISyncEntity, from a whole Universe to single Particle
*/
class RayCastQuery3D: public IQuery
{
public:
    struct Hit3D
    {
        ISyncEntity* m_pEntity;
        //geo::feature_id m_FeatureId; //\todo this should accept ANY feature and internal position on IEntity
        geo::np::RayHit3 m_RayHit;
    };

public:
    void Run( const Vec3& pos, const Vec3& dir,
              const Interval& lambda_interval,
              Real thickness,
              Flags32 flags = eRCQF_Default /*\todo filters, closes only, etc...*/ );
    //! Computes closest hit
    bool Closest( Hit3D& h3d ) const;

    //iterator_ray_hit_2d RayCast3D( ... ); //\todo Better, iterate over all hits, by now

    /*\todo Could "Add" new rays to a single RayCastQuery3D and Sync them in a single call...
      void Begin()
         ray_id Add( const Vec2& pos, const Vec2& dir, Real lambda_max, Real thickness, Flags32 flags );
      void End(bool sync)
      iterator_ray_hit_2d Begin( ray_id ) const
    */

private:
    friend class Universe;
    RayCastQuery3D( ISyncEntity* p_entity, ds::Channel* p_channel ) : m_pEntity(p_entity), m_pChannel(p_channel), m_NumHits(0) {}
    bool ProcessReturn( const ds::ReturnIt& rit );

private:
    ISyncEntity* m_pEntity;
    ds::Channel* m_pChannel;
    uint32 m_NumHits;
    Hit3D m_ClosestH3D;
};

} // namespace S2

#endif // S2_BS_UNIVERSE_H
