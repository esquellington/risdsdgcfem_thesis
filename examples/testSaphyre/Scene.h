#ifndef TEST_SAPHYRE_SCENE_H
#define TEST_SAPHYRE_SCENE_H

#include <Safra/task/IUpdaterTask.h>

// S2 Includes
#include <Saphyre2/bs/Universe.h>
#include <Saphyre2/bs/Kine2D.h>
#include <Saphyre2/bs/Kine3D.h>
#include <Saphyre2/bs/ParticleSys2D.h>
#include <Saphyre2/bs/Particle3D.h>
#include <Saphyre2/bs/Particle2D.h>
#include <Saphyre2/bs/ParticleSys2D.h>
#include <Saphyre2/bs/Fluid2D.h>
#include <Saphyre2/bs/Solid2D.h>
#include <Saphyre2/bs/Solid3D.h>

#include <map>

class Params;
class SceneRenderer;

/*! Wraps an S2::Universe adding sfr helper functionality:
  - Store S2 objects per category and name
*/
class Scene: public sfr::IUpdateTarget
{
public:
    typedef std::map<std::string, S2::ISyncEntity*> MapSyncEntity;

    //! Simple iterator that wraps MapSyncEntity container internals
    class EntityIterator
    {
    public:
        EntityIterator( const MapSyncEntity &map_sync_entity )
        : m_Container(map_sync_entity)
        , m_Iterator(map_sync_entity.begin())
        {}
        EntityIterator( const EntityIterator &eit ) : m_Container(eit.m_Container), m_Iterator(eit.m_Iterator) {}
        ~EntityIterator() {}
        inline bool IsValid() const { return m_Iterator != m_Container.end(); }
        S2::ISyncEntity *operator*() const { return m_Iterator->second; }
        void operator++() { ++m_Iterator; }

    private:
        const MapSyncEntity &m_Container;
        MapSyncEntity::const_iterator m_Iterator;
    };

public:
    Scene( const Params& params );
    ~Scene();

    void Init();

    bool Update( double dt );

    /* These methods will disappear and their functionality moved to s2sListener

    S2::Rigid *CreateRigid( const S2::Point3 &pos, const S2::Quat &rot,
                            const S2::Vec3 &vel, const S2::Vec3 &vel_rot,
                            S2::Real mass = 1.0f, const S2::Vec3 &diagonal_inertia_tensor = S2::Vec3(1,1,1) );
    S2::Particle *CreateParticleChain( unsigned int num_objects,
                                       float handle_mass,
                                       float mass_factor,
                                       const S2::Point3 &handle_pos,
                                       const S2::Vec3 &direction, float separation );

    S2::Rigid *CreateRigidChain( unsigned int num_objects,
                                 float handle_mass, const S2::Vec3 &handle_diagonal_inertia_tensor,
                                 float mass_factor,
                                 const S2::Point3 &handle_pos, const S2::Quat &handle_orientation,
                                 const S2::Vec3 &direction, float separation );

    S2::Particle *CreateParticleChainRigidUFO( unsigned int num_objects,
                                               float handle_mass,
                                               float mass_factor,
                                               const S2::Point3 &handle_pos,
                                               const S2::Vec3 &direction, float separation,
                                               S2::Rigid *p_rigid = NULL );

    S2::Particle *CreateParticleTree( unsigned int num_levels );
    */

    //! \name s2sListener Needs direct access to S2 API
    //@{
    S2::Universe *GetUniverse() { return m_pUniverse; }

    inline double GetTime() const { return m_Time; }

    bool Clear();
    bool Kill( const char *name );

    S2::ISyncEntity *Find( const char *name ) const
    {
        MapSyncEntity::const_iterator it_s2e;
        it_s2e = m_mapSyncEntity.find(std::string(name));
        return (it_s2e!=m_mapSyncEntity.end()) ? it_s2e->second : 0;
    }
    const char *Add( S2::ISyncEntity *p_sync_entity, const char *name = 0 );
    EntityIterator GetEntityIterator() const { return EntityIterator(m_mapSyncEntity); }
    //@}

private:
    const Params& m_rParams;
    double m_Time;
    S2::Universe *m_pUniverse;
    MapSyncEntity m_mapSyncEntity;
    friend class SceneRenderer;
};

#endif // TEST_SAPHYRE_SCENE_H
