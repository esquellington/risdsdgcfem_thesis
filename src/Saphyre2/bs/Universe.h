#ifndef S2_BS_UNIVERSE_H
#define S2_BS_UNIVERSE_H

#include <Saphyre2/bs/Config.h>
#include <Saphyre2/bs/ISyncEntity.h>
#include <util/GPointerContainer.h>

namespace S2 {

// Forward declarations
class Particle3D;
class Particle2D;
class ParticleSys2D;
class Fluid2D;
class Solid2D;
class Solid3D;

class Kine2D;
class Kine3D;

class IQuery;
class RayCastQuery2D;
class RayCastQuery3D;

//! API Universe.
/*! Universe public API for Saphyre2.
  Offers user-level functionality for creating, running and querying
  an animated universe.

  Factory of Animation Objects and Systems.
*/
class Universe: public ISyncEntity
{
public:

    static const unsigned int cDefaultMaxObjects = 64;
    static const unsigned int cDefaultMaxSystems = 16;

    //!< Static on-creation params, inited to defaults
    struct Params
    {
        Flags32 m_Flags;
        unsigned int m_MaxObjects;
        unsigned int m_MaxSystems;

        Real m_LengthScale;
        Real m_MassScale;
        Real m_TimeScale;

        Params()
        : m_Flags(0)
        , m_MaxObjects(cDefaultMaxObjects)
        , m_MaxSystems(cDefaultMaxSystems)
        , m_LengthScale(1.0f)
        , m_MassScale(1.0f)
        , m_TimeScale(1.0f)
        {}
    };

public:
    Universe();
    ~Universe();

    EEntityType GetType() const { return eEntity_Universe; }

    //! \name Universe Initialization/Edition Methods.
    //@{
    bool SetParams( const Params &params ); //!< Returns false if error
    inline const Params &GetParams() const { return m_Params; }

    Particle3D *CreateParticle3D(); //\todo Flags can specify on-creation preferences: ex: IsBullet
    Particle2D *CreateParticle2D();
    ParticleSys2D *CreateParticleSys2D();
    Fluid2D *CreateFluid2D();
    Solid2D *CreateSolid2D();
    Solid3D *CreateSolid3D();

    Kine2D *CreateKine2D(); //flags include eStatic,eBoundary,eSimple/eMultiple
    Kine3D *CreateKine3D(); //flags include eStatic,eBoundary,eSimple/eMultiple
    //@}

    //! \name Simulation methods
    //@{
    bool Update( Real dt ); //unconditional sync, even if universe untouched
    bool ForceSync(); //unconditional sync, even if universe untouched
    //@}

    //! \name Queries
    //@{
    RayCastQuery2D *CreateRayCastQuery2D();
    RayCastQuery3D *CreateRayCastQuery3D();
    //IntersectionQuery *CreateIntersectionQuery();
    //@}

    bool SyncMe( ISyncEntity *p_child );
    bool SyncMe( IQuery *p_query );

private:

    //! \name ISyncEntity internal protocol
    //@{
    void BeginDef_Internal() {}
    bool EndDef_Internal();
    void Lock_Internal() {}
    void Unlock_Internal() {}
    bool Sync_Internal();
    //@}

    bool ProcessUpdate( const ds::ReturnIt &rit );

    friend class BSG;

private:
    Params m_Params;
    Real m_LastTime;
    util::GPointerContainer<ISyncEntity> m_Children; //\todo Should use per-entity-type pools instead/in addition
    util::GPointerContainer<IQuery> m_Queries; //\todo Should use per-entity-type pools instead/in addition
};

} // namespace S2

#endif // S2_BS_UNIVERSE_H
