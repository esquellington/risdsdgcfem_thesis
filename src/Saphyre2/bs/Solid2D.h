#ifndef S2_BS_SOLID2D_H
#define S2_BS_SOLID2D_H

#include <Saphyre2/bs/ISyncEntity.h>
#include <Geo/shape/MeshSolidShape2.h>
#include <Geo/geo.h>

namespace S2
{

//! Solid2D class
/*! App-level Solid2D interface.

  \todo By now Solid is just the same as a Particle system, but should
  have a different API, cannot assume a fixed number of pre-defined
  particles. The API should not assume any internal
  implementation. Forces should be applied at points, not at specific
  particles.

*/
class Solid2D: public ISyncEntity
{
public:
    static const unsigned int cMaxForceApplicationsPerDT = 16;

    //!< Static on-creation params, inited to defaults
    struct Params
    {
        Flags32 m_Flags;
        // Solver params
        Real m_FixedDT;
        Real m_DiscreteLength;
        // Material params
        Real m_Density;
        Real m_Thickness;
        Real m_YoungModulus;
        Real m_PoissonRatio;
        Real m_DampingRatio;
        Real m_PlasticYield;
        Real m_PlasticMax;
        Real m_PlasticCreepPerSecond;

        Params()
        : m_Flags(0)
        , m_FixedDT(0.01f)
        , m_DiscreteLength(0.1f)
        , m_Density(10)
        , m_Thickness(0.01f)
        , m_YoungModulus(10)
        , m_PoissonRatio(0.25)
        , m_DampingRatio(1.0)
        , m_PlasticYield(0)
        , m_PlasticMax(0)
        , m_PlasticCreepPerSecond(0)
            {}
    };

    /*! Extra params... mainly for dev/debug, some may become
       first-class Params
    */
    struct Params_Extra
    {
        Vec2 m_Gravity;
        String32 m_MM;
        String32 m_RM;
        String32 m_DM;
        String32 m_Integrator;
        uint32 m_NR_Iter;
        float m_NR_RelEpsilon;
        String32 m_LS_Solver;
        uint32 m_LS_Iter;
        uint32 m_LS_Restarts;
        float m_LS_RelEpsilon;

        Params_Extra()
        : m_Gravity(0,0)
            {}
    };

    typedef geo::GObjectSS<geo::MeshSolidShape2> geo_object_type;

public:

    EEntityType GetType() const { return eEntity_Solid2D; }
    bool Update( Real dt );

    //! \name Initialization/Edition
    //@{
    bool SetParams( const Params& params ); //!< Returns false if error
    inline const Params& GetParams() const { return m_Params; }
    geo_object_type* CreateMeshGO( geo::ShapeID shape_id ); //\note Shape MUST be a MeshSolidShape2 by now
    geo_object_type* CreateEmbeddedGO( geo::ShapeID shape_id ); //\note Shape MUST be a PolygonalShape2 by now
    //\todo Once the Shape is set, the GO is created and can be retrieved to CHANGE creation transform/SDOF
    //@}

    const geo_object_type* GetMeshGO() const;
    const geo_object_type* GetEmbeddedGO() const;

    //! \name Optimized whole-object edit methods
    //@{
    bool SetPosCoM( const Vec2& pos ); //Moves CoM to specified position
    /*\todo
    void SetVel( const Vec2& vel );
    void ApplyForce( const Vec2& f_global );
    void ApplyImpulse( const Vec2& j_global );
    */
    //@}

    //! \name Continuous-solid methods (mandatory)
    //@{
    /*\todo
    void SetPos( const Point2& solid_pos, const Point2& pos );
    void SetVel( const Point2& solid_pos, const Vec2& vel );
    void SetAcc( const Point2& solid_pos, const Vec2& acc );

    void ApplyForce( const Point2& solid_pos, const Vec2& f_global );
    void ApplyImpulse( const Point2& solid_pos, const Vec2& j_global );
    */
    bool ApplyPressure( const Point2& solid_pos, Real radius, Real pressure );
    //@}

    //! \name Per-node edit methods
    //@{
    //\todo Supporting this assumes that MeshGO nodes will remain
    //unchanged, which disallows fracture, remeshing, etc... we CAN
    //support it, but
    void ApplyForce( uint32 nid, const Vec2& f_global );
    //@}

    //! \name Constraint edit methods
    //@{
    typedef uint32 kpc_id_type;
    kpc_id_type CreateKinematicPointConstraint( const Point2& solid_pos, const Point2& world_pos, const Point2& world_vel );
    void DestroyKinematicPointConstraint( kpc_id_type kpcid );
    void SetKinematicPointConstraint_Pos( kpc_id_type kpcid, const Point2& world_pos ); //KPCID MUST be sync
    void SetKinematicPointConstraint_Vel( kpc_id_type kpcid, const Point2& world_vel ); //KPCID MUST be sync
    struct KinematicPointConstraint
    {
        KinematicPointConstraint( machine_uint_type eid, const Point2& solid_pos )
        : m_EID( eid ), m_SolidPos(solid_pos) {}
        machine_uint_type m_EID;
        Point2 m_SolidPos;
    };
    std::vector< KinematicPointConstraint > m_vecKPC;
    //@}

    //! Extra params... mainly for dev/debug
    bool SetParamsExtra( const Params_Extra& params_extra ); //!< Returns false if error

private:

    //! \name ISyncEntity internal protocol
    //@{
    void BeginDef_Internal();
    bool EndDef_Internal();
    void Lock_Internal() {}
    void Unlock_Internal();
    bool ProcessReturn_Internal( const ds::ReturnIt& rit );
    //@}

    //! \name Constructor/Destructor and Methods with controlled scope
    //@{
    Solid2D( ISyncEntity* p_parent );
    ~Solid2D();
    //@}

    bool ProcessUpdate( const ds::ReturnIt& rit );

    friend class Universe; //!< Allows only Universe to create/destroy instances.

private:
    enum ETouchedAttribute {
        eTouchedNothing     = 0,
        eTouchedPosCoM      = (1<<0),
        eTouchedVel         = (1<<1),
        eTouchedForce       = (1<<2),
        eTouchedImpulse     = (1<<3),
        eTouchedPressure    = (1<<4),
        eTouchedKPC         = (1<<5),
        eTouchedParamsExtra = (1<<6),
        eTouchedAnything    = 0xFFFFFFFF
    };

private:

    //! \name Fixed Parameters
    //@{
    Params m_Params;
    geo::ShapeID m_ShapeId;
    geo_object_type* m_pMeshGO; /*\todo both the ORIGINAL
                                  discretization and the INSTANTANEOUS
                                  deformed configuration should be
                                  accessible so that the user can map
                                  it to a higher resolution
                                  representation. Consider forcing it
                                  to Polygonal2 (for surface-only) or
                                  TriSolid2 (for more expensive
                                  volumetric deformations...) */
    geo_object_type* m_pEmbeddedGO;
    //@}

    //! \name User State
    //@{
    Vec2 m_PosCoM;
    //\note Positions are stored as m_pMeshGO SDOF
    //\todo MAYBE we'll need to store velocities... but not by now...
    //@}

    //! \name Per-particle Touch Events (pos/vel/force/impulse)
    /*!\todo specific TouchedParticle2D class sucks... Use a
      std::triad<uint16 touch_type,uint16 pid, vec2 data> array for
      touch events and, if maximum capacity is hit, END the current
      Edit(), Sync it, clear touch mask, OPEN a new Edit() and write
      the exceeding touch events there!

      \todo Consider a global pool of touch events instead of
      per-entity arrays, though cache behaviour would be a lot
      worse...
    */
    //@{
    unsigned int m_MaxTP;
    unsigned int m_NumTP;
    ds::TouchedParticle2D* m_vecTP;
    //@}

    //! \name Single radial pressure application per frame, by now
    //@{
    ds::RadialPressureAtPoint2D m_RadialPressure;
    //@}
};

} // namespace S2

#endif // S2_SOLID2D_H
