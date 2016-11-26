#ifndef S2_BS_SOLID3D_H
#define S2_BS_SOLID3D_H

#include <Saphyre2/bs/ISyncEntity.h>
#include <Geo/shape/TetSolidShape3.h>
#include <Geo/geo.h>

namespace S2
{

//! Solid3D class
/*! App-level Solid3D interface.

  \todo By now Solid is just the same as a Particle system, but should
  have a different API, cannot assume a fixed number of pre-defined
  particles. The API should not assume any internal
  implementation. Forces should be applied at points, not at specific
  particles.

*/
class Solid3D: public ISyncEntity
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
        , m_Density(1)
        , m_YoungModulus(1000)
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
        Vec3f m_Gravity;
        String32 m_MM;
        String32 m_RM;
        String32 m_DM;
        String32 m_Integrator;
        uint32 m_NR_Iter;
        float32 m_NR_RelEpsilon;
        String32 m_LS_Solver;
        uint32 m_LS_Iter;
        uint32 m_LS_Restarts;
        float32 m_LS_RelEpsilon;
        //Extra MIG2015
        String32 m_CS_PosCST;
        String32 m_CS_VelCST;
        float32 m_CS_DepthMargin;
        float32 m_CS_PenaltyKs;
        float32 m_CS_PenaltyKd;
        float32 m_CS_RelaxCoeff;
        float32 m_CS_FrictionS;
        float32 m_CS_FrictionD;

        Params_Extra()
        : m_Gravity(0,0,0)
            {}
    };

    typedef geo::GObjectSS<geo::TetSolidShape3> geo_object_type;

public:

    EEntityType GetType() const { return eEntity_Solid3D; }
    bool Update( Real dt );

    //! \name Initialization/Edition
    //@{
    bool SetParams( const Params& params ); //!< Returns false if error
    inline const Params& GetParams() const { return m_Params; }
    geo_object_type* CreateMeshGO( geo::ShapeID shape_id ); //\note Shape MUST be a TetSolidShape3 by now
    geo_object_type* CreateEmbeddedGO( geo::ShapeID shape_id ); //\note Shape MUST be a TriSurfaceShape3 by now
    //\todo Once the Shape is set, the GO is created and can be retrieved to CHANGE creation transform/SDOF
    //@}
    const geo_object_type* GetMeshGO() const;
    const geo_object_type* GetEmbeddedGO() const;

    //! \name Optimized whole-object edit methods
    //@{
    bool SetPosCoM( const Vec3& pos ); //Moves CoM to specified position
    bool SetRotCoM( const Mat3x3& rot ); //Rotates around CoM
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
    //@}

    //! \name Per-node edit methods
    //@{
    //\todo NO
    //@}

    //! \name Constraint edit methods
    //@{
    typedef uint32 kpc_id_type;
    kpc_id_type CreateKinematicPointConstraint( const Point3& solid_pos, const Point3& world_pos, const Point3& world_vel );
    void DestroyKinematicPointConstraint( kpc_id_type kpcid );
    void SetKinematicPointConstraint_Pos( kpc_id_type kpcid, const Point3& world_pos ); //KPCID MUST be sync
    void SetKinematicPointConstraint_Vel( kpc_id_type kpcid, const Point3& world_vel ); //KPCID MUST be sync
    struct KinematicPointConstraint
    {
        KinematicPointConstraint( machine_uint_type eid, const Point3& solid_pos )
        : m_EID( eid ), m_SolidPos(solid_pos) {}
        machine_uint_type m_EID;
        Point3 m_SolidPos;
    };
    std::vector< KinematicPointConstraint > m_vecKPC;

    //! Extra params... mainly for dev/debug
    bool SetParamsExtra( const Params_Extra& params_extra ); //!< Returns false if error

    //TEMP: FOR DCLFEM
    bool RayCastNode( const Point3& ray_pos, const Vec3& ray_dir, Real thickness, Point3& solid_pos ); //TEMP to pick a node
    bool RayCastFace( const Point3& ray_pos, const Vec3& ray_dir, Vec3& hit_pos, Vec3& hit_normal ); //TEMP to pick a face
    void FixAllNodesInBoundaryPlane( const Point3& ray_pos, const Vec3& ray_dir, Real thickness ); //TEMP to create KNC FOR ALL nodes in a plane
    //@}

    void Dbg_SaveState( util::ItemStream& its ) const;

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
    Solid3D( ISyncEntity* p_parent );
    ~Solid3D();
    //@}

    bool ProcessUpdate( const ds::ReturnIt& rit );

    friend class Universe; //!< Allows only Universe to create/destroy instances.

private:
    enum ETouchedAttribute {
        eTouchedNothing     = 0,
        eTouchedPosCoM      = (1<<0),
        eTouchedRotCoM      = (1<<1),
        eTouchedVel         = (1<<2),
        eTouchedForce       = (1<<3),
        eTouchedImpulse     = (1<<4),
        eTouchedKPC         = (1<<5),
        eTouchedParamsExtra = (1<<6),
        eTouchedAnything = 0xFFFFFFFF
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
                                  representation. */
    geo_object_type* m_pEmbeddedGO;
    //@}

    //! \name User State
    //@{
    //\note Positions are stored as m_pMeshGO SDOF
    //\todo MAYBE we'll need to store velocities... but not by now...
    //@}

    //! \name Per-particle Touch Events (pos/vel/force/impulse)
    /*!\todo specific TouchedParticle3D class sucks... Use a
      std::triad<uint16 touch_type,uint16 pid, vec2 data> array for
      touch events and, if maximum capacity is hit, END the current
      Edit(), Sync it, clear touch mask, OPEN a new Edit() and write
      the exceeding touch events there!

      \todo Consider a global pool of touch events instead of
      per-entity arrays, though cache behaviour would be a lot
      worse...
    */
    //@{
    /*\todo
    unsigned int m_MaxTP;
    unsigned int m_NumTP;
    */
    //@}

    Vec3 m_PosCoM;
    Mat3x3 m_RotCoM;
};

} // namespace S2

#endif // S2_SOLID3D_H
