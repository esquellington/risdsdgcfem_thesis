#ifndef S2_DS_COMMANDS_H
#define S2_DS_COMMANDS_H

#include <Saphyre2/ds/Config.h>

namespace S2 { namespace ds {

typedef int32 CommandID;

//! DS Command Types
enum ECommandType {

    /* Control commands \todo
    eCtrl_Log,
    eCtrl_Viz,
    eCtrl_Prof,
    */

    // Entity commands
    eCmd_Create,         //!< Par = {Id,ParentEID,Type,Def}; Ret = NewEntity {Id,EID,Type,[Def]} / Error
    eCmd_Destroy,        //!< Par = {Id,EID}        Ret = None | Error
    eCmd_Edit,           //!< Par = {Id,EID}        Ret = None | Error
    eCmd_Update,         //!< Par = {Id,EID,dt}     Ret = Flat list of {Update {Id,Data}}
    eCmd_Internal,       //!< Par = {EID,...}       Ret = {...}

    // Queries
    eQuery_RayCast,
    eQuery_Intersection,

    // Debug commands
    eDbg_Conf,           //!< Par = {Id,EID,op={on/off},flags}; Ret = Log/Viz/Prof Items

    eDbg_QueryStats,     //!< Par = {Id,EID,flags};             Ret = Stats Items

    eDbg_QueryParams,    //!< Par = {Id,EID,flags};             Ret = Params Items
    eDbg_SyncParams,     //!< Par = {Id,EID};                   Ret = Params Items

    eNumCommandTypes
};

//! DS Return Types
enum EReturnType {

    // Simpe answers (May be complex and contain data such as error-id or requested values)
    eRet_Ok = 0,
    eRet_Error,

    // Complex answers
    eRet_NewEntity,
    eRet_KillEntity,
    eRet_Update,
    eRet_Internal,

    // Query answers
    eRet_RayCast,
    eRet_Intersection,

    // Debug answers
    eRet_Stats,
    eRet_Params,

    eNumReturnTypes
};

enum EEntityType {

    eEntity_DSH = 0,
    eEntity_Connector,
    eEntity_Constraint,
    eEntity_Interaction,
    eEntity_Geom,

    eNumEntityTypes
};

enum EDSHType {
    eDSH_Generic = 0,

    // Top-level
    eDSH_Multiverse,
    eDSH_Universe2D,
    eDSH_Universe3D,

    // Containers
    eDSH_Aggregate2D,
    eDSH_Aggregate3D,

    // Leafs
    eDSH_Leaf_Particle2D,
    eDSH_Leaf_ParticleSystem2D,
    eDSH_Leaf_Rigid2D,
    eDSH_Leaf_RigidSystem2D,
    eDSH_Leaf_Curve2D,
    eDSH_Leaf_Solid2D,
    eDSH_Leaf_Fluid2D,

    eDSH_Leaf_Particle3D,
    eDSH_Leaf_ParticleSystem3D,
    eDSH_Leaf_Rigid3D,
    eDSH_Leaf_RigidSystem3D,
    eDSH_Leaf_Curve3D,
    eDSH_Leaf_Surface3D,
    eDSH_Leaf_Solid3D,
    eDSH_Leaf_Fluid3D,

    eNumDSHTypes
};


enum EGeomType {
    eGeom_Simple = 0,
    eGeom_Multiple,

    eNumGeomTypes
};


//! Debug Modules
enum EDbgModule {
    eDbgModule_Log,
    eDbgModule_Viz,
    eDbgModule_Prof,
    eDbgModule_Stats,
    eDbgModule_Params
};

/*! Debug Flags per Entity. Used independently for both Log and Viz.
  Total of 16 flags at most:
  - Entity-generic flags use bits [15..0]
  - Entity-specific flags use bits [31..16]
*/
enum EDbgFlags {
    //---- Entity-generic
    eDbg_Nothing       = 0,
    eDbg_State         = (1<<16),     //!< Pos/Vel/etc...
    eDbg_Last          = eDbg_State,  //!< Signals last generic flag

    //---- DSH-Specific (start at 2)
    eDbg_DSH_Shape         = eDbg_Last << 1, //!< Bounding Volume (Viz, mostly)
    eDbg_DSH_Energy        = eDbg_Last << 2, //!< T,U i E=T+U separatedly
    eDbg_DSH_SimScheme     = eDbg_Last << 3,
    eDbg_DSH_BV            = eDbg_Last << 4, //!< Bounding Volume (Viz, mostly)

    /* \todo INADEQUATE flags. Always Recursive by now
    eDbg_DSH_Children      = eDbg_Last << 5,
    eDbg_DSH_Connectors    = eDbg_Last << 6,
    eDbg_DSH_Constraints   = eDbg_Last << 7,
    eDbg_DSH_Interactions  = eDbg_Last << 8,
    eDbg_DSH_Geoms         = eDbg_Last << 9,
    eDbg_DSH_Recursive     = (eDbg_DSH_Children
                              | eDbg_DSH_Connectors
                              | eDbg_DSH_Constraints
                              | eDbg_DSH_Interactions
                              | eDbg_DSH_Geoms),
    */

    //---- Connector-Specific
    eDbg_Connector_XXX1    = eDbg_Last << 1,
    eDbg_Connector_XXX2    = eDbg_Last << 2,

    //---- Constraint-Specific
    eDbg_Constraint_Error  = eDbg_Last << 1,
    eDbg_Constraint_DOF    = eDbg_Last << 2,
    eDbg_Constraint_Limits = eDbg_Last << 3,

    //---- Interaction-Specific

    //---- Geom-Specific
    eDbg_Geom_BV           = eDbg_Last << 1, //!< Bounding Volume (Viz, mostly)

    eDbg_All           = 0xFFFF0000
};

//! Profile Flags per Entity
enum EProfFlags {
    //---- Entity-generic
    eProf_Nothing       = 0,
    eProf_Create        = (1<<16),
    eProf_Edit          = (1<<17),
    eProf_Update        = (1<<18),

    eProf_All           = 0xFFFF0000
};

//---- Specific DSH subclass flags
enum EDbgFlags_Model {

    //---- Model-generic
    eDbg_Model_Nothing        = 0,
    eDbg_Model_Forces         = 1,
    eDbg_Model_Impulses       = (1<<1),
    eDbg_Model_Last           = eDbg_Model_Impulses,

    // Particle-specific
    eDbg_Particle_XXX1        = eDbg_Model_Last << 1,

    // SPH-specific
    eDbg_Fluid_SPH_SCG        = eDbg_Model_Last << 1,
    eDbg_Fluid_SPH_Neighbours = eDbg_Model_Last << 2,
    eDbg_Fluid_SPH_Field      = eDbg_Model_Last << 3,        //!< Kernels and Density
    eDbg_Fluid_SPH_Collisions = eDbg_Model_Last << 4,

    eDbg_Model_All            = 0x0000FFFF
};


//---- Common data types used in Channels (both Commands and Returns)
//TEMPORAL: This is really ugly, consider sending Complex items instead of custom structs...
struct ForceAtPoint3D { Point3 m_Pos; Vec3 m_Force; };
struct ImpulseAtPoint3D { Point3 m_Pos; Vec3 m_Impulse; };
struct RadialPressureAtPoint2D { Point2 m_Pos; Real m_Radius; Real m_Pressure; };

template< unsigned int N >
struct GTouchedParticle {
    enum EType { eTouchedPos=0,
                 eTouchedVel=(1<<30),
                 eTouchedForce=(2<<30),
                 eTouchedImpulse=(3<<30) }; //! We use 2 upper bits
    static const unsigned int cMaskUpper2b = 0xC0000000;
    static const unsigned int cMaskLower30b = ~0xC0000000;
    finline GTouchedParticle() {}
    finline GTouchedParticle( unsigned int pid, EType type, const mal::GVec<Real,N> &vec )
    : m_Bits( (cMaskLower30b & pid) | (cMaskUpper2b & uint32(type)) ), m_Vec(vec) {}
    finline uint32 GetParticleId() const { return m_Bits & cMaskLower30b; }
    finline EType GetType() const { return EType(m_Bits & cMaskUpper2b); }
    uint32 m_Bits; mal::GVec<Real,N> m_Vec;
};
typedef GTouchedParticle<2> TouchedParticle2D;
typedef GTouchedParticle<3> TouchedParticle3D;

//! Enumeration of additional channel data types
enum EChannelDataPlaTypes {
    eType_ForceAtPoint3D = cNumPlaTypes + 1,
    eType_ImpulseAtPoint3D,
    eType_RadialPressureAtPoint2D,
    eType_TouchParticle2D,

    eNumChannelDataTypes //!< Number of standard types + basic types + channel data types
};


} } // namespace S2::ds

//!\name Traits of additional channel data types (MUST BE IN GLOBAL NAMESPACE)
//@{
template <> struct pla_type_id< S2::ds::ForceAtPoint3D > { static const uint16 value = S2::ds::eType_ForceAtPoint3D; };
template <> struct pla_type_id< S2::ds::ImpulseAtPoint3D > { static const uint16 value = S2::ds::eType_ImpulseAtPoint3D; };
template <> struct pla_type_id< S2::ds::RadialPressureAtPoint2D > { static const uint16 value = S2::ds::eType_RadialPressureAtPoint2D; };
template <> struct pla_type_id< S2::ds::TouchedParticle2D > { static const uint16 value = S2::ds::eType_TouchParticle2D; };
//@}

#endif // S2_DS_COMMANDS_H
