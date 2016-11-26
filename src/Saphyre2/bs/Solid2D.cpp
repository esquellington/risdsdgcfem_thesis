#include "Solid2D.h"
#include <Saphyre2/bs/Universe.h>
#include <Saphyre2/ds/Commands.h>
#include <string.h>

namespace S2
{

// Constructor/Destructor and Methods with controlled scope
Solid2D::Solid2D( ISyncEntity *p_parent )
: ISyncEntity(p_parent)
, m_pMeshGO(0)
, m_pEmbeddedGO(0)
, m_PosCoM(0,0)
, m_MaxTP(0)
, m_NumTP(0)
, m_vecTP(0)
{
}

Solid2D::~Solid2D()
{
    if( m_pMeshGO ) S2::BSG::GetObjectFactory().Release( m_pMeshGO );
    if( m_pEmbeddedGO ) S2::BSG::GetObjectFactory().Release( m_pEmbeddedGO );
    if( m_vecTP ) delete[] m_vecTP;
}

bool Solid2D::Update( Real dt )
{
    if( !IsValid() ) return false;

    BS_BCMD("Solid2D::Update");

    GetChannel()->BeginCommand( ds::eCmd_Update );
    {
        GetChannel()->Write( "eid", GetEID() );
        GetChannel()->Write( "dt", (Real)dt );
    }
    GetChannel()->EndCommand();

    // Imperative sync
    bool bResult = Sync_Internal();
    BS_ECMD( bResult );
    return bResult;
}

//---- ISyncEntity internal protocol

//! \pre !IsCreated()
void Solid2D::BeginDef_Internal()
{
    // Default params already set
}

//! \pre IsBeingDefined()
bool Solid2D::EndDef_Internal()
{
    BS_BCMD("Solid2D::EndDef_Internal");

    //---- Send creation command
    GetChannel()->BeginCommand( ds::eCmd_Create );
    {
        GetChannel()->Write( "uid", GetUID() );
        GetChannel()->Write( "peid", GetParent()->GetEID() );
        GetChannel()->Write( "entity_type", (uint32)ds::eEntity_DSH );
        GetChannel()->Write( "dsh_type", (uint32)ds::eDSH_Leaf_Solid2D );

        // Params
        GetChannel()->Write( "flags", m_Params.m_Flags );
        GetChannel()->Write( "fixed_dt", m_Params.m_FixedDT );
        GetChannel()->Write( "density", m_Params.m_Density );
        GetChannel()->Write( "thickness", m_Params.m_Thickness );
        GetChannel()->Write( "young_modulus", m_Params.m_YoungModulus );
        GetChannel()->Write( "poisson_ratio", m_Params.m_PoissonRatio );
        GetChannel()->Write( "damping_ratio", m_Params.m_DampingRatio );
        GetChannel()->Write( "plastic_yield", m_Params.m_PlasticYield );
        GetChannel()->Write( "plastic_max", m_Params.m_PlasticMax );
        GetChannel()->Write( "plastic_creep_per_second", m_Params.m_PlasticCreepPerSecond );

        // Init pose \todo NO, we DO NOT SET CoM at init, must be EXPLICITLY set afterwards
        //GetChannel()->Write( "pos_com", m_PosCoM );

        // Shape Init
        BS_ASSERT( BSG::GetShapeLibrary().Lookup(m_ShapeId).IsValid() );
        if( false ) //\todo accept m_EShapeId and pass it here... geo::cInvalidShapeId != m_ShapeId )
            GetChannel()->Write( "shape_id", (uint32)m_ShapeId );
        else
            GetChannel()->WriteItem( "shape_def", BSG::GetShapeLibrary().Lookup(m_ShapeId) );

        // Vars \todo ALLOW setting initial transform and SDOF from m_pMeshGO?
        //GetChannel()->WriteArray( "pos_i", m_vecPos, m_Params.m_NumParticles );
        //GetChannel()->Write( "transform", transform );
    }
    GetChannel()->EndCommand();

    bool bResult(true);
    BS_ECMD( bResult );
    return bResult;
}

//! \pre IsTouched() and was IsLocked()
void Solid2D::Unlock_Internal()
{
    //---- Send edit command
    GetChannel()->BeginCommand( ds::eCmd_Edit );
    {
        GetChannel()->Write( "eid", GetEID() );
        if( IsTouched(eTouchedPressure) ) GetChannel()->Write( "pressure", m_RadialPressure );
        if( IsTouched(eTouchedPosCoM) ) GetChannel()->Write( "pos_com", m_PosCoM );
        //\todo TouchedNodes...
        if( m_NumTP > 0 )
        {
            GetChannel()->WriteArray( "vec_tp", m_vecTP, m_NumTP );
            m_NumTP = 0;
        }
    }
    GetChannel()->EndCommand();
}

//! \pre ???
bool Solid2D::ProcessReturn_Internal( const ds::ReturnIt &rit )
{
    ds::ReturnIt create_kpc_it = rit.GetSubItem().Find("create_kpc");
    //get internal return type
    if( create_kpc_it.IsValid() )
    {
        kpc_id_type kpc_uid = create_kpc_it.GetSubItem().Find("kpc_uid").Get<kpc_id_type>();
        machine_uint_type kpc_eid = create_kpc_it.GetSubItem().Find("kpc_eid").Get<machine_uint_type>();
        BS_ASSERT( machine_uint_type(-1) == m_vecKPC[kpc_uid].m_EID );
        m_vecKPC[kpc_uid].m_EID = kpc_eid;
        BS_INFO("ISyncEntity::ProcessReturn_Internal() SubEntity[" << GetEID() << "/" << kpc_eid << "] created!");
    }
    //\todo Other internal command responses would be handled here too
    return true;
}



//---- Solid2D Initialization/Edition Methods.
bool Solid2D::SetParams( const Params &params )
{
    if( !IsBeingDefined() )
        return false;

    // Check params \todo Check ranges, not only NaN
    BS_ASSERT( !mal::IsNaN( params.m_FixedDT ) );
    BS_ASSERT( !mal::IsNaN( params.m_DiscreteLength ) );
    BS_ASSERT( !mal::IsNaN( params.m_Density ) );
    BS_ASSERT( !mal::IsNaN( params.m_Thickness ) );
    BS_ASSERT( !mal::IsNaN( params.m_YoungModulus ) );
    BS_ASSERT( !mal::IsNaN( params.m_PoissonRatio ) );
    BS_ASSERT( !mal::IsNaN( params.m_DampingRatio ) );
    BS_ASSERT( !mal::IsNaN( params.m_PlasticYield ) );
    BS_ASSERT( !mal::IsNaN( params.m_PlasticMax ) );
    BS_ASSERT( !mal::IsNaN( params.m_PlasticCreepPerSecond ) );

    m_Params = params;

    // \todo For large systems, this is an overkill... usually only a
    // few particles will be touched in a user-frame. Use maxtp <<
    // numparticles.
    /*
    m_MaxTP = params.m_NumParticles;
    m_vecTP = new ds::TouchedParticle2D[m_MaxTP];
    */

    //\todo Consider return false if params incorrect, ignoring offending values
    return true;
}

Solid2D::geo_object_type *Solid2D::CreateMeshGO( geo::ShapeID shape_id )
{
    if( !IsBeingDefined() ) return 0;
    BS_ASSERT( BSG::GetShapeLibrary().Lookup(shape_id).IsValid() );
    m_ShapeId = shape_id;
    m_pMeshGO = static_cast<geo_object_type*>( S2::BSG::GetObjectFactory().CreateSS( shape_id ) );

    // \todo For large systems, this is an overkill... usually only a
    // few particles will be touched in a user-frame. Use maxtp <<
    // numparticles.
    uint32 num_nodes = m_pMeshGO->GetShape()->GetNumV();
    m_MaxTP = num_nodes;
    m_vecTP = new ds::TouchedParticle2D[num_nodes];

    return m_pMeshGO;
}

const Solid2D::geo_object_type *Solid2D::GetMeshGO() const
{
    return m_pMeshGO;
}

Solid2D::geo_object_type *Solid2D::CreateEmbeddedGO( geo::ShapeID shape_id )
{
    if( !IsBeingDefined() ) return 0;
    BS_ASSERT( BSG::GetShapeLibrary().Lookup(shape_id).IsValid() );
    m_pEmbeddedGO = static_cast<geo_object_type*>( S2::BSG::GetObjectFactory().CreateSS( shape_id ) );
    return m_pEmbeddedGO;
}

const Solid2D::geo_object_type *Solid2D::GetEmbeddedGO() const
{
    return m_pEmbeddedGO;
}

//---- Optimized for-each-particle methods
bool Solid2D::SetPosCoM( const Vec2 &pos )
{
    if( !IsLocked() && !IsBeingDefined() ) return false;
    Touch( eTouchedPosCoM );
    m_PosCoM = pos;
    return true;
}

/*
void Solid2D::SetVel( const Vec2 &vel )
{
    if( !IsLocked() && !IsBeingDefined() ) return;
}

void Solid2D::ApplyForce( const Vec2 &f_global )
{
    if( !IsLocked() ) return;
}

void Solid2D::ApplyImpulse( const Vec2 &j_global )
{
    if( !IsLocked() ) return;
}

//---- Per-Particle methods
const Point2 &Solid2D::GetPos( int pid ) const
{
    return m_vecPos[pid];
}

Vec2 Solid2D::GetVel( int pid ) const
{
    return Vec2(0,0);//m_vecVel[pid];
}

void Solid2D::SetPos( int pid, const Point2 &pos )
{
    if( !IsLocked() ) return;
    BS_ASSERT( m_NumTP < m_MaxTP ); //!\todo Force Sync() and add TP if full
    m_vecTP[m_NumTP++] = ds::TouchedParticle2D( pid, ds::TouchedParticle2D::eTouchedPos, pos );
    Touch(eTouchedPos);
}
void Solid2D::SetVel( int pid, const Vec2 &vel )
{
    if( !IsLocked() ) return;
    BS_ASSERT( m_NumTP < m_MaxTP ); //!\todo Force Sync() and add TP if full
    m_vecTP[m_NumTP++] = ds::TouchedParticle2D( pid, ds::TouchedParticle2D::eTouchedVel, vel );
    Touch(eTouchedVel);
}

void Solid2D::ApplyImpulse( int pid, const Vec2 &j_global )
{
    if( !IsLocked() ) return;
    BS_ASSERT( m_NumTP < m_MaxTP ); //!\todo Force Sync() and add TP if full
    m_vecTP[m_NumTP++] = ds::TouchedParticle2D( pid, ds::TouchedParticle2D::eTouchedImpulse, j_global );
    Touch(eTouchedImpulse);
}
*/

void Solid2D::ApplyForce( uint32 nid, const Vec2& f_global )
{
    BS_BCMD( "Solid2D::ApplyForce" );
    BS_ASSERT( !mal::IsNaN( f_global ) );
    BS_ASSERT( IsLocked() );
    BS_ASSERT( m_NumTP < m_MaxTP ); //!\todo Force Sync() and add TP if full
    m_vecTP[m_NumTP++] = ds::TouchedParticle2D( nid, ds::TouchedParticle2D::eTouchedForce, f_global );
    Touch(eTouchedForce);
    BS_ECMD(true);
}

bool Solid2D::ApplyPressure( const Point2 &solid_pos, Real radius, Real pressure )
{
    BS_BCMD( "Solid2D::ApplyPressure" );
    BS_ASSERT( !mal::IsNaN( solid_pos ) );
    BS_ASSERT( !mal::IsNaN( radius ) );
    BS_ASSERT( !mal::IsNaN( pressure ) );
    BS_ASSERT( IsLocked() );
    Touch( eTouchedPressure );
    m_RadialPressure.m_Pos = solid_pos;
    m_RadialPressure.m_Radius = radius;
    m_RadialPressure.m_Pressure = pressure;
    BS_ECMD(true);
    return true;
}

Solid2D::kpc_id_type Solid2D::CreateKinematicPointConstraint( const Point2 &solid_pos, const Point2 &world_pos, const Point2 &world_vel )
{
    BS_BCMD( "Solid2D::CreateKinematicPointConstraint" );
    BS_ASSERT( !mal::IsNaN( solid_pos ) );
    BS_ASSERT( !mal::IsNaN( world_pos ) );
    BS_ASSERT( !mal::IsNaN( world_vel ) );
    BS_ASSERT( IsLocked() );
    m_vecKPC.push_back( KinematicPointConstraint( machine_uint_type(-1), solid_pos ) );
    kpc_id_type kpc_uid( m_vecKPC.size()-1 );
    GetChannel()->BeginCommand( ds::eCmd_Internal );
    {
        GetChannel()->Write( "uid", GetUID() );
        GetChannel()->Write( "eid", GetEID() );
        GetChannel()->BeginComplex( "create_kpc", eType_Unknown );
        {
            GetChannel()->Write( "kpc_uid", kpc_uid );
            GetChannel()->Write( "solid_pos", solid_pos );
            GetChannel()->Write( "world_pos", world_pos );
            GetChannel()->Write( "world_vel", world_vel );
        }
        GetChannel()->EndComplex();
    }
    GetChannel()->EndCommand();
    BS_ECMD(true);
    Touch(eTouchedKPC);
    return kpc_uid;
}

void Solid2D::DestroyKinematicPointConstraint( kpc_id_type kpcid )
{
    BS_BCMD( "Solid2D::DestroyKinematicPointConstraint" );
    BS_INFO( "KPC = " << kpcid );
    BS_ASSERT( IsLocked() );
    if( machine_uint_type(-1) == m_vecKPC[kpcid].m_EID )
    {
        BS_LOG_ERROR( "KPC %d has no valid EID, it's either unborn or destroyed", kpcid );
        BS_ECMD(false);
        return;
    }
    BS_ASSERT( kpcid < m_vecKPC.size() && machine_uint_type(-1) != m_vecKPC[kpcid].m_EID );
    GetChannel()->BeginCommand( ds::eCmd_Internal );
    {
        GetChannel()->Write( "uid", GetUID() );
        GetChannel()->Write( "eid", GetEID() );
        GetChannel()->BeginComplex( "destroy_kpc", eType_Unknown );
        {
            GetChannel()->Write( "kpc_eid", m_vecKPC[kpcid].m_EID );
        }
        GetChannel()->EndComplex();
    }
    GetChannel()->EndCommand();
    // Erase local KPC entry NOOO this changes all other indices!!
    //m_vecKPC.erase( m_vecKPC.begin() + kpcid );
    //TEMP: mark it deleted... should remove it but then other KPC indices would change!!!...
    m_vecKPC[kpcid].m_EID = machine_uint_type(-1);
    BS_ECMD(true);
    Touch(eTouchedKPC);
}

//KPCID MUST be sync
void Solid2D::SetKinematicPointConstraint_Pos( kpc_id_type kpcid, const Point2 &world_pos )
{
    BS_BCMD( "Solid2D::SetKinematicPointConstraint_Pos" );
    BS_INFO( "KPC = " << kpcid );
    BS_ASSERT( !mal::IsNaN( world_pos ) );
    BS_ASSERT( IsLocked() );
    if( machine_uint_type(-1) == m_vecKPC[kpcid].m_EID )
    {
        BS_LOG_ERROR( "KPC %d has no valid EID, it's either unborn or destroyed", kpcid );
        BS_ECMD(false);
        return;
    }
    BS_ASSERT( kpcid < m_vecKPC.size() && machine_uint_type(-1) != m_vecKPC[kpcid].m_EID );
    GetChannel()->BeginCommand( ds::eCmd_Internal );
    {
        GetChannel()->Write( "uid", GetUID() );
        GetChannel()->Write( "eid", GetEID() );
        GetChannel()->BeginComplex( "edit_kpc", eType_Unknown );
        {
            GetChannel()->Write( "kpc_eid", m_vecKPC[kpcid].m_EID );
            GetChannel()->Write( "world_pos", world_pos );
        }
        GetChannel()->EndComplex();
    }
    GetChannel()->EndCommand();
    BS_ECMD(true);
    Touch(eTouchedKPC);
}
//KPCID MUST be sync
void Solid2D::SetKinematicPointConstraint_Vel( kpc_id_type kpcid, const Point2 &world_vel )
{
    BS_BCMD( "Solid2D::SetKinematicPointConstraint_Vel" );
    BS_INFO( "KPC = " << kpcid );
    BS_ASSERT( !mal::IsNaN( world_vel ) );
    BS_ASSERT( IsLocked() );
    if( machine_uint_type(-1) == m_vecKPC[kpcid].m_EID )
    {
        BS_LOG_ERROR( "KPC %d has no valid EID, it's either unborn or destroyed", kpcid );
        BS_ECMD(false);
        return;
    }
    BS_ASSERT( kpcid < m_vecKPC.size() && machine_uint_type(-1) != m_vecKPC[kpcid].m_EID );
    GetChannel()->BeginCommand( ds::eCmd_Internal );
    {
        GetChannel()->Write( "uid", GetUID() );
        GetChannel()->Write( "eid", GetEID() );
        GetChannel()->BeginComplex( "edit_kpc", eType_Unknown );
        {
            GetChannel()->Write( "kpc_eid", m_vecKPC[kpcid].m_EID );
            GetChannel()->Write( "world_vel", world_vel );
        }
        GetChannel()->EndComplex();
    }
    GetChannel()->EndCommand();
    BS_ECMD(true);
    Touch(eTouchedKPC);
}


bool Solid2D::SetParamsExtra( const Params_Extra &params_extra )
{
    BS_BCMD( "Solid2D::SetParamsExtra" );
    BS_ASSERT( IsLocked() );

    // Check params \todo Check ranges, not only NaN
    BS_ASSERT( !mal::IsNaN( params_extra.m_Gravity ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_NR_Iter ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_NR_RelEpsilon ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_LS_Iter ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_LS_Restarts ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_LS_RelEpsilon ) );

    GetChannel()->BeginCommand( ds::eCmd_Internal );
    {
        GetChannel()->Write( "uid", GetUID() );
        GetChannel()->Write( "eid", GetEID() );
        GetChannel()->BeginComplex( "set_params_extra", eType_Unknown );
        {
            GetChannel()->Write( "gravity", params_extra.m_Gravity );
            GetChannel()->Write( "mm", params_extra.m_MM );
            GetChannel()->Write( "rm", params_extra.m_RM );
            GetChannel()->Write( "dm", params_extra.m_DM );
            GetChannel()->Write( "integrator", params_extra.m_Integrator );
            GetChannel()->Write( "nr_iter", params_extra.m_NR_Iter );
            GetChannel()->Write( "nr_rel_epsilon", params_extra.m_NR_RelEpsilon );
            GetChannel()->Write( "ls_solver", params_extra.m_LS_Solver );
            GetChannel()->Write( "ls_iter", params_extra.m_LS_Iter );
            GetChannel()->Write( "ls_restarts", params_extra.m_LS_Restarts );
            GetChannel()->Write( "ls_rel_epsilon", params_extra.m_LS_RelEpsilon );
        }
        GetChannel()->EndComplex();
    }
    GetChannel()->EndCommand();
    BS_ECMD(true);
    Touch(eTouchedParamsExtra);
    return true;
}

//---- Internal ISyncEntity methods
bool Solid2D::ProcessUpdate( const ds::ReturnIt &rit )
{
    const Real *p_state = rit.Find("state").GetArrayPtr<Real>();
    memcpy( m_pMeshGO->GetVecDOF_WriteOnly(),
            p_state,
            m_pMeshGO->GetShape()->GetNumV()*sizeof(Point2) );
    return true;
}

} // namespace S2
