#include "Solid3D.h"
#include <Saphyre2/bs/Universe.h>
#include <Saphyre2/ds/Commands.h>
#include <string.h>

#include <Geo/np/RayCast.h> //TEMP FOR RayCastFace...

namespace S2
{

// Constructor/Destructor and Methods with controlled scope
Solid3D::Solid3D( ISyncEntity* p_parent )
: ISyncEntity(p_parent)
, m_pMeshGO(0)
, m_pEmbeddedGO(0)
, m_PosCoM(0,0,0)
, m_RotCoM(Mat3x3::Identity())
  /*
, m_MaxTP(0)
, m_NumTP(0)
, m_vecTP(0)
  */
{
}

Solid3D::~Solid3D()
{
    if( m_pMeshGO ) S2::BSG::GetObjectFactory().Release( m_pMeshGO );
    if( m_pEmbeddedGO ) S2::BSG::GetObjectFactory().Release( m_pEmbeddedGO );
    //if( m_vecTP ) delete m_vecTP;
}

bool Solid3D::Update( Real dt )
{
    if( !IsValid() ) return false;

    BS_BCMD("Solid3D::Update()");

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
void Solid3D::BeginDef_Internal()
{
    // Default params already set
}

//! \pre IsBeingDefined()
bool Solid3D::EndDef_Internal()
{
    BS_BCMD("Solid3D::EndDef()");

    //---- Send creation command
    GetChannel()->BeginCommand( ds::eCmd_Create );
    {
        GetChannel()->Write( "uid", GetUID() );
        GetChannel()->Write( "peid", GetParent()->GetEID() );
        GetChannel()->Write( "entity_type", (uint32)ds::eEntity_DSH );
        GetChannel()->Write( "dsh_type", (uint32)ds::eDSH_Leaf_Solid3D );

        // Params
        GetChannel()->Write( "flags", m_Params.m_Flags );
        GetChannel()->Write( "fixed_dt", m_Params.m_FixedDT );
        GetChannel()->Write( "density", m_Params.m_Density );
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
void Solid3D::Unlock_Internal()
{
    //---- Send edit command
    GetChannel()->BeginCommand( ds::eCmd_Edit );
    {
        GetChannel()->Write( "eid", GetEID() );
        if( IsTouched(eTouchedPosCoM) ) GetChannel()->Write( "pos_com", m_PosCoM );
        if( IsTouched(eTouchedRotCoM) ) GetChannel()->Write( "rot_com", m_RotCoM );
    }
    GetChannel()->EndCommand();
}

//! \pre ???
bool Solid3D::ProcessReturn_Internal( const ds::ReturnIt& rit )
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



//---- Solid3D Initialization/Edition Methods.
bool Solid3D::SetParams( const Params& params )
{
    if( !IsBeingDefined() )
        return false;

    //\todo CHECK PARAMS CORRECT, IN RANGE, not NaN, etc... BS_ASSERT(  );

    m_Params = params;

    // \todo For large systems, this is an overkill... usually only a
    // few particles will be touched in a user-frame. Use maxtp <<
    // numparticles.
    /*
    m_MaxTP = params.m_NumParticles;
    m_vecTP = new ds::TouchedParticle3D[m_MaxTP];
    */

    return true;
}

Solid3D::geo_object_type* Solid3D::CreateMeshGO( geo::ShapeID shape_id )
{
    if( !IsBeingDefined() ) return 0;
    BS_ASSERT( BSG::GetShapeLibrary().Lookup(shape_id).IsValid() );
    m_ShapeId = shape_id;
    m_pMeshGO = static_cast<geo_object_type*>( S2::BSG::GetObjectFactory().CreateSS( shape_id ) );
    return m_pMeshGO;
}

const Solid3D::geo_object_type* Solid3D::GetMeshGO() const
{
    return m_pMeshGO;
}

Solid3D::geo_object_type *Solid3D::CreateEmbeddedGO( geo::ShapeID shape_id )
{
    if( !IsBeingDefined() ) return 0;
    BS_ASSERT( BSG::GetShapeLibrary().Lookup(shape_id).IsValid() );
    m_pEmbeddedGO = static_cast<geo_object_type*>( S2::BSG::GetObjectFactory().CreateSS( shape_id ) );
    return m_pEmbeddedGO;
}

const Solid3D::geo_object_type *Solid3D::GetEmbeddedGO() const
{
    return m_pEmbeddedGO;
}

//---- Optimized for-each-particle methods
bool Solid3D::SetPosCoM( const Vec3& pos )
{
    if( !IsLocked() && !IsBeingDefined() ) return false;
    Touch( eTouchedPosCoM );
    m_PosCoM = pos;
    return true;
}

bool Solid3D::SetRotCoM( const Mat3x3& rot )
{
    if( !IsLocked() && !IsBeingDefined() ) return false;
    Touch( eTouchedRotCoM );
    m_RotCoM = rot;
    return true;
}

/*
void Solid3D::SetVel( const Vec2& vel )
{
    if( !IsLocked()& & !IsBeingDefined() ) return;
}

void Solid3D::ApplyForce( const Vec2& f_global )
{
    if( !IsLocked() ) return;
}

void Solid3D::ApplyImpulse( const Vec2& j_global )
{
    if( !IsLocked() ) return;
}

//---- Per-Particle methods
const Point2& Solid3D::GetPos( int pid ) const
{
    return m_vecPos[pid];
}

Vec2 Solid3D::GetVel( int pid ) const
{
    return Vec2(0,0);//m_vecVel[pid];
}

void Solid3D::SetPos( int pid, const Point2& pos )
{
    if( !IsLocked() ) return;
    BS_ASSERT( m_NumTP < m_MaxTP ); //!\todo Force Sync() and add TP if full
    m_vecTP[m_NumTP++] = ds::TouchedParticle3D( pid, ds::TouchedParticle3D::eTouchedPos, pos );
    Touch(eTouchedPos);
}
void Solid3D::SetVel( int pid, const Vec2& vel )
{
    if( !IsLocked() ) return;
    BS_ASSERT( m_NumTP < m_MaxTP ); //!\todo Force Sync() and add TP if full
    m_vecTP[m_NumTP++] = ds::TouchedParticle3D( pid, ds::TouchedParticle3D::eTouchedVel, vel );
    Touch(eTouchedVel);
}

void Solid3D::ApplyForce( int pid, const Vec2& f_global )
{
    if( !IsLocked() ) return;
    BS_ASSERT( m_NumTP < m_MaxTP ); //!\todo Force Sync() and add TP if full
    m_vecTP[m_NumTP++] = ds::TouchedParticle3D( pid, ds::TouchedParticle3D::eTouchedForce, f_global );
    Touch(eTouchedForce);
}
void Solid3D::ApplyImpulse( int pid, const Vec2& j_global )
{
    if( !IsLocked() ) return;
    BS_ASSERT( m_NumTP < m_MaxTP ); //!\todo Force Sync() and add TP if full
    m_vecTP[m_NumTP++] = ds::TouchedParticle3D( pid, ds::TouchedParticle3D::eTouchedImpulse, j_global );
    Touch(eTouchedImpulse);
}
*/

Solid3D::kpc_id_type Solid3D::CreateKinematicPointConstraint( const Point3& solid_pos, const Point3& world_pos, const Point3& world_vel )
{
    BS_BCMD( "Creating KPC" );
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

void Solid3D::DestroyKinematicPointConstraint( kpc_id_type kpcid )
{
    BS_BCMD( "Destroying KPC" );
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
void Solid3D::SetKinematicPointConstraint_Pos( kpc_id_type kpcid, const Point3& world_pos )
{
    BS_BCMD( "Editing KPC Pos" );
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
void Solid3D::SetKinematicPointConstraint_Vel( kpc_id_type kpcid, const Point3& world_vel )
{
    BS_BCMD( "Editing KPC Vel" );
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

bool Solid3D::SetParamsExtra( const Params_Extra& params_extra )
{
    BS_BCMD( "Setting Params Extra" );
    BS_ASSERT( IsLocked() );

    // Check params \todo Check ranges, not only NaN
    BS_ASSERT( !mal::IsNaN( params_extra.m_Gravity ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_NR_Iter ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_NR_RelEpsilon ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_LS_Iter ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_LS_Restarts ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_LS_RelEpsilon ) );

    BS_ASSERT( !mal::IsNaN( params_extra.m_CS_DepthMargin ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_CS_PenaltyKs ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_CS_PenaltyKd ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_CS_RelaxCoeff ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_CS_FrictionS ) );
    BS_ASSERT( !mal::IsNaN( params_extra.m_CS_FrictionD ) );

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
            //Extra MIG2015
            GetChannel()->Write( "cs_posm", params_extra.m_CS_PosCST );
            GetChannel()->Write( "cs_velm", params_extra.m_CS_VelCST );
            GetChannel()->Write( "cs_depth_margin", params_extra.m_CS_DepthMargin );
            GetChannel()->Write( "cs_penalty_ks", params_extra.m_CS_PenaltyKs );
            GetChannel()->Write( "cs_penalty_kd", params_extra.m_CS_PenaltyKd );
            GetChannel()->Write( "cs_relax_coeff", params_extra.m_CS_RelaxCoeff );
            GetChannel()->Write( "cs_friction_s", params_extra.m_CS_FrictionS );
            GetChannel()->Write( "cs_friction_d", params_extra.m_CS_FrictionD );
        }
        GetChannel()->EndComplex();
    }
    GetChannel()->EndCommand();
    BS_ECMD(true);
    Touch(eTouchedParamsExtra);
    return true;
}

bool Solid3D::RayCastNode( const Point3& ray_pos, const Vec3& ray_dir, Real thickness, Point3& solid_pos )
{
    Real thickness_sq( mal::Sq(thickness) );
    Real min_dist_sq( thickness_sq );
    Point3 min_dist_node_pos(0,0,0);
    uint32 num_nodes( m_pMeshGO->GetShape()->GetNumV() );
    for( unsigned int it_node=0; it_node<num_nodes; it_node++ )
    {
        Point3 node_pos = m_pMeshGO->GetSDOF(it_node);
        Vec3 diff( node_pos - ray_pos );
        Real dist_sq( mal::NormSq( diff - mal::Dot(diff,ray_dir)*ray_dir ) );
        if( dist_sq < min_dist_sq )
        {
            min_dist_sq = dist_sq;
            min_dist_node_pos = node_pos;
        }
    }
    if( min_dist_sq < thickness_sq )
    {
        solid_pos = min_dist_node_pos;
        return true;
    }
    else
        return false;
}

bool Solid3D::RayCastFace( const Point3& ray_pos, const Vec3& ray_dir, Vec3& hit_pos, Vec3& hit_normal )
{
    const geo::TetSolidShape3* pTSS3( m_pMeshGO->GetShape() );
    uint32 num_bf( pTSS3->GetNumBF() );
    geo::np::RayHit3 rh;
    Real lambda_min(1000);
    for( unsigned int it_bf=0; it_bf<num_bf; it_bf++ )
    {
        Vec3 p0 = m_pMeshGO->GetSDOF( pTSS3->BF_VID(it_bf,0) );
        Vec3 p1 = m_pMeshGO->GetSDOF( pTSS3->BF_VID(it_bf,1) );
        Vec3 p2 = m_pMeshGO->GetSDOF( pTSS3->BF_VID(it_bf,2) );
        if( geo::np::RayCast_Triangle3_SingleSided( ray_pos, ray_dir, Interval(0,1000), //clip at lambda_min?!
                                                    p0, p1, p2,
                                                    rh ) )
        {
            /*
            BS_LOG("RayCastFace %d\n\tlambda = (%f,%f)\n\tp = (%f,%f,%f)\n\tnormal = (%f,%f,%f)",
                   it_bf,
                   rh.m_Interval.Min(), rh.m_Interval.Max(),
                   rh.m_Point[0], rh.m_Point[1], rh.m_Point[2],
                   rh.m_Normal[0], rh.m_Normal[1], rh.m_Normal[2] );
            */
            if( rh.m_Interval.Min() < lambda_min )
            {
                lambda_min = rh.m_Interval.Min();
                hit_normal  = rh.m_Normal;
            }
        }
    }
    hit_pos = ray_pos + lambda_min * ray_dir;
    return lambda_min < 1000;
}

void Solid3D::FixAllNodesInBoundaryPlane( const Point3& ray_pos, const Vec3& ray_dir, Real thickness )
{
    Vec3 hit_pos(0,0,0), hit_normal(0,0,0);
    if( RayCastFace( ray_pos, ray_dir, hit_pos, hit_normal ) )
    {
        /*
        BS_LOG("RayCastFaceHIT!\n\tp = (%f,%f,%f)\n\tnormal = (%f,%f,%f)",
               hit_pos[0], hit_pos[1], hit_pos[2],
               hit_normal[0], hit_normal[1], hit_normal[2] );
        */
        Lock();
        {
            // Create a KNC for ALL nodes in ALL faces coplanar to F
            uint32 num_nodes( m_pMeshGO->GetShape()->GetNumV() );
            for( unsigned int it_node=0; it_node<num_nodes; it_node++ )
            {
                Point3 node_pos = m_pMeshGO->GetSDOF(it_node);
                if( mal::Abs( mal::Dot( node_pos - hit_pos, hit_normal ) ) < thickness  )
                    CreateKinematicPointConstraint( node_pos, node_pos, Vec3::Zero() );
            }
        }
        Unlock(true);
    }
}

//---- Internal ISyncEntity methods
bool Solid3D::ProcessUpdate( const ds::ReturnIt& rit )
{
    const Real* p_state = rit.Find("state").GetArrayPtr<Real>();
    memcpy( m_pMeshGO->GetVecDOF_WriteOnly(),
            p_state,
            m_pMeshGO->GetShape()->GetNumV()*sizeof(Point3) );
    return true;
}

void Solid3D::Dbg_SaveState( util::ItemStream& its ) const
{
    its.WriteArray( "state", m_pMeshGO->GetVecDOF(), m_pMeshGO->GetShape()->GetNumV() );
}

} // namespace S2
