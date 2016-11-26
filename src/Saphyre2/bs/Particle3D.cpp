#include <Saphyre2/bs/Particle3D.h>
#include <Saphyre2/bs/Universe.h>
#include <Saphyre2/ds/Commands.h>

namespace S2
{

// Constructor/Destructor and Methods with controlled scope
Particle3D::Particle3D( ISyncEntity *p_parent )
: ISyncEntity(p_parent)
, m_ShapeId(geo::cInvalidShapeId)
, m_Pos0(Vec3::Zero())
, m_Pos1(Vec3::Zero())
, m_Vel0(Vec3::Zero())
, m_Vel1(Vec3::Zero())
, m_AccForce1(Vec3::Zero())
, m_AccImpulse1(Vec3::Zero())
{
}

Particle3D::~Particle3D()
{
}

bool Particle3D::Update( Real dt )
{
    if( !IsValid() ) return false;
    
    BS_BCMD("Particle3D::Update()");

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
void Particle3D::BeginDef_Internal()
{
    // Default params already set
}

//! \pre IsBeingDefined()
bool Particle3D::EndDef_Internal()
{
    BS_BCMD("Particle3D::EndDef()");

    //---- Send creation command
    GetChannel()->BeginCommand( ds::eCmd_Create );
    {
        GetChannel()->Write( "uid", GetUID() );
        GetChannel()->Write( "peid", GetParent()->GetEID() );
        GetChannel()->Write( "entity_type", (uint32)ds::eEntity_DSH );
        GetChannel()->Write( "dsh_type", (uint32)ds::eDSH_Leaf_Particle3D );
        
        // Params
        GetChannel()->Write( "flags", m_Params.m_Flags );
        GetChannel()->Write( "mass", m_Params.m_Mass );
        GetChannel()->Write( "radius", m_Params.m_Radius );
        GetChannel()->Write( "coeff_restitution", m_Params.m_CoeffRestitution );
        
        // Vars
        GetChannel()->Write( "pos", m_Pos1 );
        GetChannel()->Write( "vel", m_Vel1 );

        // Shape Init
        BS_ASSERT( BSG::GetShapeLibrary().Lookup(m_ShapeId).IsValid() );
        if( false ) //\todo accept m_InternalShapeId and pass it here... geo::cInvalidShapeId != m_ShapeId )
            GetChannel()->Write( "shape_id", (uint32)m_ShapeId );
        else
            GetChannel()->WriteItem( "shape_def", BSG::GetShapeLibrary().Lookup(m_ShapeId) );
    }
    GetChannel()->EndCommand();
    
    bool bResult(true);
    BS_ECMD( bResult );
    return bResult;
}

//! \pre IsTouched() and was IsLocked()
void Particle3D::Unlock_Internal()
{
    BS_BCMD("Particle3D::EndDef()");
    
    //---- Send edit command
    GetChannel()->BeginCommand( ds::eCmd_Edit );
    {
        GetChannel()->Write( "eid", GetEID() );
        if( IsTouched(eTouchedPos) ) GetChannel()->Write( "pos", m_Pos1 );
        if( IsTouched(eTouchedVel) ) GetChannel()->Write( "vel", m_Vel1 );
        if( IsTouched(eTouchedForce) ) GetChannel()->Write( "force", m_AccForce1 );
        if( IsTouched(eTouchedImpulse) ) GetChannel()->Write( "impulse", m_AccImpulse1 );
    }
    GetChannel()->EndCommand();

    BS_ECMD(true);    
}

//---- Particle3D Initialization/Edition Methods
bool Particle3D::SetParams( const Params &params )
{
    if( !IsBeingDefined() ) return false;
    m_Params = params;
    return true;
}

bool Particle3D::AttachShape( geo::ShapeID shape_id )
{
    if( !IsBeingDefined() ) return false;
    BS_ASSERT( BSG::GetShapeLibrary().Lookup(shape_id).IsValid() );
    return true;
}

bool Particle3D::AttachShape( const geo::IShape &r_shape )
{
    if( !IsBeingDefined() ) return false;
    m_ShapeId = BSG::GetShapeLibrary().Register(&r_shape);
    return true;
}

bool Particle3D::SetPos( const Point3 &pos )
{    
    if( !IsLocked() && !IsBeingDefined() ) return false;
    m_Pos1 = pos;
    m_Pos0 = pos;
    Touch( eTouchedPos );
    return true;    
}

bool Particle3D::SetVel( const Vec3 &vel )
{
    if( !IsLocked() && !IsBeingDefined() ) return false;
    m_Vel1 = vel;
    m_Vel0 = vel;    
    Touch( eTouchedVel );
    return true;
}

bool Particle3D::ApplyForce( const Vec3 &f )
{
    if( !IsLocked() ) return false;    
    m_AccForce1 += f;   
    Touch( eTouchedForce );
    return true;
}

bool Particle3D::ApplyImpulse( const Vec3 &j )
{
    if( !IsLocked() ) return false;    
    m_AccImpulse1 += j;
    Touch( eTouchedForce );
    return true;
}


//---- Internal ISyncEntity methods
bool Particle3D::ProcessUpdate( const ds::ReturnIt &rit )
{
    const Real *p_state = rit.Find("state").GetArrayPtr<Real>();
    m_Pos0.FromArray( &p_state[0] );
    m_Vel0.FromArray( &p_state[3] );
    
    // Reset user forces
    m_AccForce1 = Vec3::Zero();
    m_AccImpulse1 = Vec3::Zero();
    return true;
}
    
} // namespace S2
