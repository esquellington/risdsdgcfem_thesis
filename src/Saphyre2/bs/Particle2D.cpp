#include <Saphyre2/bs/Particle2D.h>
#include <Saphyre2/bs/Universe.h>
#include <Saphyre2/ds/Commands.h>

namespace S2
{

// Constructor/Destructor and Methods with controlled scope
Particle2D::Particle2D( ISyncEntity *p_parent )
: ISyncEntity(p_parent)
, m_ShapeId(geo::cInvalidShapeId)
, m_Pos0(Vec2::Zero())
, m_Pos1(Vec2::Zero())
, m_Vel0(Vec2::Zero())
, m_Vel1(Vec2::Zero())
, m_AccForce1(Vec2::Zero())
, m_AccImpulse1(Vec2::Zero())
{
}

Particle2D::~Particle2D()
{
}

bool Particle2D::Update( Real dt )
{
    if( !IsValid() ) return false;
    
    BS_BCMD("Particle2D::Update()");

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
void Particle2D::BeginDef_Internal()
{
    // Default params already set
}

//! \pre IsBeingDefined()
bool Particle2D::EndDef_Internal()
{
    BS_BCMD("Particle2D::EndDef()");

    //---- Send creation command
    GetChannel()->BeginCommand( ds::eCmd_Create );
    {
        GetChannel()->Write( "uid", GetUID() );
        GetChannel()->Write( "peid", GetParent()->GetEID() );
        GetChannel()->Write( "entity_type", (uint32)ds::eEntity_DSH );
        GetChannel()->Write( "dsh_type", (uint32)ds::eDSH_Leaf_Particle2D );
        
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
void Particle2D::Unlock_Internal()
{
    BS_BCMD("Particle2D::EndDef()");
    
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

//---- Particle2D Initialization/Edition Methods
bool Particle2D::SetParams( const Params &params )
{
    if( !IsBeingDefined() ) return false;
    m_Params = params;    
    return true;
}

bool Particle2D::AttachShape( geo::ShapeID shape_id )
{
    if( !IsBeingDefined() ) return false;
    BS_ASSERT( BSG::GetShapeLibrary().Lookup(shape_id).IsValid() );
    m_ShapeId = shape_id;
    return true;
}

bool Particle2D::AttachShape( const geo::IShape &r_shape )
{
    if( !IsBeingDefined() ) return false;
    m_ShapeId = BSG::GetShapeLibrary().Register(&r_shape);
    return true;
}

bool Particle2D::SetPos( const Point2 &pos )
{    
    if( !IsLocked() && !IsBeingDefined() ) return false;
    m_Pos1 = pos;
    m_Pos0 = pos;
    Touch( eTouchedPos );
    return true;    
}

bool Particle2D::SetVel( const Vec2 &vel )
{
    if( !IsLocked() && !IsBeingDefined() ) return false;
    m_Vel1 = vel;
    m_Vel0 = vel;    
    Touch( eTouchedVel );
    return true;
}

bool Particle2D::ApplyForce( const Vec2 &f )
{
    if( !IsLocked() ) return false;    
    m_AccForce1 += f;   
    Touch( eTouchedForce );
    return true;
}

bool Particle2D::ApplyImpulse( const Vec2 &j )
{
    if( !IsLocked() ) return false;    
    m_AccImpulse1 += j;
    Touch( eTouchedForce );
    return true;
}


//---- Internal ISyncEntity methods
bool Particle2D::ProcessUpdate( const ds::ReturnIt &rit )
{
    const Real *p_state = rit.Find("state").GetArrayPtr<Real>();
    m_Pos0.FromArray( &p_state[0] );
    m_Vel0.FromArray( &p_state[2] );
    
    // Reset user forces
    m_AccForce1 = Vec2::Zero();
    m_AccImpulse1 = Vec2::Zero();
    return true;
}
    
} // namespace S2
