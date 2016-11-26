#include <Saphyre2/bs/Fluid2D.h>
#include <Saphyre2/bs/Universe.h>
#include <Saphyre2/ds/Commands.h>
#include <string.h>

namespace S2
{

//---- Constructor
/*! Creates BS object, but not DS object
  \sa Object initialization.
*/
Fluid2D::Fluid2D( ISyncEntity *p_parent )
: ISyncEntity(p_parent)
, m_vecPos(0)
{
}

Fluid2D::~Fluid2D()
{
    if( m_vecPos ) delete m_vecPos;
}

bool Fluid2D::Update( Real dt )
{
    BS_BCMD("Fluid2D::Update()");

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
void Fluid2D::BeginDef_Internal()
{
    // Default params already set
}

//! \pre IsBeingDefined()
bool Fluid2D::EndDef_Internal()
{
    BS_BCMD("Fluid2D::EndDef()");

    //---- Send creation command
    GetChannel()->BeginCommand( ds::eCmd_Create );
    {
        GetChannel()->Write( "uid", GetUID() );
        GetChannel()->Write( "peid", GetParent()->GetEID() );
        GetChannel()->Write( "entity_type", (uint32)ds::eEntity_DSH );
        GetChannel()->Write( "dsh_type", (uint32)ds::eDSH_Leaf_Fluid2D );
        
        // Params        
        GetChannel()->Write( "flags", m_Params.m_Flags );
        GetChannel()->Write( "num_particles", (uint32)m_Params.m_NumParticles );
        GetChannel()->Write( "density", m_Params.m_Density );
        GetChannel()->Write( "thickness", m_Params.m_Thickness );
        GetChannel()->Write( "init_aabb_posmin", m_Params.m_InitShape_AABB_PosMin );
        GetChannel()->Write( "init_aabb_posmax", m_Params.m_InitShape_AABB_PosMax );
        GetChannel()->Write( "bounds_aabb_posmin", m_Params.m_Bounds_AABB_PosMin );
        GetChannel()->Write( "bounds_aabb_posmax", m_Params.m_Bounds_AABB_PosMax );
    }
    GetChannel()->EndCommand();
    
    bool bResult(true);
    BS_ECMD( bResult );
    return bResult;
}

//! \pre IsTouched() and was IsLocked()
void Fluid2D::Unlock_Internal()
{
    BS_BCMD("Fluid2D::EndDef()");

    //---- Send edit command
    GetChannel()->BeginCommand( ds::eCmd_Edit );
    {
        GetChannel()->Write( "eid", GetEID() );

        if( IsTouched(eTouchedPressure) ) GetChannel()->Write( "pressure", m_RadialPressure );

        //! \todo Reset accumulators...
    }
    GetChannel()->EndCommand();
}

//---- Fluid2D Initialization/Edition Methods.
bool Fluid2D::SetParams( const Params &params )
{
    if( !IsBeingDefined() ) return false;
    m_Params = params;
    m_vecPos = new Point2[m_Params.m_NumParticles];
    return true;
}

//---- \name Query
const Point2 &Fluid2D::GetPos( int idx ) const
{
    return m_vecPos[idx];
}

bool Fluid2D::ApplyPressure( const Point2 &fluid_pos, Real radius, Real pressure )
{
    BS_ASSERT( IsLocked() );
    Touch( eTouchedPressure );
    m_RadialPressure.m_Pos = fluid_pos;
    m_RadialPressure.m_Radius = radius;
    m_RadialPressure.m_Pressure = pressure;
    return true;
}

//---- Internal ISyncEntity methods
bool Fluid2D::ProcessUpdate( const ds::ReturnIt &rit )
{
    const Real *p_state = rit.Find("state").GetArrayPtr<Real>();
    // Update
    memcpy( &m_vecPos[0], p_state, m_Params.m_NumParticles*sizeof(Point2) );
    // \todo Reset user forces
    return true;
}

} // namespace S2
