#include "ParticleSys2D.h"
#include <Saphyre2/bs/Universe.h>
#include <Saphyre2/ds/Commands.h>
#include <string.h>

namespace S2
{

// Constructor/Destructor and Methods with controlled scope
ParticleSys2D::ParticleSys2D( ISyncEntity *p_parent )
: ISyncEntity(p_parent)
, m_vecPos(0)
, m_MaxTP(0)
, m_NumTP(0)
, m_vecTP(0)
{
}

ParticleSys2D::~ParticleSys2D()
{
    if( m_vecPos ) delete[] m_vecPos;
    if( m_vecTP ) delete[] m_vecTP;
}

bool ParticleSys2D::Update( Real dt )
{
    if( !IsValid() ) return false;

    BS_BCMD("ParticleSys2D::Update()");

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
void ParticleSys2D::BeginDef_Internal()
{
    // Default params already set
}

//! \pre IsBeingDefined()
bool ParticleSys2D::EndDef_Internal()
{
    BS_BCMD("ParticleSys2D::EndDef()");

    //---- Send creation command
    GetChannel()->BeginCommand( ds::eCmd_Create );
    {
        GetChannel()->Write( "uid", GetUID() );
        GetChannel()->Write( "peid", GetParent()->GetEID() );
        GetChannel()->Write( "entity_type", (uint32)ds::eEntity_DSH );
        GetChannel()->Write( "dsh_type", (uint32)ds::eDSH_Leaf_ParticleSystem2D );

        // Params
        GetChannel()->Write( "flags", m_Params.m_Flags );
        GetChannel()->Write( "count", m_Params.m_NumParticles );
        GetChannel()->Write( "mass", m_Params.m_TotalMass );
        GetChannel()->Write( "radius", m_Params.m_ParticleRadius );
        GetChannel()->Write( "coeff_restitution", m_Params.m_ParticleCoeffRestitution );

        // Vars
        GetChannel()->WriteArray( "pos_i", m_vecPos, m_Params.m_NumParticles );
    }
    GetChannel()->EndCommand();

    bool bResult(true);
    BS_ECMD( bResult );
    return bResult;
}

//! \pre IsTouched() and was IsLocked()
void ParticleSys2D::Unlock_Internal()
{
    // Send changes command if any
    if( m_NumTP > 0 )
    {
        GetChannel()->BeginCommand( ds::eCmd_Edit );
        {
            GetChannel()->Write( "eid", GetEID() );
            GetChannel()->WriteArray("vec_tp",m_vecTP,m_NumTP);
        }
        GetChannel()->EndCommand();

        // clear user touch events
        m_NumTP = 0;
    }
}

//---- ParticleSys2D Initialization/Edition Methods.
bool ParticleSys2D::SetParams( const Params &params )
{
    if( !IsBeingDefined() )
        return false;

    BS_ASSERT( params.m_NumParticles < (1<<16) );

    m_Params = params;
    m_vecPos = new Point2[params.m_NumParticles];

    // \todo For large systems, this is an overkill... usually only a
    // few particles will be touched in a user-frame. Use maxtp <<
    // numparticles.
    m_MaxTP = params.m_NumParticles;
    m_vecTP = new ds::TouchedParticle2D[m_MaxTP];

    return true;
}

void ParticleSys2D::DefineParticle( unsigned int pid, const Point2 &pos )
{
    if( !IsBeingDefined() )
        return;
    BS_ASSERT( pid < m_Params.m_NumParticles );
    m_vecPos[pid] = pos;
}

//---- Optimized for-each-particle methods
void ParticleSys2D::SetVel( const Vec2 &vel )
{
    if( !IsLocked() && !IsBeingDefined() ) return;
}

void ParticleSys2D::ApplyForce( const Vec2 &f_global )
{
    if( !IsLocked() ) return;
}

void ParticleSys2D::ApplyImpulse( const Vec2 &j_global )
{
    if( !IsLocked() ) return;
}

//---- Per-Particle methods
const Point2 &ParticleSys2D::GetPos( int pid ) const
{
    return m_vecPos[pid];
}

Vec2 ParticleSys2D::GetVel( int pid ) const
{
    return Vec2(0,0);//m_vecVel[pid];
}

void ParticleSys2D::SetPos( int pid, const Point2 &pos )
{
    if( !IsLocked() ) return;
    BS_ASSERT( m_NumTP < m_MaxTP ); //!\todo Force Sync() and add TP if full
    m_vecTP[m_NumTP++] = ds::TouchedParticle2D( pid, ds::TouchedParticle2D::eTouchedPos, pos );
    Touch(eTouchedPos);
}
void ParticleSys2D::SetVel( int pid, const Vec2 &vel )
{
    if( !IsLocked() ) return;
    BS_ASSERT( m_NumTP < m_MaxTP ); //!\todo Force Sync() and add TP if full
    m_vecTP[m_NumTP++] = ds::TouchedParticle2D( pid, ds::TouchedParticle2D::eTouchedVel, vel );
    Touch(eTouchedVel);
}

void ParticleSys2D::ApplyForce( int pid, const Vec2 &f_global )
{
    if( !IsLocked() ) return;
    BS_ASSERT( m_NumTP < m_MaxTP ); //!\todo Force Sync() and add TP if full
    m_vecTP[m_NumTP++] = ds::TouchedParticle2D( pid, ds::TouchedParticle2D::eTouchedForce, f_global );
    Touch(eTouchedForce);
}
void ParticleSys2D::ApplyImpulse( int pid, const Vec2 &j_global )
{
    if( !IsLocked() ) return;
    BS_ASSERT( m_NumTP < m_MaxTP ); //!\todo Force Sync() and add TP if full
    m_vecTP[m_NumTP++] = ds::TouchedParticle2D( pid, ds::TouchedParticle2D::eTouchedImpulse, j_global );
    Touch(eTouchedImpulse);
}


//---- Internal ISyncEntity methods
bool ParticleSys2D::ProcessUpdate( const ds::ReturnIt &rit )
{
    const Real *p_state = rit.Find("state").GetArrayPtr<Real>();
    memcpy( &m_vecPos[0], p_state, m_Params.m_NumParticles*sizeof(Point2) );
    return true;
}

} // namespace S2
