#include <Saphyre2/bs/Universe.h>

#include <Saphyre2/bs/Particle3D.h>
#include <Saphyre2/bs/Particle2D.h>
#include <Saphyre2/bs/ParticleSys2D.h>
#include <Saphyre2/bs/Fluid2D.h>
#include <Saphyre2/bs/Solid2D.h>
#include <Saphyre2/bs/Solid3D.h>
#include <Saphyre2/bs/Kine2D.h>
#include <Saphyre2/bs/Kine3D.h>

#include <Saphyre2/bs/RayCastQuery.h>

#include <Saphyre2/bs/BSG.h>
#include <Saphyre2/ds/Commands.h>

namespace S2 {

Universe::Universe()
: ISyncEntity(0)
, m_LastTime(0)
{
}

Universe::~Universe()
{

}

Particle3D *Universe::CreateParticle3D()
{
    BS_BCMD("Universe::CreateParticle3D()");

    bool bResult(true);
    if( IsWaitingCreation() )
    {
        BS_INFO("Calling Universe::Sync() automatically...");
        bResult = Sync();
    }

    Particle3D *pP3D = 0;
    if( bResult )
    {
        pP3D = new Particle3D(this);
        pP3D->SetChannel(GetChannel());
        m_Children.Add( pP3D );
    }

    BS_ECMD(bResult);
    return pP3D;
}

Particle2D *Universe::CreateParticle2D()
{
    BS_BCMD("Universe::CreateParticle2D()");

    bool bResult(true);
    if( IsWaitingCreation() )
    {
        BS_INFO("Calling Universe::Sync() automatically...");
        bResult = Sync();
    }

    Particle2D *pP2D = 0;
    if( bResult )
    {
        pP2D = new Particle2D(this);
        pP2D->SetChannel(GetChannel());
        m_Children.Add(pP2D);
    }

    BS_ECMD(bResult);
    return pP2D;
}

ParticleSys2D *Universe::CreateParticleSys2D()
{
    BS_BCMD("Universe::CreateParticleSys2D()");

    bool bResult(true);
    if( IsWaitingCreation() )
    {
        BS_INFO("Calling Universe::Sync() automatically...");
        bResult = Sync();
    }

    ParticleSys2D *pPSys2D = 0;
    if( bResult )
    {
        pPSys2D = new ParticleSys2D(this);
        pPSys2D->SetChannel(GetChannel());
        m_Children.Add(pPSys2D);
    }

    BS_ECMD(bResult);
    return pPSys2D;
}

Fluid2D *Universe::CreateFluid2D()
{
    BS_BCMD("Universe::CreateFluid2D()");

    bool bResult(true);
    if( IsWaitingCreation() )
    {
        BS_INFO("Calling Universe::Sync() automatically...");
        bResult = Sync();
    }

    Fluid2D *pF2D = 0;
    if( bResult )
    {
        pF2D = new Fluid2D(this);
        pF2D->SetChannel(GetChannel());
        m_Children.Add(pF2D);
    }

    BS_ECMD(bResult);
    return pF2D;
}

Solid2D *Universe::CreateSolid2D()
{
    BS_BCMD("Universe::CreateSolid2D()");

    bool bResult(true);
    if( IsWaitingCreation() )
    {
        BS_INFO("Calling Universe::Sync() automatically...");
        bResult = Sync();
    }

    Solid2D *pS2D = 0;
    if( bResult )
    {
        pS2D = new Solid2D(this);
        pS2D->SetChannel(GetChannel());
        m_Children.Add(pS2D);
    }

    BS_ECMD(bResult);
    return pS2D;
}

Solid3D *Universe::CreateSolid3D()
{
    BS_BCMD("Universe::CreateSolid3D()");

    bool bResult(true);
    if( IsWaitingCreation() )
    {
        BS_INFO("Calling Universe::Sync() automatically...");
        bResult = Sync();
    }

    Solid3D *pS3D = 0;
    if( bResult )
    {
        pS3D = new Solid3D(this);
        pS3D->SetChannel(GetChannel());
        m_Children.Add(pS3D);
    }

    BS_ECMD(bResult);
    return pS3D;
}

Kine2D *Universe::CreateKine2D()
{
    BS_BCMD("Universe::CreateKine2D()");

    bool bResult(true);
    if( IsWaitingCreation() )
    {
        BS_INFO("Calling Universe::Sync() automatically...");
        bResult = Sync();
    }

    Kine2D *pK2D = 0;
    if( bResult )
    {
        pK2D = new Kine2D(this);
        pK2D->SetChannel(GetChannel());
        m_Children.Add(pK2D);
    }

    BS_ECMD(bResult);
    return pK2D;
}


Kine3D *Universe::CreateKine3D()
{
    BS_BCMD("Universe::CreateKine3D()");

    bool bResult(true);
    if( IsWaitingCreation() )
    {
        BS_INFO("Calling Universe::Sync() automatically...");
        bResult = Sync();
    }

    Kine3D *pK3D = 0;
    if( bResult )
    {
        pK3D = new Kine3D(this);
        pK3D->SetChannel(GetChannel());
        m_Children.Add(pK3D);
    }

    BS_ECMD(bResult);
    return pK3D;
}

RayCastQuery2D *Universe::CreateRayCastQuery2D()
{
    BS_BCMD("Universe::CreateRayCastQuery2D()");

    bool bResult(true);
    if( IsWaitingCreation() )
    {
        BS_INFO("Calling Universe::Sync() automatically...");
        bResult = Sync();
    }

    RayCastQuery2D *pRCQ(0);
    if( bResult )
    {
        pRCQ = new RayCastQuery2D( this, GetChannel() );
        m_Queries.Add(pRCQ);
    }

    BS_ECMD(bResult);
    return pRCQ;
}

RayCastQuery3D *Universe::CreateRayCastQuery3D()
{
    BS_BCMD("Universe::CreateRayCastQuery3D()");

    bool bResult(true);
    if( IsWaitingCreation() )
    {
        BS_INFO("Calling Universe::Sync() automatically...");
        bResult = Sync();
    }

    RayCastQuery3D *pRCQ(0);
    if( bResult )
    {
        pRCQ = new RayCastQuery3D( this, GetChannel() );
        m_Queries.Add(pRCQ);
    }

    BS_ECMD(bResult);
    return pRCQ;
}

bool Universe::Update( Real dt )
{
    if( !IsValid() ) return false;

    BS_BCMD("Universe::Update()");

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

bool Universe::ForceSync()
{
    if( !IsValid() ) return false;
    BS_BCMD("Universe::ForceSync()");
    bool bResult = Sync_Internal();
    BS_ECMD( bResult );
    return bResult;
}

//---- ISyncEntity internal protocol
//! \pre IsBeingDefined()
bool Universe::EndDef_Internal()
{
    BS_BCMD("Universe::EndDef()");

    //---- Send creation command
    GetChannel()->BeginCommand( ds::eCmd_Create );
    {
        GetChannel()->Write( "uid", GetUID() );
        GetChannel()->Write( "peid", (machine_uint_type)0 );
        GetChannel()->Write( "entity_type", (uint32)ds::eEntity_DSH );
        GetChannel()->Write( "dsh_type", (uint32)ds::eDSH_Universe3D );
        GetChannel()->Write( "flags", m_Params.m_Flags );
        GetChannel()->Write( "length_scale", m_Params.m_LengthScale );
        GetChannel()->Write( "mass_scale", m_Params.m_MassScale );
        GetChannel()->Write( "time_scale", m_Params.m_TimeScale );
    }
    GetChannel()->EndCommand();

    bool bResult(true);
    BS_ECMD( bResult );
    return bResult;
}

//! \pre IsTouched() or IsBeingCreated() or called from SyncMe(child) or from Universe::Update(dt)
bool Universe::Sync_Internal()
{
    BS_BCMD("Universe::Sync()");
    // Flush
    bool bResult = GetChannel()->FlushCommands();
    // Process returns directly or dispatch to children
    if( bResult )
    {
        for( ds::ReturnIt rit=GetChannel()->GetReturnIt();
             bResult && rit.IsValid();
             ++rit )
        {
            if( rit.IsSimple() )
            {
                switch( rit.GetType() )
                {
                case ds::eRet_Ok: bResult = true; break;
                case ds::eRet_Error: bResult = false; break;
                default: bResult = false; break;
                }
            }
            else // IsComplex
            {
                switch( rit.GetType() )
                {
                case ds::eRet_RayCast: //\todo ds::eRet_Query or process all as IQuery using fall-through switch-case with ds::eRet_Intersection, etc...??
                    //\todo case ds::eRet_Intersection:
                    {
                        IQuery *pQuery = reinterpret_cast<IQuery*>( rit.GetSubItem().Find("qid").Get<machine_uint_type>() );
                        //\todo CHECK IF VALID IQuery!!
                        bResult = bResult && pQuery->ProcessReturn( rit );
                    }
                    break;
                default:
                    {
                        ISyncEntity *pSE = reinterpret_cast<ISyncEntity*>( rit.GetSubItem().Find("uid").Get<machine_uint_type>() );
                        //\todo CHECK IF VALID ISyncEntity!!
                        bResult = bResult && pSE->ProcessReturn( rit );
                    }
                    break;
                }
            }
        }
    }
    // Discard processed returns
    GetChannel()->DiscardReturn();
    BS_ECMD(bResult);
    return bResult;
}

bool Universe::ProcessUpdate( const ds::ReturnIt &rit )
{
    return true;
}

/*! Even if the Universe is not out-of-sync, when a child requests so,
  we perform a whole Universe Sync() unconditionally. Most often, only
  child's commands will be sent, but occasionally other commands will
  be sent preserving chronological order at the Universe-level.
*/
bool Universe::SyncMe( ISyncEntity *p_child )
{
    return Sync_Internal();
}

/*! Even if the Universe is not out-of-sync, when a query requests so,
  we perform a whole Universe Sync() unconditionally. Most often, only
  query's commands will be sent, but occasionally other commands will
  be sent preserving chronological order at the Universe-level.
*/
bool Universe::SyncMe( IQuery *p_query )
{
    return Sync_Internal();
}

} // namespace S2
