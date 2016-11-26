#include "ISyncEntity.h"

namespace S2 {

ISyncEntity::ISyncEntity( ISyncEntity *p_parent )
: m_pParent(p_parent)
, m_State(eState_Allocated)
, m_Touched(0)
, m_pChannel(0)
, m_LogFlags(0), m_VizFlags(0), m_ProfFlags(0)
, m_UserData(0)
#ifdef __S2_BS_ENABLE_DS_STATS
, m_pStatsIS(0), m_StatsFlags(0)
#endif
#ifdef __S2_BS_ENABLE_DS_PARAMS
, m_pParamsIS(0), m_ParamsFlags(0)
#endif
{}

ISyncEntity::~ISyncEntity()
{
#ifdef __S2_BS_ENABLE_DS_STATS
    if( 0 != m_pStatsIS ) delete m_pStatsIS;
#endif
#ifdef __S2_BS_ENABLE_DS_PARAMS
    if( 0 != m_pParamsIS ) delete m_pParamsIS;
#endif
}

//---- Lifetime/Sync API
bool ISyncEntity::BeginDef()
{
    switch( m_State )
    {
        // Already being defined (redundant but quiet)
    case eState_BeingDefined:
        return true;

        // Success
    case eState_Allocated:
        m_State = eState_BeingDefined;
        BeginDef_Internal();
        return true;

        // Failure i all other cases
    default:
        return false;
    }
}

bool ISyncEntity::EndDef( bool b_force_sync )
{
    if( !IsBeingDefined() )
        return false;

    if( EndDef_Internal() )
    {
        m_State = eState_WaitingCreation;
        if( b_force_sync )
            Sync();
        return true;
    }
    else
    {
        m_State = eState_Invalid;
        return false;
    }
}

bool ISyncEntity::Lock()
{
    switch( m_State )
    {
        // Already locked (redundant but quiet)
    case eState_Locked:
        return true;

        // Success
    case eState_InSync:
    case eState_OutSync:
        m_State = eState_Locked;
        Lock_Internal();
        return true;

        // Failure in all other cases
    default:
        return false;
    }
}

bool ISyncEntity::Unlock( bool b_force_sync )
{
    if( !IsLocked() )
        return false;

    //if touched, set out of sync
    if( IsTouched() )
    {
        m_State = eState_OutSync;
        Unlock_Internal();
        if( b_force_sync )
            Sync();
    }
    else
        m_State = eState_InSync;

    return true;
}

bool ISyncEntity::Sync()
{
    if( !IsTouched() && !IsWaitingCreation() )
        return true;

    if( Sync_Internal() )
    {
        Untouch();
        m_State = eState_InSync;
        return true;
    }
    else
    {
        m_State = eState_Invalid;
        return false;
    }
}

bool ISyncEntity::Destroy()
{
    if( IsValid() && Destroy_Internal() )
    {
        m_State = eState_Destroyed;
        Sync_Internal();
        return true;
    }
    else
    {
        m_State = eState_Invalid;
        return false;
    }
}

//---- Default Internal methods
bool ISyncEntity::Sync_Internal()
{
    // Sync up-hierarchy if required
    if( GetParent() ) GetParent()->SyncMe( this );
    else BS_LOG_ERROR("ISyncEntity::Sync_Internal() called default implementation with no parent!");
    return true;
}

bool ISyncEntity::Destroy_Internal()
{
    GetChannel()->BeginCommand( ds::eCmd_Destroy );
    {
        GetChannel()->Write( "eid", GetEID() );
    }
    GetChannel()->EndCommand();
    return true;
}

bool ISyncEntity::ProcessReturn_Internal( const ds::ReturnIt &rit )
{
    BS_LOG_ERROR("ISyncEntity::ProcessReturn_Internal() called default implementation!");
    return false;
}

//---- Common return processing
bool ISyncEntity::ProcessReturn( const ds::ReturnIt &rit )
{
    switch( rit.GetType() )
    {
    case ds::eRet_Ok:
        {
            //todo it's a Complex OK!
            return true;
        }
        break;
    case ds::eRet_Error:
        {
            //todo it's a Complex Error! May contain messages!!
            return false;
        }
        break;
    case ds::eRet_NewEntity:
        {
            ds::ReturnIt dit = rit.GetSubItem();
            SetEID( dit.Find("eid").Get<machine_uint_type>() );
            BS_INFO("ISyncEntity::ProcessReturn() Entity[" << GetEID() << "] created!");
            return true;
        }
        break;
    case ds::eRet_KillEntity:
        {
            m_State = eState_Destroyed;
            //\todo DESTROY IT!
            BS_INFO("ISyncEntity::ProcessReturn() Entity[" << GetEID() << "] destroyed!");
            return true;
        }
        break;
    case ds::eRet_Internal: return ProcessReturn_Internal( rit ); break;
    case ds::eRet_Update: return ProcessUpdate( rit.GetSubItem() ); break;
    default: BS_INFO("Unknown Return!"); return false; break;
    }
}

bool ISyncEntity::DbgSetFlags_Log( Flags32 flags )
{
    m_LogFlags = flags;
    GetChannel()->BeginCommand( ds::eDbg_Conf );
    {
        GetChannel()->Write( "eid", GetEID() );
        GetChannel()->Write( "dbg_module", (uint32)ds::eDbgModule_Log );
        GetChannel()->Write( "flags", flags );
    }
    GetChannel()->EndCommand();
    return true;
}

bool ISyncEntity::DbgSetFlags_Viz( Flags32 flags )
{
    m_VizFlags = flags;
    GetChannel()->BeginCommand( ds::eDbg_Conf );
    {
        GetChannel()->Write( "eid", GetEID() );
        GetChannel()->Write( "dbg_module", (uint32)ds::eDbgModule_Viz );
        GetChannel()->Write( "flags", flags );
    }
    GetChannel()->EndCommand();
    return true;
}

bool ISyncEntity::DbgSetFlags_Prof( Flags32 flags )
{
    m_ProfFlags = flags;
    GetChannel()->BeginCommand( ds::eDbg_Conf );
    {
        GetChannel()->Write( "eid", GetEID() );
        GetChannel()->Write( "dbg_module", (uint32)ds::eDbgModule_Prof );
        GetChannel()->Write( "flags", flags );
    }
    GetChannel()->EndCommand();
    return true;
}

util::ItemStream::ItemItRW ISyncEntity::Dbg_QueryStats( Flags32 flags )
{
#ifdef __S2_BS_ENABLE_DS_STATS
    m_StatsFlags = flags;
    if( 0 == m_pStatsIS ) m_pStatsIS = new util::ItemStream(1<<13,1<<11,false); //This IS cannot be realloc, because we store pointers to its contents
    BSG::ClearStats(); //\todo IMPORTANT, to ensure that the new eRet_Stats is the only element in the stream
    Sync();
    GetChannel()->BeginCommand( ds::eDbg_QueryStats );
    {
        GetChannel()->Write( "eid", GetEID() );
        GetChannel()->Write( "flags", flags );
    }
    GetChannel()->EndCommand();
    GetChannel()->FlushCommands();
    // Retrieve first (and only) result from StatsIS and store it into m_pStatsIS
    //\todo Set a proper ID, cannot use GetEID() because it's 64b!
    m_pStatsIS->Clear();
    util::ItemStream::ItemItRW it = m_pStatsIS->WriteItem( 666, BSG::GetStatIS()->Begin() );
    return it.GetSubItem().Find("stats");
#else
    return util::ItemStream::ItemItRW();
#endif
}

util::ItemStream::ItemItRW ISyncEntity::Dbg_SyncStats()
{
#ifdef __S2_BS_ENABLE_DS_STATS
    BSG::ClearStats(); //\todo IMPORTANT, to ensure that the new eRet_Stats is the only element in the stream
    Sync();
    GetChannel()->BeginCommand( ds::eDbg_QueryStats );
    {
        GetChannel()->Write( "eid", GetEID() );
        GetChannel()->Write( "flags", m_StatsFlags );
    }
    GetChannel()->EndCommand();
    GetChannel()->FlushCommands();
    // Retrieve first (and only) result from StatsIS and store it into m_pStatsIS
    util::ItemStream::ItemItRW old_it = m_pStatsIS->BeginRW();
    m_pStatsIS->Clear();
    util::ItemStream::ItemItRW it = m_pStatsIS->WriteItem( 666, BSG::GetStatIS()->Begin() );
    it.TouchRecursive(); //\todo This is required, otherwise PropertyTreeControl does NOT show updated properties...
    BS_ASSERT( it == old_it );
    return it.GetSubItem().Find("stats");
#else
    return util::ItemStream::ItemItRW();
#endif
}

 //\note Get current stats, read-only
util::ItemStream::ItemIt ISyncEntity::Dbg_GetStats() const
{
#ifdef __S2_BS_ENABLE_DS_STATS
    if( m_pStatsIS->IsEmpty() ) return util::ItemStream::ItemIt();
    // Retrieve first (and only) result from StatsIS and store it into m_pStatsIS
    util::ItemStream::ItemIt old_it = m_pStatsIS->Begin();
    return old_it.GetSubItem().Find("stats");
#else
    return util::ItemStream::ItemIt();
#endif
}

util::ItemStream::ItemItRW ISyncEntity::Dbg_QueryParams( Flags32 flags )
{
#ifdef __S2_BS_ENABLE_DS_PARAMS
    m_ParamsFlags = flags;
    if( 0 == m_pParamsIS ) m_pParamsIS = new util::ItemStream(1<<13,1<<11,false); //This IS cannot be realloc, because we store pointers to its contents
    BSG::ClearParams(); //\todo IMPORTANT, to ensure that the new eRet_Params is the only element in the stream
    Sync();
    GetChannel()->BeginCommand( ds::eDbg_QueryParams );
    {
        GetChannel()->Write( "eid", GetEID() );
        GetChannel()->Write( "flags", flags );
    }
    GetChannel()->EndCommand();
    GetChannel()->FlushCommands();
    // Retrieve first (and only) result from ParamsIS and store it into m_pParamsIS
    //\todo Set a proper ID, cannot use GetEID() because it's 64b!
    m_pParamsIS->Clear();
    util::ItemStream::ItemItRW it = m_pParamsIS->WriteItem( 666, BSG::GetParamIS()->Begin() );
    return it.GetSubItem().Find("params");
#else
    return util::ItemStream::ItemItRW();
#endif
}

util::ItemStream::ItemItRW ISyncEntity::Dbg_SyncParams()
{
#ifdef __S2_BS_ENABLE_DS_PARAMS
    BSG::ClearParams(); //\todo IMPORTANT, to ensure that the new eRet_Params is the only element in the stream
    Sync();
    GetChannel()->BeginCommand( ds::eDbg_SyncParams );
    {
        GetChannel()->Write( "eid", GetEID() );
        GetChannel()->WriteItem( "params", m_pParamsIS->Begin().GetSubItem().Find("params") );
    }
    GetChannel()->EndCommand();
    GetChannel()->FlushCommands();
    //\todo ForceSync(); to ensure that all changed params have been propagated to BS cached state
    // Retrieve first (and only) result from ParamsIS and store it into m_pParamsIS
    util::ItemStream::ItemItRW old_it = m_pParamsIS->BeginRW();
    // \todo Update old_it with new values and Touch anything that has changed
    m_pParamsIS->Clear();
    util::ItemStream::ItemItRW it = m_pParamsIS->WriteItem( 666, BSG::GetParamIS()->Begin() );
    BS_ASSERT( it == old_it ); //TEMPORAL: It's redundant, but for safety we'll leave it
    return it.GetSubItem().Find("params");
#else
    return util::ItemStream::ItemItRW();
#endif
}

} //namespace S2
