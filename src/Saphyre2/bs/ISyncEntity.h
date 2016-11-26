#ifndef S2_BS_ISYNCENTITY_H
#define S2_BS_ISYNCENTITY_H

#include <Saphyre2/bs/BSG.h>

namespace S2 {

/*! S2 Entity base interface */
class IEntity
{
public:
    virtual EEntityType GetType() const = 0;
    virtual const char *GetName() const = 0;
protected:
    IEntity() {}
    virtual ~IEntity() {}
};

//! S2 Query base interfacce (\todo placeholder by now, should add EQueryType, etc...)
class IQuery;

/*! API Generic Syncrhonized Entity.

ISyncEntity life-cycle is expected to be:
0- Allocation { new }
1- Definition { BeginDef(), EndDef() }
2- Creation   { Sync() }
3- Query { None } and Modification { Lock(), Unlock(), Sync() }
4- Update { ProcessReturnIt() }
5- Destruction { Destroy() }
*/
class ISyncEntity: public IEntity
{
private:
    //! ISyncEntity states sorted in expected lifetime order
    enum EState {
        eState_Invalid = 0,
        eState_Allocated,
        eState_BeingDefined,
        eState_WaitingCreation,
        eState_InSync,
        eState_Locked,    //!< In Locked it can be Touched()
        eState_OutSync,   //!< Unlocked but out of sync
        eState_Destroyed, //!< Invalid!!
        eNumStates
    };

public:

    /*! \name Lifetime/Sync API
        \note All prototol calls Return true if success or false if
        error in sync protocol OR in specific call execution
    */
    //@{
    inline bool IsValid() const { return (eState_Invalid != m_State && eState_Destroyed != m_State); }
    inline bool IsBeingDefined() const { return (eState_BeingDefined == m_State); }
    inline bool IsWaitingCreation() const { return (eState_WaitingCreation == m_State); }
    inline bool IsLocked() const { return (eState_Locked == m_State); }
    inline bool IsSync() const { return (eState_InSync == m_State); }

    bool BeginDef();
    bool EndDef( bool b_force_sync = false );
    bool Lock();
    bool Unlock( bool b_force_sync = false );
    bool Sync();
    bool Destroy();
    //@}

    //!\name Bookkeeping API
    //@{
    inline machine_uint_type GetUID() const { return reinterpret_cast<machine_uint_type>(this); }
    inline machine_uint_type GetEID() const { return m_EntityId; }
    inline ISyncEntity *GetParent() const { return m_pParent; }
    inline void SetUserData( void *user_data ) { m_UserData = user_data; }
    inline void *GetUserData() const { return m_UserData; }
    inline void SetName( const char *name ) { m_Name.Set(name); }
    inline const char *GetName() const { return m_Name; }
    //@}

    //! \name Debug API
    //@{
    bool DbgSetFlags_Log( Flags32 flags );
    bool DbgSetFlags_Viz( Flags32 flags );
    bool DbgSetFlags_Prof( Flags32 flags );

    Flags32 DbgGetFlags_Log() const { return m_LogFlags; }
    Flags32 DbgGetFlags_Viz() const { return m_VizFlags; }
    Flags32 DbgGetFlags_Prog() const { return m_ProfFlags; }

    // Synchronous Stats/Params Query/Sync. Requires ISyncEntity::Sync()
    util::ItemStream::ItemItRW Dbg_QueryStats( Flags32 flags = 0xFFFFFFFF );
    util::ItemStream::ItemItRW Dbg_SyncStats();
    util::ItemStream::ItemIt Dbg_GetStats() const; //\note Get current stats, read-only
    util::ItemStream::ItemItRW Dbg_QueryParams( Flags32 flags = 0xFFFFFFFF );
    util::ItemStream::ItemItRW Dbg_SyncParams();

    //! Save raw state
    virtual void Dbg_SaveState( util::ItemStream& its ) const {}
    //@}

    //! A child politely requests to be Synced.
    virtual bool SyncMe( ISyncEntity *p_child ) { return false; }
    //! A query politely requests to be Synced.
    virtual bool SyncMe( IQuery *p_query ) { return false; }

protected:
    friend class BSG;
    friend class Universe;

    ISyncEntity( ISyncEntity *p_parent );
    virtual ~ISyncEntity();

    //! \name ISyncEntity internal protocol
    //@{
    virtual void BeginDef_Internal() = 0;
    virtual bool EndDef_Internal() = 0;
    virtual void Lock_Internal() = 0;
    virtual void Unlock_Internal() = 0;
    virtual bool Sync_Internal();
    virtual bool Destroy_Internal();
    virtual bool ProcessReturn_Internal( const ds::ReturnIt &rit );
    //@}

    //! \name Command sending and result processing
    //@{
    inline void SetEID( machine_uint_type eid ) { m_EntityId = eid; }
    inline void SetChannel( ds::Channel *p_channel ) { m_pChannel = p_channel; }
    inline ds::Channel *GetChannel() const { return m_pChannel; }

    bool ProcessReturn( const ds::ReturnIt &rit );
    virtual bool ProcessUpdate( const ds::ReturnIt &rit ) = 0;
    //@}

    //! \name Touched mask manipulation for subclasses
    //@{
    inline void Touch( int32 tm ) { m_Touched |= tm; }
    inline bool IsTouched( int32 tm = 0xFFFFFFFF ) const { return ( 0 != (m_Touched & tm)); }
    inline int32 GetTouched() const { return m_Touched; }
    inline void Untouch() { m_Touched = 0; }
    //@}

private:
    ISyncEntity *m_pParent;
    EState m_State;
    int32 m_Touched;
    ds::Channel *m_pChannel;
    machine_uint_type m_EntityId;
    Flags32 m_LogFlags, m_VizFlags, m_ProfFlags;

    void *m_UserData;
    String32 m_Name;

#ifdef __S2_BS_ENABLE_DS_STATS
    util::ItemStream *m_pStatsIS;
    Flags32 m_StatsFlags;
#endif
#ifdef __S2_BS_ENABLE_DS_PARAMS
    util::ItemStream *m_pParamsIS;
    Flags32 m_ParamsFlags;
#endif
};

} // namespace S2

#endif // S2_BS_ISYNCENTITY_H
