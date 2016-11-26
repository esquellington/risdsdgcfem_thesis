#ifndef S2_DS_DSH_IENTITY_H_
#define S2_DS_DSH_IENTITY_H_

#include "../Config.h"
#include "../Commands.h" //for EEntityType
#include <util/ItemStream.h>
#include <util/VizStream.h>

namespace S2 { namespace ds {

typedef util::ItemStream ReturnStream;
typedef util::ItemStream::ItemIt ParamIt;

// Fwd declarations
class IDynamicSystemHierarchy;

//! Entity interface. Provides type and parent DSH retrieval.
class IEntity
{
public:
    IEntity( machine_uint_type uid, IDynamicSystemHierarchy* p_parent )
    : m_pParentDSH(p_parent)
    , m_UserID(uid)
    , m_LogFlags(0), m_VizFlags(0), m_ProfFlags(0)
    {}

    //! \name Entity interface
    //@{
    virtual EEntityType GetEntityType() const = 0;
    virtual bool Create( const ParamIt& pit ) = 0;
    virtual bool Edit( const ParamIt& pit ) = 0;
    virtual bool Destroy( ReturnStream& rets ) { return Destroy_Internal(rets,true); }
    virtual bool Internal( const ParamIt& pit, ReturnStream& rets ) { return false; }
    virtual void QueryState( ReturnStream& rets ) const = 0;
    //@}

    //! \name Debug interface
    //@{
    finline void SetLogFlags( Flags32 flags ) { m_LogFlags = flags; }
    finline void SetVizFlags( Flags32 flags ) { m_VizFlags = flags; }
    finline void SetProfFlags( Flags32 flags ) { m_ProfFlags = flags; }
    finline const Flags32& GetLogFlags() const { return m_LogFlags; }
    finline const Flags32& GetVizFlags() const { return m_VizFlags; }
    finline const Flags32& GetProfFlags() const { return m_ProfFlags; }

    //virtual void DoLog( util::LogStream& ls ) const {}//= 0;
    virtual void DoViz( util::VizStream& vs ) const {}//= 0;

    virtual void QueryStats( Flags32 flags, ReturnStream& rets ) const {}//= 0;

    virtual void QueryParams( Flags32 flags, ReturnStream& rets ) const {}//= 0;
    virtual void SyncParams( const ParamIt& pit, ReturnStream& rets ) {}//= 0;
    //@}

    //! \name Bookkeeping
    //@{
    finline IDynamicSystemHierarchy* GetParent() const { return m_pParentDSH; }
    finline machine_uint_type GetUID() const { return m_UserID; }
    finline bool IsDSH() const { return eEntity_DSH == GetEntityType(); }
    finline bool IsGeom() const { return eEntity_Geom == GetEntityType(); }
    //@}

protected:
    //! \name Controlled self-destruction (dangerous?)
    //@{
    virtual ~IEntity() {}
    void Suicide() { delete this; }
    virtual bool Destroy_Internal( ReturnStream& rets, bool b_notify_parent ) = 0;
    //@}

private:
    IDynamicSystemHierarchy* m_pParentDSH;
    machine_uint_type m_UserID;

    //! \name Debug attribs
    //@{
    Flags32 m_LogFlags;
    Flags32 m_VizFlags;
    Flags32 m_ProfFlags;
    //@}
};

} } // namespace S2::ds

#endif // S2_DS_DSH_IENTITY_H_
