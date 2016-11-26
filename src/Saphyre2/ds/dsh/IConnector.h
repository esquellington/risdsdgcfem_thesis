#ifndef S2_DS_DSH_ICONNECTOR_H
#define S2_DS_DSH_ICONNECTOR_H

#include "IEntity.h"
#include <Saphyre2/ms/IConnector.h>

namespace S2 { namespace ds {

// Fwd declarations
class IDynamicSystemHierarchy;

/*! Connector Entity.
  
  Wraps an IConnector implementation and adds external access through
  standard IEntity interface.
*/
class IConnector: public IEntity
{
public:
    finline IConnector( uint32 uid, IDynamicSystemHierarchy *parent ) : IEntity(uid,parent) {}

    //!\name IConnector interface
    //@{
    virtual ms::IConnector *GetImpl() = 0;
    virtual const ms::IConnector *GetImpl() const = 0;
    //@}

    //\todo Implement IEntity interface: Create,Edit,Destroy...
};

class IBaseConnector: public IConnector
{
public:
    finline IBaseConnector( uint32 uid, IDynamicSystemHierarchy *parent ) : IConnector(uid,parent) , m_pImpl(0) {}

    //!\name IConnector implementation
    //@{
    ms::IConnector *GetImpl() { return m_pImpl; }
    const ms::IConnector *GetImpl() const { return m_pImpl; }
    //@}
    
    //!\name IEntity implementation
    //@{
    EEntityType GetEntityType() const { return eEntity_Connector; }
    bool Destroy_Internal( ReturnStream &rets, bool b_notify_parent );
    //@}
    
protected:
    ms::IConnector *m_pImpl;
};
    
}} // namespace S2::ds

#endif // S2_DS_DSH_ICONNECTOR_H
