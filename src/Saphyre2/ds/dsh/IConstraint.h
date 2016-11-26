#ifndef S2_DS_DSH_ICONSTRAINT_H
#define S2_DS_DSH_ICONSTRAINT_H

#include "IEntity.h"
#include <Saphyre2/ms/IConstraint.h>

namespace S2 { namespace ds {

// Fwd declarations
class IDynamicSystemHierarchy;

/*! Constraint Entity.
  
  Wraps an IConstraint implementation and adds external access through
  standard IEntity interface.
*/
class IConstraint: public IEntity
{
public:
    finline IConstraint( uint32 uid, IDynamicSystemHierarchy *parent ) : IEntity(uid,parent) {}

    virtual ms::IConstraint *GetImpl() = 0;
    virtual const ms::IConstraint *GetImpl() const = 0;
};
    
class IBaseConstraint: public IConstraint
{
public:
    finline IBaseConstraint( uint32 uid, IDynamicSystemHierarchy *parent ) : IConstraint(uid,parent), m_pImpl(0) {}

    //!\name IConstraint implementation
    //@{
    ms::IConstraint *GetImpl() { return m_pImpl; }
    const ms::IConstraint *GetImpl() const { return m_pImpl; }
    //@}
    
    //!\name IEntity implementation
    //@{
    EEntityType GetEntityType() const { return eEntity_Constraint; }
    bool Destroy_Internal( ReturnStream &rets, bool b_notify_parent );
    //@}
    
protected:
    ms::IConstraint *m_pImpl;
};    

}} // namespace S2::ds

#endif // S2_DS_DSH_ICONSTRAINT_H
