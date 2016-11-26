#include "IConstraint.h"
#include "IDynamicSystemHierarchy.h"
#include <Saphyre2/ds/DSG.h>

namespace S2 { namespace ds {

bool IBaseConstraint::Destroy_Internal( ReturnStream &rets, bool b_notify_parent )
{
    // Notify destruction if it's a user-known entity
    if( GetUID() > 0 )
    {
        rets.BeginComplex( 666, (uint32)eRet_KillEntity ); //\todo cid = 666??!!
        {
            rets.Write("uid",GetUID());
        }
        rets.EndComplex();
    }        
    // Tell parent
    if( b_notify_parent && GetParent() ) GetParent()->RemoveConstraint(this);
    // Suicide
    IEntity::Suicide();
    return true;
}

}} //namespace S2::ds
