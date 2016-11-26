#include "IDynamicSystemHierarchy.h"

// Required by IBaseDSH implementation
#include "IGeom.h"
#include "IConnector.h"
#include "IConstraint.h"
#include <Saphyre2/ds/DSG.h>

namespace S2 {
namespace ds {

void IBaseDSH::QueryState( ReturnStream &rets ) const
{
    for( ContainerDSH::iterator it_dsh=m_ChildrenDSH.Begin(); it_dsh.IsValid(); ++it_dsh )
        it_dsh->QueryState(rets);
    
    /*
      for( ContainerConnectors::iterator it_connector=m_Connectors.Begin(); it_connector.IsValid(); ++it_connector )
      it_connector->QueryState(rets);
      
      for( ContainerConstraints::iterator it_constraint=m_Constraints.Begin(); it_constraint.IsValid(); ++it_constraint )
      it_constraint->QueryState(rets);
      
      for( ContainerGeoms::iterator it_geom=m_Geoms.Begin(); it_geom.IsValid(); ++it_geom )
      it_geom->QueryState(rets);
    */
}
    
void IBaseDSH::DoViz( util::VizStream &vs ) const
{
    //todo It's NOT a good idea to offer a flag to filter DoViz
    // recursivity, as it's easy to ignore a lot of viz-enabled
    // children by just missing to enable recursivity on any
    // ancestor. Requiring such fine-grained control from the API
    // is not appropiate by now.    
    //if( GetVizFlags().IsSet(eDbg_DSH_Children) )

    // Viz SS, Children and Geoms
    if( GetVizFlags().Test(eDbg_DSH_SimScheme) && m_pSS )
        m_pSS->DoViz(vs);
    for( ContainerDSH::iterator it_dsh=m_ChildrenDSH.Begin(); it_dsh.IsValid(); ++it_dsh )
        it_dsh->DoViz(vs);
    for( ContainerGeoms::iterator it_geom=m_Geoms.Begin(); it_geom.IsValid(); ++it_geom )
        it_geom->DoViz(vs);
}
//@}

bool IBaseDSH::Destroy_Internal( ReturnStream &rets, bool b_notify_parent )
{
    DS_INFO("IDynamicSystemHierarchy::Destroy_Internal() Entity[" << GetUID() << "] destroyed!");
    
    // Destroy recursively
    for( ContainerDSH::iterator it_dsh=m_ChildrenDSH.Begin(); it_dsh.IsValid(); ++it_dsh )
        it_dsh->Destroy_Internal(rets,false);

    /*
      for( ContainerConnectors::iterator it_connector=m_Connectors.Begin(); it_connector.IsValid(); ++it_connector )
      it_connector->Destroy_Internal(rets,false);

      for( ContainerConstraints::iterator it_constraint=m_Constraints.Begin(); it_constraint.IsValid(); ++it_constraint )
      it_constraint->Destroy_Internal(rets,false);

      for( ContainerGeoms::iterator it_geom=m_Geoms.Begin(); it_geom.IsValid(); ++it_geom )
      it_geom->Destroy_Internal(rets,false);
    */

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
    if( b_notify_parent && GetParent() ) GetParent()->RemoveChild(this);

    // Suicide
    IEntity::Suicide();

    return true;
}

}} // namespace S2::ds
