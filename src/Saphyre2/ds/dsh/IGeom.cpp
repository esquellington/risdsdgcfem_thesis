#include "IGeom.h"
#include "IDynamicSystemHierarchy.h"
#include <Geo/util/Viz.h>

namespace S2 {
namespace ds {

//---- IBaseGeom implementation
bool IBaseGeom::Destroy_Internal( ReturnStream &rets, bool b_notify_parent )
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
    if( b_notify_parent && GetParent() ) GetParent()->RemoveGeom(this);
    // Suicide
    IEntity::Suicide();
    return true;
}

void IBaseGeom::DoViz( util::VizStream &vs ) const
{
    if( GetVizFlags().Test(eDbg_State) )
        geo::VizObject(GetGO(), vs, geo::eODDF_Shape | geo::eSDDF_DCR );//| geo::eSDDF_Default ); // ); //DCR is EXTREMELY SLOW, specially if __ENABLE_DCR3_DRAW_ERRORS is defined in Geo/Viz.cpp
    // if( GetVizFlags().Test(eDbg_Geom_BV) && 0 != GetBV() )
    //     geo::VizBoundingVolume(GetBV(),vs);
}

}} // namespace S2::ds
